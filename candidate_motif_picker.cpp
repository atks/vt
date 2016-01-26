/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include "candidate_motif_picker.h"

/**
 * Constructor.
 */
CandidateMotifPicker::CandidateMotifPicker(bool debug)
{
    max_mlen = 10;
    mt = new MotifTree(max_mlen, debug);

    this->debug = debug;
};

/**
 * Destructor.
 */
CandidateMotifPicker::~CandidateMotifPicker()
{
    delete mt;
}

/**
 * Pick candidate motifs in different modes.
 * Invokes motif tree and the candidate motifs are stored in a
 * heap within the motif tree.
 */
void CandidateMotifPicker::generate_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, Variant& variant)
{
    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "PICK CANDIDATE MOTIFS\n\n";
    }

    this->v = v;

    bool alt_is_longest_allele = false;
    char** alleles = bcf_get_allele(v);
    uint32_t n_allele = bcf_get_n_allele(v);
    uint32_t longest_allele_index = 0;
    uint32_t longest_allele_length = strlen(alleles[0]);

    for (uint32_t i=1; i<n_allele; ++i)
    {
        uint32_t len = strlen(alleles[i]);
        if (len>longest_allele_length)
        {
            longest_allele_index = i;
            longest_allele_length = len;
        }
    }

    if (debug)
    {
        std::cerr << "Alleles              : ";
        for (uint32_t i=0; i<n_allele; ++i)
        {
            if (i) std::cerr << ",";
            std::cerr << alleles[i];
        }
        std::cerr << "\n";

        std::cerr << "Repeat Tract Position : [" << variant.vntr.exact_rbeg1 << "," << variant.vntr.exact_rend1 << "]\n";
        std::cerr << "Repeat Tract          : " << variant.vntr.exact_repeat_tract << "\n";
        std::cerr << "Longest Allele        : " << alleles[longest_allele_index] << "\n";
        std::cerr << "Longest Allele index  : " << longest_allele_index << "\n";
            
        if (longest_allele_index)
        {
            std::string spiked_seq = variant.vntr.exact_repeat_tract;
            std::cerr << "\tspos1          " << variant.vntr.exact_rbeg1-bcf_get_pos1(v) << "\n";
            std::cerr << "\treplace length " << strlen(alleles[0]) << "\n";
                   
            spiked_seq.replace(variant.vntr.exact_rbeg1-bcf_get_pos1(v), strlen(alleles[0]), alleles[longest_allele_index]);
            spiked_seq.insert(variant.vntr.exact_rbeg1-bcf_get_pos1(v), 1, '[');
            spiked_seq.insert(variant.vntr.exact_rbeg1-bcf_get_pos1(v)+ strlen(alleles[longest_allele_index])+1, 1, ']');
            std::cerr << "Spiked Longest Allele : "   << spiked_seq << "\n";
        }   
    }

    if (longest_allele_index)
    {
        std::string spiked_seq = variant.vntr.exact_repeat_tract;
        spiked_seq.replace(variant.vntr.exact_rbeg1-bcf_get_pos1(v), strlen(alleles[0]), alleles[longest_allele_index]);
        mt->detect_candidate_motifs(spiked_seq);
    }
    else
    {
        mt->detect_candidate_motifs(variant.vntr.exact_repeat_tract);
    }
}

/**
 * Initialize candidate motif from VCF record.
 */
void CandidateMotifPicker::set_motif_from_info_field(Variant& variant)
{
    VNTR& vntr = variant.vntr;
    char *motif = NULL;
    int32_t n = 0;
    if (bcf_get_info_string(variant.h, variant.v, "MOTIF", &motif, &n)>0)
    {
        vntr.motif.assign(motif);
        vntr.mlen = vntr.motif.size();
        free(motif);

        vntr.basis = vntr.get_basis(vntr.motif);
    }
    else
    {
        vntr.motif = "";
    }   
}
    
/**
 * Choose the next best motif.
 */
bool CandidateMotifPicker::next_motif(bcf_hdr_t* h, bcf1_t* v, Variant& variant)
{
    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "PICKING NEXT BEST MOTIF\n\n";
    }

    while (!mt->pcm.empty())
    {
        CandidateMotif cm = mt->pcm.top();

        //check for existence of pattern in indel sequence
        if (is_in_indel_fragment(cm.motif))
        {
            if (debug)
            {
                printf("selected: %10s %.2f %.2f\n", mt->pcm.top().motif.c_str(),
                                                                mt->pcm.top().score,
                                                                mt->pcm.top().fit);
            }
            variant.vntr.motif = cm.motif;
            variant.vntr.mlen = cm.motif.size();
            variant.vntr.motif_score = cm.score;
            mt->pcm.pop();
            return true;
        }
        else
        {
            if (debug)
            {
                printf("rejected: %10s %.2f %.2f (not in indel fragment)\n",
                                                                mt->pcm.top().motif.c_str(),
                                                                mt->pcm.top().score,
                                                                mt->pcm.top().fit);
            }
            mt->pcm.pop();
        }
    }

    return false;
}

/**
 * Checks if motif is in indel fragment.
 */
bool CandidateMotifPicker::is_in_indel_fragment(std::string motif)
{
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string shifted_motif = motif.substr(i) + motif.substr(0,i);

        char** alleles = bcf_get_allele(v);
        uint32_t n_allele = bcf_get_n_allele(v);
        for (uint32_t i=0; i<n_allele; ++i)
        {
            if (strstr(alleles[i], shifted_motif.c_str()))
            {
                return true;
            }
        }
    }

    return false;
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
std::string CandidateMotifPicker::choose_repeat_unit(std::string& ref, std::string& motif)
{
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = mt->shift_str(motif, i);
        if (ref.compare(0, smotif.size(), smotif)==0)
        {
            return smotif;
        }
    }

    return motif;
}