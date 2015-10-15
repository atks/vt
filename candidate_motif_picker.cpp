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

    if (variant.ins)
    {
        char** alleles = bcf_get_allele(v);

        if (debug)
        {
            const char* repeat_tract = variant.vntr.exact_repeat_tract.c_str();
            std::cerr << "Longest Allele : "   << alleles[0][0] << "[" <<  &alleles[1][1]  << "]" << &repeat_tract[1] << "\n";
        }

        //spike in inserted allele
        std::string spiked_seq(alleles[1]);
        std::string insertion = variant.vntr.exact_repeat_tract.substr(strlen(alleles[0]), variant.vntr.exact_repeat_tract.size()-strlen(alleles[0]));
        spiked_seq.append(insertion);
        mt->detect_candidate_motifs(spiked_seq);
        
        indel_sequence.assign(&alleles[1][1]);
    }
    else
    {
        mt->detect_candidate_motifs(variant.vntr.exact_repeat_tract);

        char** alleles = bcf_get_allele(v);
        indel_sequence.assign(&alleles[0][1]);
    }
    
    if (debug)
    {
        std::cerr << "Indel : "  << indel_sequence << "\n";
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
        
        if (indel_sequence.find(shifted_motif.c_str())!=std::string::npos)
        {
            return true;
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

