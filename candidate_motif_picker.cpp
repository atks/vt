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
 *
 * 1. gets reference sequence from the REF column of a VCF record
 * 2. examines this sequence for candidate motifs that are stored in the priority queue in the MotifTree class.
 * 3. the candidate motifs are then to be accessed via next_motif()
 */
void CandidateMotifPicker::generate_candidate_motifs(Variant& variant)
{
    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "PICK CANDIDATE MOTIFS\n\n";
    }

    this->v = variant.v;

    bool alt_is_longest_allele = false;
    char** alleles = bcf_get_allele(v);
    uint32_t pos1 = bcf_get_pos1(v);
    uint32_t n_allele = bcf_get_n_allele(v);
    uint32_t longest_allele_index = 0;
    uint32_t ref_len = strlen(alleles[0]);
    uint32_t longest_allele_length = ref_len;

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

        std::cerr << "Repeat Tract Position : [" << variant.vntr.exact_beg1 << "," << variant.vntr.exact_end1 << "]\n";
        std::cerr << "Repeat Tract          : " << variant.vntr.exact_repeat_tract << "\n";
        std::cerr << "Longest Allele        : " << alleles[longest_allele_index] << "\n";
        std::cerr << "Longest Allele index  : " << longest_allele_index << "\n";
    }

    if (longest_allele_index)
    {
        //in the case of multiallelics, because the repeat tract is actually obtain by merging a regions from
        //pairwise left and right alignment of the alternative alleles with the reference allele, it is possible that the
        //tract occurs prior to the position of the multiallelic variant.
        int32_t offset = 0;
        if (variant.vntr.exact_beg1<pos1)
        {
            offset = bcf_get_pos1(v) - variant.vntr.exact_beg1;
        }

        std::string spiked_seq = variant.vntr.exact_repeat_tract;
        spiked_seq.replace(offset, ref_len, alleles[longest_allele_index]);
        mt->detect_candidate_motifs(spiked_seq);

        if (debug)
        {
            //this implictly requires that the variants are left aligned.
            std::string spiked_seq = variant.vntr.exact_repeat_tract;
            std::cerr << "\texact repeat tract " << variant.vntr.exact_repeat_tract << "\n";
            std::cerr << "\trbeg1              " << variant.vntr.exact_beg1 << "\n";
            std::cerr << "\tpos1               " << bcf_get_pos1(v) << "\n";
            std::cerr << "\toffset             " << offset << "\n";
            std::cerr << "\treplace length     " << ref_len << "\n";

//            spiked_seq.replace(variant.vntr.exact_beg1-bcf_get_pos1(v), strlen(alleles[0]), alleles[longest_allele_index]);
//            spiked_seq.insert(variant.vntr.exact_beg1-bcf_get_pos1(v), 1, '[');
//            spiked_seq.insert(variant.vntr.exact_beg1-bcf_get_pos1(v)+ strlen(alleles[longest_allele_index])+1, 1, ']');
//            std::cerr << "Spiked Longest Allele : "   << spiked_seq << "\n";
        }

    }
    else
    {
        mt->detect_candidate_motifs(variant.vntr.exact_repeat_tract);
    }
}

/**
 * Generate candidate motifs from a repeat tract.
 *
 * 1. assigns the repeat_tract to the fuzzy_repeat_tract field of the variants VNTR.
 * 2. examines this sequence for candidate motifs that are stored in the priority queue in the MotifTree class.
 * 3. the candidate motifs are then to be accessed via next_motif()
 */
void CandidateMotifPicker::generate_candidate_motifs(char* repeat_tract, Variant& variant)
{
//    std::cerr << "INVOKED??\n" ;
    variant.vntr.fuzzy_repeat_tract.assign(repeat_tract);
//    std::cerr << "assigned repeat tract"  << repeat_tract << "\n";

    mt->detect_candidate_motifs(variant.vntr.fuzzy_repeat_tract);
//    std::cerr << mt->pcm.size() << "\n";

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
        vntr.basis = vntr.get_basis(vntr.motif);
        vntr.mlen = vntr.motif.size();
        vntr.blen = vntr.basis.size();
        free(motif);
    }
    else
    {
        vntr.motif = "";
    }
}

/**
 * Gets inserted or deleted allele of a biallelic indel.
 *
 * Returns true if indel allele is simple and false if ambiguous
 */
bool CandidateMotifPicker::get_indel(std::string ref, std::string alt, std::string& indel)
{
    //trim right
    while (ref.size()!=0 && alt.size()!=0 && ref.at(ref.size()-1)==alt.at(alt.size()-1))
    {
        ref.erase(ref.size()-1, 1);
        alt.erase(alt.size()-1, 1);
    }

    //trim left
    while (ref.size()!=0 && alt.size()!=0 && ref.at(0)==alt.at(0))
    {
        ref.erase(0, 1);
        alt.erase(0, 1);
    }

    if (ref.size()==0)
    {
        indel.assign(alt);
        return true;
    }

    if (alt.size()==0)
    {
        indel.assign(ref);
        return true;
    }

    //not a simple indel
    if (ref.size()>alt.size())
    {
        indel.assign(ref.substr(alt.size(),ref.size()-alt.size()));
    }
    else if (ref.size()<alt.size())
    {
        indel.assign(alt.substr(ref.size(),alt.size()-ref.size()));
    }
    else
    {
        fprintf(stderr, "[%s:%d %s] Not an indel!!!!  ref:%s alt:%s\n", __FILE__, __LINE__, __FUNCTION__, ref.c_str(), alt.c_str());
        exit(1);
    }

    return false;
}

/**
 * Updates the motif of an indel allele.
 *
 * Returns true if the motif is from a simple indel.
 * Returns false if the motif is ambiguous.
 */
void CandidateMotifPicker::update_exact_repeat_unit(Variant& variant)
{
    VNTR& vntr = variant.vntr;
    bcf1_t* v = variant.v;
    char** alleles = bcf_get_allele(v);
    int32_t n_allele = bcf_get_n_allele(v);
    std::string ref;
    std::string alt;
    std::string indel;
    std::string indel_repeat_unit;

    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "UPDATE EXACT REPEAT UNIT\n\n";
    }

    if (n_allele==2)
    {
        ref.assign(alleles[0]);
        alt.assign(alleles[1]);

        bool is_simple_indel = get_indel(ref, alt, indel);

        if (is_simple_indel)
        {
            int32_t i = VNTR::is_periodic(indel);
            if (i)
            {
                indel_repeat_unit.assign(indel.substr(0, i));
            }
            else
            {
                indel_repeat_unit.assign(indel);
            }

            if (debug) std::cerr << "simple indel: " << ref << "/" << alt << " => " << indel_repeat_unit << "\n";
        }
        else
        {
            indel_repeat_unit.assign(indel);

            if (debug) std::cerr << "NON simple indel: " << ref << "/" << alt << " => "  << indel_repeat_unit << "\n";
        }

        vntr.exact_ru_ambiguous = !is_simple_indel;
    }
    else
    {
        if (debug)
        {
            std::cerr << "multiallelic: examine all alleles\n";
        }

        bool all_are_simple_indels = true;
        std::map<std::string, int32_t> indels;

        for (int32_t i=1; i<n_allele; ++i)
        {
            ref.assign(alleles[0]);
            alt.assign(alleles[i]);

            if (debug) std::cerr << "allele: "  << i << "\n";


            if (ref.size()==alt.size())
            {
                if (debug) std::cerr << "\tnon indel: "  << ref << "/" << alt << "\n";
                continue;
            }

            bool is_simple_indel =  get_indel(ref, alt, indel);
            all_are_simple_indels = all_are_simple_indels && is_simple_indel;

            if (is_simple_indel)
            {
                int32_t i = VNTR::is_periodic(indel);
                if (i)
                {
                    indel_repeat_unit.assign(indel.substr(0, i));
                }
                else
                {
                    indel_repeat_unit.assign(indel);
                }

                if (debug) std::cerr << "\tsimple indel: " << ref << "/" << alt << " => " << indel_repeat_unit << "\n";
            }
            else
            {
                indel_repeat_unit.assign(indel);
                if (debug) std::cerr << "\tNON simple indel: " << ref << "/" << alt << " => " << indel_repeat_unit << "\n";
            }

            ++indels[indel_repeat_unit];
        }

        if (debug) std::cerr << "\narbitrating ...\n";

        std::map<std::string, int32_t>::iterator i = indels.begin();
        std::string best_repeat_unit = "";
        int32_t best_count = 0;
        while (i!=indels.end())
        {
             if (debug) std::cerr << "\t" << i->first << " : " << i->second << "\n";
            
            if (i->second>best_count)
            {
                best_repeat_unit = i->first;
                best_count = i->second;
            }

            ++i;
        }

        if (debug) std::cerr << "select : " <<  best_repeat_unit << " (" << best_count << ")\n";

        vntr.exact_ru_ambiguous = !all_are_simple_indels;
        indel_repeat_unit = best_repeat_unit;
    }

    vntr.exact_ru = indel_repeat_unit;
    vntr.exact_motif = VNTR::canonicalize2(indel_repeat_unit);
    vntr.exact_basis = VNTR::get_basis(vntr.exact_motif);
    vntr.exact_mlen = vntr.exact_motif.size();
    vntr.exact_blen = vntr.exact_basis.size();
}

/**
 * Iterates through the candidate motifs detected in the motif tree.
 *
 *  1. examines it if the motif is represented in the motif tree.
 *     and updated the motif field in the VNTR object of variant
 *     with the candidate motif.
 *  2.
 *
 */
bool CandidateMotifPicker::next_motif(Variant& variant, int32_t mode)
{
    if (mode==CHECK_MOTIF_PRESENCE_IN_ALLELE)
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
                
                VNTR& vntr = variant.vntr;
                vntr.fuzzy_motif = VNTR::canonicalize2(cm.motif);
                vntr.fuzzy_basis = VNTR::get_basis(cm.motif);
                vntr.fuzzy_mlen = vntr.fuzzy_motif.size();
                vntr.fuzzy_blen = vntr.fuzzy_basis.size();
  
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
    else if (mode==NO_REQUIREMENT)
    {
        if (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();
            variant.vntr.motif = cm.motif;
            variant.vntr.basis = VNTR::get_basis(cm.motif);
            variant.vntr.mlen = cm.motif.size();
            variant.vntr.blen = variant.vntr.basis.size();
            variant.vntr.fuzzy_motif = variant.vntr.motif;
            variant.vntr.fuzzy_basis = variant.vntr.basis;
            variant.vntr.fuzzy_mlen = variant.vntr.mlen;
            variant.vntr.fuzzy_blen = variant.vntr.blen;
            mt->pcm.pop();

//            std::cerr << variant.vntr.motif << " " << variant.vntr.motif_score  << "\n";

            return true;
        }

        return false;
    }
    else
    {
        fprintf(stderr, "[E:%s:%d %s] Motif picking mode not recognized : %d\n", __FILE__, __LINE__, __FUNCTION__, mode);
        exit(1);
        return false;
    }
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