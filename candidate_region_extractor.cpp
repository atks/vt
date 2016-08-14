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

#include "candidate_region_extractor.h"

/**
 * Constructor.
 */
CandidateRegionExtractor::CandidateRegionExtractor(std::string& ref_fasta_file, bool debug)
{
    vm = new VariantManip(ref_fasta_file.c_str());
    rs = new ReferenceSequence(ref_fasta_file);

    this->debug = debug; 
};

/**
 * Destructor.
 */
CandidateRegionExtractor::~CandidateRegionExtractor()
{
    delete vm;
    delete rs;
}

/**
 * Pick candidate region.
 *
 * @mode - REFERENCE     use refence field
 *       - ALLELE_EXACT  by exact alignment
 *       - ALLELE_FUZZY  by fuzzy alignment
 */
void CandidateRegionExtractor::pick_candidate_region(Variant& variant, int32_t mode, int32_t amode)
{
    bcf_hdr_t* h = variant.h;
    bcf1_t* v = variant.v;

    if (mode==REFERENCE)
    {
        if (amode&FINAL)
        {
            VNTR& vntr = variant.vntr;
            vntr.repeat_tract.assign(bcf_get_ref(v));
            vntr.beg1 = bcf_get_pos1(v);
            vntr.end1 = bcf_get_end1(v);
            vntr.rl = vntr.end1-vntr.beg1+1;
            vntr.ll = vntr.rl; //??
        }

        if (amode&EXACT)
        {
            VNTR& vntr = variant.vntr;
            vntr.exact_repeat_tract.assign(bcf_get_ref(v));
            vntr.exact_beg1 = bcf_get_pos1(v);
            vntr.exact_end1 = bcf_get_end1(v);
            vntr.exact_rl = vntr.exact_end1-vntr.exact_beg1+1;
            vntr.exact_ll = vntr.exact_rl; //??
        }

        if (amode&FUZZY)
        {
            VNTR& vntr = variant.vntr;
            vntr.fuzzy_repeat_tract.assign(bcf_get_ref(v));
            vntr.fuzzy_beg1 = bcf_get_pos1(v);
            vntr.fuzzy_end1 = bcf_get_end1(v);
            vntr.fuzzy_rl = vntr.fuzzy_end1-vntr.fuzzy_beg1+1;
            vntr.fuzzy_ll = vntr.fuzzy_rl; //??
        }
    }
    else if (mode==EXACT_LEFT_RIGHT_ALIGNMENT)
    {
        if (amode==EXACT)
        {
            extract_regions_by_exact_alignment(variant);
        }
        else
        {
            fprintf(stderr, "[E:%s] EXACT_LEFT_RIGHT_ALIGNMENT cannot be updated for any other attribute types beside EXACT.\n", __FUNCTION__);
            exit(1);
        }
    }
    else if (mode==FUZZY_LEFT_RIGHT_ALIGNMENT)
    {
        if (amode==FUZZY)
        {
            extract_regions_by_fuzzy_alignment(variant);
        }
        else
        {
            fprintf(stderr, "[E:%s] EXACT_LEFT_RIGHT_ALIGNMENT cannot be updated for any other attribute types beside FUZZY.\n", __FUNCTION__);
            exit(1);
        }
    }
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
std::string CandidateRegionExtractor::choose_repeat_unit(std::string& ref, std::string& motif)
{
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = VNTR::shift_str(motif, i);
        if (ref.compare(0, smotif.size(), smotif)==0)
        {
            return smotif;
        }
    }

    return motif;
}

/**
 * Checks if a vntr is a homopolymer.
 */
bool CandidateRegionExtractor::is_homopolymer(bcf_hdr_t* h, bcf1_t* v)
{
    bool is_homopolymer = false;
    uint32_t ref_len = strlen(bcf_get_ref(v));
    for (size_t i=1; i<bcf_get_n_allele(v); ++i)
    {
        std::string ref(bcf_get_alt(v, 0));
        std::string alt(bcf_get_alt(v, i));
        int32_t pos1 = bcf_get_pos1(v);
    }

    return is_homopolymer;
}

/**
 * Extract reference sequence region for motif discovery.
 *
 * The input is a VCF record that contains an indel.
 *
 * If the the indel has multiple alleles, it will examine all
 * alleles.
 *
 * todo: is might be a good idea to combine this step with motif detection
 *       since there seems to be a need to have an iterative process here
 *       to ensure a good candidate motif is chosen. *
 */
void CandidateRegionExtractor::extract_regions_by_exact_alignment(Variant& variant)
{
    if (debug)
    {
        if (debug) std::cerr << "********************************************\n";
        std::cerr << "EXTRACTIING REGION BY EXACT LEFT AND RIGHT ALIGNMENT\n\n";
    }

    bcf_hdr_t* h = variant.h;
    bcf1_t* v = variant.v;

    VNTR& vntr = variant.vntr;
    std::string& chrom = variant.chrom;

    int32_t min_beg1 = bcf_get_pos1(v);
    int32_t max_end1 = min_beg1;

    if (debug)
    {
       bcf_print_liten(h, v);
    }

    //merge candidate search region
    for (size_t i=1; i<bcf_get_n_allele(v); ++i)
    {
        std::string ref(bcf_get_alt(v, 0));
        std::string alt(bcf_get_alt(v, i));
        int32_t pos1 = bcf_get_pos1(v);

        //this prevents introduction of flanks that do not harbour the repeat unit
        trim(pos1, ref, alt);

        int32_t end1 = pos1 + ref.size() - 1;
        right_align(chrom, end1, ref, alt);

        int32_t beg1 = end1 - ref.size() + 1;
        left_align(chrom, beg1, ref, alt);

        min_beg1 = beg1<min_beg1 ? beg1 : min_beg1;
        max_end1 = end1>max_end1 ? end1 : max_end1;

        int32_t seq_len;
        char* seq = rs->fetch_seq(chrom, min_beg1-1, max_end1-1);

        if (debug)
        {
            std::cerr << "EXACT REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") from " << pos1 << ":" << ref << ":" << alt << "\n";
            std::cerr << "             " << seq << "\n";
        }

        if (seq) free(seq);
    }

    char* seq = rs->fetch_seq(chrom, min_beg1-1, max_end1-1);

    if (debug)
    {
        std::cerr << "FINAL EXACT REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
        std::cerr << "                   " << seq << "\n";
    }

    vntr.exact_repeat_tract = seq;
    vntr.rid = bcf_get_rid(v);
    vntr.exact_beg1 = min_beg1;
    vntr.exact_end1 = max_end1;

    if (seq) free(seq);
}

/**
 * Left align alleles.
 */
void CandidateRegionExtractor::left_align(std::string& chrom, int32_t& pos1, std::string& ref, std::string& alt)
{
    int32_t seq_len;
    while (ref.at(ref.size()-1)==alt.at(alt.size()-1) && pos1>1)
    {
        char base = rs->fetch_base(chrom, pos1-2);
        ref.erase(ref.size()-1,1);
        alt.erase(alt.size()-1,1);
        ref.insert(0, 1, base);
        alt.insert(0, 1, base);
        --pos1;
    }
}

/**
 * Right align alleles.
 */
void CandidateRegionExtractor::right_align(std::string& chrom, int32_t& pos1, std::string& ref, std::string& alt)
{
    while (ref.at(0)==alt.at(0))
    {
        char base = rs->fetch_base(chrom, pos1);
        ref.erase(0,1);
        alt.erase(0,1);
        ref.push_back(base);
        alt.push_back(base);
        ++pos1;
    }
}

/**
 * Extract reference sequence region for motif discovery in a fuzzy fashion.
 */
void CandidateRegionExtractor::extract_regions_by_fuzzy_alignment(Variant& variant)
{
    if (debug)
    {
        if (debug) std::cerr << "********************************************\n";
        std::cerr << "EXTRACTIING REGION BY FUZZY ALIGNMENT\n\n";
    }

    bcf_hdr_t* h = variant.h;
    bcf1_t* v = variant.v;
    std::string&  chrom = variant.chrom;

    int32_t min_beg1 = bcf_get_pos1(v);
    int32_t max_end1 = min_beg1;

    //merge candidate search region
    for (size_t i=1; i<bcf_get_n_allele(v); ++i)
    {
        std::string ref(bcf_get_alt(v, 0));
        std::string alt(bcf_get_alt(v, i));
        int32_t pos1 = bcf_get_pos1(v);

        trim(pos1, ref, alt);

        if (debug)
        {
            std::cerr << "indel fragment : " << (ref.size()<alt.size()? alt : ref) << "\n";
            std::cerr << "               : " << ref << ":" << alt << "\n";
        }

        min_beg1 = fuzzy_left_align(chrom, pos1, ref, alt, 3);
        max_end1 = fuzzy_right_align(chrom, pos1 + ref.size() - 1, ref, alt, 3);

        char* seq = rs->fetch_seq(chrom, min_beg1-1, max_end1-1);
        if (debug)
        {
            std::cerr << "FUZZY REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
            std::cerr << "             " << seq << "\n";
        }

        if (seq) free(seq);
    }

    char* seq = rs->fetch_seq(chrom, min_beg1-1, max_end1-1);

    if (debug)
    {
        std::cerr << "FINAL FUZZY REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
        std::cerr << "                   " << seq << "\n";
    }

    VNTR& vntr = variant.vntr;
    vntr.exact_repeat_tract.assign(seq);
    vntr.exact_beg1 = min_beg1;

    if (seq) free(seq);
}

/**
 * Fuzzy left align alleles allowing for mismatches and indels defined by penalty.
 *
 * @chrom   - chromosome
 * @pos1    - 1 based position
 * @ref     - reference sequence
 * @alt     - alternative sequence
 * @penalty - mismatch/indels allowed
 *
 * Returns left aligned position.
 */
uint32_t CandidateRegionExtractor::fuzzy_left_align(std::string& chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
{
    if (ref==alt)
    {
        return pos1;
    }

    //std::cerr << "fuzzy left alignment: " << chrom << ":" << pos1 << ":" << ref << ":" << alt << " (" <<  penalty << ")\n";
    while (ref.at(ref.size()-1)==alt.at(alt.size()-1) && pos1>1)
    {
        char base = rs->fetch_base(chrom, pos1-2);
        ref.erase(ref.size()-1,1);
        alt.erase(alt.size()-1,1);
        ref.insert(0, 1, base);
        alt.insert(0, 1, base);
        --pos1;
    }

    if (penalty)
    {
        uint32_t pos1_sub = pos1;
        uint32_t pos1_del = pos1;
        uint32_t pos1_ins = pos1;

        //substitution
        char base = rs->fetch_base(chrom, pos1-2);
        std::string new_ref = ref;
        std::string new_alt = alt;
        new_ref.erase(new_ref.size()-1,1);
        new_alt.erase(new_alt.size()-1,1);
        new_ref.insert(0, 1, base);
        new_alt.insert(0, 1, base);
//      std::cerr << "\tsub: " << chrom << ":" << pos1-1 << ":" << new_ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
        pos1_sub = fuzzy_left_align(chrom, pos1-1, new_ref, new_alt, penalty-1);

        //deletion
        if (ref.size()>1)
        {
            std::string new_ref = ref;
            new_ref.erase(new_ref.size()-1,1);
//            std::cerr << "\tdel: " << chrom << ":" << pos1 << ":" << new_ref << ":" << alt << " (" <<  penalty-1 << ")\n";
            pos1_del = fuzzy_left_align(chrom, pos1, new_ref, alt, penalty-1);
        }

        //insertion
        if (alt.size()>1)
        {
            std::string new_alt = alt;
            new_alt.erase(new_alt.size()-1,1);
//            std::cerr << "\tins: " << chrom << ":" << pos1 << ":" << ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_ins = fuzzy_left_align(chrom, pos1, ref, new_alt, penalty-1);
        }

        pos1 = std::min(pos1_sub, std::min(pos1_del, pos1_ins));
    }

    return pos1;
}

/**
 * Fuzzy right align alleles allowing for mismatches and indels defined by penalty.
 *
 * @chrom   - chromosome
 * @pos1    - 1 based position
 * @ref     - reference sequence
 * @alt     - alternative sequence
 * @penalty - mismatch/indels allowed
 *
 * Returns right aligned position.
 */
uint32_t CandidateRegionExtractor::fuzzy_right_align(std::string& chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
{
    if (ref==alt)
    {
        return pos1;
    }

    while (ref.at(0)==alt.at(0))
    {
        char base = rs->fetch_base(chrom, pos1);
        ref.erase(0,1);
        alt.erase(0,1);
        ref.push_back(base);
        alt.push_back(base);
        ++pos1;
    }

    if (penalty)
    {
        uint32_t pos1_sub = pos1;
        uint32_t pos1_del = pos1;
        uint32_t pos1_ins = pos1;

        //substitution
        char base = rs->fetch_base(chrom, pos1);
        std::string new_ref = ref;
        std::string new_alt = alt;
        new_ref.erase(0,1);
        new_alt.erase(0,1);
        new_ref.push_back(base);
        new_alt.push_back(base);
        //std::cerr << "\tsub: " << chrom << ":" << pos1+1 << ":" << new_ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
        pos1_sub = fuzzy_right_align(chrom, pos1+1, new_ref, new_alt, penalty-1);

        //deletion
        if (ref.size()>1)
        {
            std::string new_ref = ref;
            new_ref.erase(0,1);
            //std::cerr << "\tdel: " << chrom << ":" << pos1 << ":" << new_ref << ":" << alt << " (" <<  penalty-1 << ")\n";
            pos1_del = fuzzy_right_align(chrom, pos1, new_ref, alt, penalty-1);
        }

        //insertion
        if (alt.size()>1)
        {
            std::string new_alt = alt;
            new_alt.erase(0,1);
            //std::cerr << "\tins: " << chrom << ":" << pos1 << ":" << ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_ins = fuzzy_right_align(chrom, pos1, ref, new_alt, penalty-1);
        }

        pos1 = std::max(pos1_sub, std::max(pos1_del, pos1_ins));
    }

    return pos1;
}

/**
 * Trim alleles.
 */
void CandidateRegionExtractor::trim(int32_t& pos1, std::string& ref, std::string& alt)
{
    while (true)
    {
        if (ref.size()==1 || alt.size()==1)
        {
            break;
        }
        else if (ref.at(0)!=alt.at(0) && ref.at(ref.size()-1)!=alt.at(alt.size()-1))
        {
            break;
        }
        else
        {
            //trim from the right side
            if (ref.at(ref.size()-1)==alt.at(alt.size()-1))
            {
                ref.erase(ref.size()-1,1);
                alt.erase(alt.size()-1,1);
            }
            //trim from the left side
            else if (ref.at(0)==alt.at(0))
            {
                ref.erase(0,1);
                alt.erase(0,1);
                ++pos1;
            }

            //we choose one side to trim at a time to ensure that we do not accidentally end up with an empty allele
        }
    }
}