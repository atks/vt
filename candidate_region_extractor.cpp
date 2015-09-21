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
    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }

    max_mlen = 10;
    
   this->debug = debug;
    
    mt = new MotifTree(max_mlen, debug);
};

/**
 * Destructor.
 */
CandidateRegionExtractor::~CandidateRegionExtractor()
{
    delete vm;
    fai_destroy(fai);
    delete rfhmm;
    delete lfhmm;

    if (factors)
    {
        for (size_t i=1; i<=max_len; ++i)
        {
            free(factors[i]);
        }
        free(factors);
    }
}

/**
 * Pick candidate region.
 *
 * @mode - REFERENCE     use refence field
 *       - ALLELE_EXACT  by exact alignment
 *       - ALLELE_FUZZY  by fuzzy alignment
 */
void CandidateRegionExtractor::pick_candidate_region(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr, uint32_t mode)
{
    if (mode==REFERENCE)
    {
        vntr.repeat_tract.assign(bcf_get_ref(v));
        vntr.rbeg1 = bcf_get_pos1(v);
    }
    else if (mode==EXACT_LEFT_RIGHT_ALIGNMENT)
    {
        extract_regions_by_exact_alignment(h, v, vntr);
    }
    else if (mode==FUZZY_LEFT_RIGHT_ALIGNMENT)
    {
        extract_regions_by_fuzzy_alignment(h, v, vntr);
    }
    else if (mode==FUZZY_LEFT_RIGHT_ALIGNMENT_WITH_PENALTY)
    {
        extract_regions_by_fuzzy_alignment_with_penalty(h, v, vntr);
    }
}

/**
 * Detect repeat region.
 */
void CandidateRegionExtractor::detect_repeat_region(bcf_hdr_t* h, bcf1_t *v, Variant& variant, uint32_t mode)
{
    VNTR& vntr = variant.vntr;
    
    //simple single base pair clipping of ends
    if (mode==CLIP_ENDS)
    {
        if (debug)
        {
            std::cerr << "********************************************\n";
            std::cerr << "CLIP ENDS\n";
            std:: cerr << "\n";
        }

        if (vntr.repeat_tract.size()>2)
        {
            vntr.repeat_tract = vntr.repeat_tract.substr(1, vntr.repeat_tract.size()-2);
            ++vntr.rbeg1;
        }

        vntr.ru = choose_repeat_unit(vntr.repeat_tract, vntr.motif);
        vntr.rl = (float)vntr.repeat_tract.size()/vntr.ru.size();

        if (debug)
        {
            vntr.print();
        }
    }
    //simple single base pair clipping of ends
    else if (mode==CLIP_1L2R)
    {
        if (debug)
        {
            std::cerr << "********************************************\n";
            std::cerr << "CLIP ENDS\n";
            std:: cerr << "\n";
        }

        if (vntr.repeat_tract.size()>3)
        {
            vntr.repeat_tract = vntr.repeat_tract.substr(1, vntr.repeat_tract.size()-3);
            ++vntr.rbeg1;
        }

        vntr.ru = choose_repeat_unit(vntr.repeat_tract, vntr.motif);
        vntr.rl = (float)vntr.repeat_tract.size()/vntr.ru.size();

        if (debug)
        {
            vntr.print();
        }
    }
    
    //fill in flanks
    const char* chrom = variant.chrom.c_str();
    uint32_t pos1 = vntr.rbeg1;
    int32_t len = 0;
    faidx_fetch_seq(fai, chrom, pos1-10, pos1-1, &len);
    
};

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
std::string CandidateRegionExtractor::choose_repeat_unit(std::string& ref, std::string& motif)
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
 *
 *  
 */
void CandidateRegionExtractor::extract_regions_by_exact_alignment(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr)
{
    if (debug)
    {
        if (debug) std::cerr << "********************************************\n";
        std::cerr << "EXTRACTIING REGION BY EXACT LEFT AND RIGHT ALIGNMENT\n\n";
    }

    const char* chrom = bcf_get_chrom(h, v);

    int32_t min_beg1 = bcf_get_pos1(v);
    int32_t max_end1 = min_beg1;

    //merge candidate search region
    for (size_t i=1; i<bcf_get_n_allele(v); ++i)
    {
        std::string ref(bcf_get_alt(v, 0));
        std::string alt(bcf_get_alt(v, i));
        int32_t pos1 = bcf_get_pos1(v);

        //why do this??
        trim(pos1, ref, alt);

        if (debug)
        {
           bcf_print_liten(h, v);
        }

        int32_t end1 = pos1 + ref.size() - 1;
        right_align(chrom, end1, ref, alt);

        int32_t beg1 = end1 - ref.size() + 1;
        left_align(chrom, beg1, ref, alt);

        min_beg1 = beg1<min_beg1 ? beg1 : min_beg1;
        max_end1 = end1>max_end1 ? end1 : max_end1;

        int32_t seq_len;
        char* seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);

        if (debug)
        {
            std::cerr << "EXACT REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
            std::cerr << "             " << seq << "\n";
        }

        if (seq_len) free(seq);
    }

    int32_t seq_len;
    char* seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);

    if (debug)
    {
        std::cerr << "FINAL EXACT REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
        std::cerr << "                   " << seq << "\n";
    }

    vntr.repeat_tract = seq;
    vntr.rid = bcf_get_rid(v);
    vntr.rbeg1 = min_beg1;
    vntr.rend1 = max_end1;
    
    if (seq_len) free(seq);
}

/**
 * Left align alleles.
 */
void CandidateRegionExtractor::left_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt)
{
    int32_t seq_len;
    char* seq;
    while (ref.at(ref.size()-1)==alt.at(alt.size()-1) && pos1>1)
    {
        seq = faidx_fetch_seq(fai, chrom, pos1-2, pos1-2, &seq_len);
        if (seq_len)
        {
            ref.erase(ref.size()-1,1);
            alt.erase(alt.size()-1,1);
            ref.insert(0, 1, seq[0]);
            alt.insert(0, 1, seq[0]);
            free(seq);
            --pos1;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }
}

/**
 * Right align alleles.
 */
void CandidateRegionExtractor::right_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt)
{
    int32_t seq_len;
    char* seq;
    while (ref.at(0)==alt.at(0))
    {
        seq = faidx_fetch_seq(fai, chrom, pos1, pos1, &seq_len);
        if (seq_len)
        {
            ref.erase(0,1);
            alt.erase(0,1);
            ref.push_back(seq[0]);
            alt.push_back(seq[0]);
            free(seq);
            ++pos1;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }
}

/**
 * Extract reference sequence region for motif discovery in a fuzzy fashion.
 */
void CandidateRegionExtractor::extract_regions_by_fuzzy_alignment(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr)
{
    if (debug)
    {
        if (debug) std::cerr << "********************************************\n";
        std::cerr << "EXTRACTIING REGION BY FUZZY ALIGNMENT\n\n";
    }

    const char* chrom = bcf_get_chrom(h, v);

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

        int32_t seq_len;
        char* seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);
        if (debug)
        {
            std::cerr << "FUZZY REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
            std::cerr << "             " << seq << "\n";
        }

        if (seq_len) free(seq);
    }

    int32_t seq_len;
    char* seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);

    if (debug)
    {
        std::cerr << "FINAL FUZZY REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
        std::cerr << "                   " << seq << "\n";
    }

    vntr.repeat_tract = seq;
    vntr.rbeg1 = min_beg1;

    if (seq_len) free(seq);
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
uint32_t CandidateRegionExtractor::fuzzy_left_align(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
{
    if (ref==alt)
    {
        return pos1;
    }

    //std::cerr << "fuzzy left alignment: " << chrom << ":" << pos1 << ":" << ref << ":" << alt << " (" <<  penalty << ")\n";
    int32_t seq_len;
    char* seq;
    while (ref.at(ref.size()-1)==alt.at(alt.size()-1) && pos1>1)
    {
        seq = faidx_fetch_seq(fai, chrom, pos1-2, pos1-2, &seq_len);
        if (seq_len)
        {
            ref.erase(ref.size()-1,1);
            alt.erase(alt.size()-1,1);
            ref.insert(0, 1, seq[0]);
            alt.insert(0, 1, seq[0]);
            free(seq);
            --pos1;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }

    if (penalty)
    {
        uint32_t pos1_sub = pos1;
        uint32_t pos1_del = pos1;
        uint32_t pos1_ins = pos1;

        //substitution
        seq = faidx_fetch_seq(fai, chrom, pos1-2, pos1-2, &seq_len);
        if (seq_len)
        {
            std::string new_ref = ref;
            std::string new_alt = alt;
            new_ref.erase(new_ref.size()-1,1);
            new_alt.erase(new_alt.size()-1,1);
            new_ref.insert(0, 1, seq[0]);
            new_alt.insert(0, 1, seq[0]);
//            std::cerr << "\tsub: " << chrom << ":" << pos1-1 << ":" << new_ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_sub = fuzzy_left_align(chrom, pos1-1, new_ref, new_alt, penalty-1);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }

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
uint32_t CandidateRegionExtractor::fuzzy_right_align(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
{
    if (ref==alt)
    {
        return pos1;
    }

    int32_t seq_len;
    char* seq;
    while (ref.at(0)==alt.at(0))
    {
        seq = faidx_fetch_seq(fai, chrom, pos1, pos1, &seq_len);
        if (seq_len)
        {
            ref.erase(0,1);
            alt.erase(0,1);
            ref.push_back(seq[0]);
            alt.push_back(seq[0]);
            free(seq);
            ++pos1;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }

    if (penalty)
    {
        uint32_t pos1_sub = pos1;
        uint32_t pos1_del = pos1;
        uint32_t pos1_ins = pos1;

        //substitution
        seq = faidx_fetch_seq(fai, chrom, pos1, pos1, &seq_len);
        if (seq_len)
        {
            std::string new_ref = ref;
            std::string new_alt = alt;
            new_ref.erase(0,1);
            new_alt.erase(0,1);
            new_ref.push_back(seq[0]);
            new_alt.push_back(seq[0]);
            //std::cerr << "\tsub: " << chrom << ":" << pos1+1 << ":" << new_ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_sub = fuzzy_right_align(chrom, pos1+1, new_ref, new_alt, penalty-1);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }

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
 * Extract reference sequence region for motif discovery in a fuzzy fashion.
 */
void CandidateRegionExtractor::extract_regions_by_fuzzy_alignment_with_penalty(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr)
{
    if (debug)
    {
        if (debug) std::cerr << "********************************************\n";
        std::cerr << "EXTRACTIING REGION BY FUZZY ALIGNMENT\n\n";
    }

    const char* chrom = bcf_get_chrom(h, v);

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

        min_beg1 = fuzzy_left_align_with_penalty(chrom, pos1, ref, alt, 3);
        max_end1 = fuzzy_right_align_with_penalty(chrom, pos1 + ref.size() - 1, ref, alt, 3);

        int32_t seq_len;
        char* seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);
        if (debug)
        {
            std::cerr << "FUZZY REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
            std::cerr << "             " << seq << "\n";
        }

        if (seq_len) free(seq);
    }

    int32_t seq_len;
    char* seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);

    if (debug)
    {
        std::cerr << "FINAL FUZZY REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
        std::cerr << "                   " << seq << "\n";
    }

    vntr.repeat_tract = seq;
    vntr.rbeg1 = min_beg1;

    if (seq_len) free(seq);
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
uint32_t CandidateRegionExtractor::fuzzy_left_align_with_penalty(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
{
    if (ref==alt)
    {
        return pos1;
    }

    //std::cerr << "fuzzy left alignment: " << chrom << ":" << pos1 << ":" << ref << ":" << alt << " (" <<  penalty << ")\n";
    int32_t seq_len;
    char* seq;
    int32_t shift = 0;
    while (ref.at(ref.size()-1)==alt.at(alt.size()-1) && pos1>1)
    {
        seq = faidx_fetch_seq(fai, chrom, pos1-2, pos1-2, &seq_len);
        if (seq_len)
        {
            ref.erase(ref.size()-1,1);
            alt.erase(alt.size()-1,1);
            ref.insert(0, 1, seq[0]);
            alt.insert(0, 1, seq[0]);
            free(seq);
            --pos1;
            ++shift;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }

    if (shift>1 && penalty)
    {
        uint32_t pos1_sub = pos1;
        uint32_t pos1_del = pos1;
        uint32_t pos1_ins = pos1;

        //substitution
        seq = faidx_fetch_seq(fai, chrom, pos1-2, pos1-2, &seq_len);
        if (seq_len)
        {
            std::string new_ref = ref;
            std::string new_alt = alt;
            new_ref.erase(new_ref.size()-1,1);
            new_alt.erase(new_alt.size()-1,1);
            new_ref.insert(0, 1, seq[0]);
            new_alt.insert(0, 1, seq[0]);
//            std::cerr << "\tsub: " << chrom << ":" << pos1-1 << ":" << new_ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_sub = fuzzy_left_align_with_penalty(chrom, pos1-1, new_ref, new_alt, penalty-1);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }

        //deletion
        if (ref.size()>1)
        {
            std::string new_ref = ref;
            new_ref.erase(new_ref.size()-1,1);
//            std::cerr << "\tdel: " << chrom << ":" << pos1 << ":" << new_ref << ":" << alt << " (" <<  penalty-1 << ")\n";
            pos1_del = fuzzy_left_align_with_penalty(chrom, pos1, new_ref, alt, penalty-1);
        }

        //insertion
        if (alt.size()>1)
        {
            std::string new_alt = alt;
            new_alt.erase(new_alt.size()-1,1);
//            std::cerr << "\tins: " << chrom << ":" << pos1 << ":" << ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_ins = fuzzy_left_align_with_penalty(chrom, pos1, ref, new_alt, penalty-1);
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
uint32_t CandidateRegionExtractor::fuzzy_right_align_with_penalty(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
{
    if (ref==alt)
    {
        return pos1;
    }

    int32_t seq_len;
    char* seq;
    int32_t shift = 0;
    while (ref.at(0)==alt.at(0))
    {
        seq = faidx_fetch_seq(fai, chrom, pos1, pos1, &seq_len);
        if (seq_len)
        {
            ref.erase(0,1);
            alt.erase(0,1);
            ref.push_back(seq[0]);
            alt.push_back(seq[0]);
            free(seq);
            ++pos1;
            ++shift;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }

    if (penalty==3) shift--;;

    if (shift>1 && penalty)
    {
        uint32_t pos1_sub = pos1;
        uint32_t pos1_del = pos1;
        uint32_t pos1_ins = pos1;

        //substitution
        seq = faidx_fetch_seq(fai, chrom, pos1, pos1, &seq_len);
        if (seq_len)
        {
            std::string new_ref = ref;
            std::string new_alt = alt;
            new_ref.erase(0,1);
            new_alt.erase(0,1);
            new_ref.push_back(seq[0]);
            new_alt.push_back(seq[0]);
            //std::cerr << "\tsub: " << chrom << ":" << pos1+1 << ":" << new_ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_sub = fuzzy_right_align_with_penalty(chrom, pos1+1, new_ref, new_alt, penalty-1);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }

        //deletion
        if (ref.size()>1)
        {
            std::string new_ref = ref;
            new_ref.erase(0,1);
            //std::cerr << "\tdel: " << chrom << ":" << pos1 << ":" << new_ref << ":" << alt << " (" <<  penalty-1 << ")\n";
            pos1_del = fuzzy_right_align_with_penalty(chrom, pos1, new_ref, alt, penalty-1);
        }

        //insertion
        if (alt.size()>1)
        {
            std::string new_alt = alt;
            new_alt.erase(0,1);
            //std::cerr << "\tins: " << chrom << ":" << pos1 << ":" << ref << ":" << new_alt << " (" <<  penalty-1 << ")\n";
            pos1_ins = fuzzy_right_align_with_penalty(chrom, pos1, ref, new_alt, penalty-1);
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
            if (ref.at(ref.size()-1)==alt.at(alt.size()-1))
            {
                ref.erase(ref.size()-1,1);
                alt.erase(alt.size()-1,1);
            }
            else if (ref.at(0)==alt.at(0))
            {
                ref.erase(0,1);
                alt.erase(0,1);
                ++pos1;
            }
        }
    }
}