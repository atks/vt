/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#include "vntr_annotator.h"

/**
 * Constructor.
 */
VNTRAnnotator::VNTRAnnotator(std::string& ref_fasta_file, bool debug)
{
    vm = new VariantManip(ref_fasta_file.c_str());
    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }

    max_mlen = 10;
    mt = new MotifTree(max_mlen, debug);


//    //update factors
//    factors = NULL;
//    initialize_factors(32);

    this->debug = debug;

    qual.assign(256, 'K');

    ahmm = new AHMM();
    rfhmm = new RFHMM();
    lfhmm = new LFHMM();
};

/**
 * Destructor.
 */
VNTRAnnotator::~VNTRAnnotator()
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
 * Annotates VNTR characteristics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
void VNTRAnnotator::annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant, std::string mode)
{
    if (variant.type==VT_VNTR)
    {
        if (debug) std::cerr << "ANNOTATING VNTR/STR \n";

        //1. detect candidate motifs from a reference seqeuence
        pick_candidate_motifs(h, v, REFERENCE);

        //2. choose the best candidate motif
        choose_best_motif(h, v, mt, REFERENCE);

        return;
    }
    else if (variant.type&VT_INDEL)
    {
        //1. selects candidate region by exact left and right alignment
        //2. detects motif
        //3. detects left and right flank
        if (mode=="e")
        {

            //1. detect candidate motifs from a reference seqeuence
            pick_candidate_motifs(h, v, ALLELE_EXACT);

            //2. choose the best candidate motif


            if (debug) std::cerr << "============================================\n";
            return;
        }
        //1. selects candidate region by fuzzy left and right alignment
        //2. detects motif
        //3. detects left and right flank
        else if (mode=="f")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL FUZZILY\n";


            pick_candidate_motifs(h, v, ALLELE_FUZZY);
            
            if (!mt->pcm.empty())
            {
                CandidateMotif cm = mt->pcm.top();
                variant.vntr.motif = cm.motif;
                variant.vntr.motif_score = cm.score;
//                variant.vntr.pos1 = region.beg1;
//                variant.vntr.ref = region.ref;

                if (debug)
                {
                    uint32_t n =0;
                    while(!mt->pcm.empty())
                    {
                        CandidateMotif cm = mt->pcm.top();
                        std::cerr << cm.motif << " " << cm.score << "\n";
                        mt->pcm.pop();
                        if (n==10)
                            break;

                        ++n;
                    }
                }
            }


            return;
        }
        //1. selects candidate region by fuzzy left and right alignment
        //2. detects motif
        //3. detects left anf right flank
        else if (mode=="x")
        {
            //treat homopolymers as a special case
            if (is_homopolymer(h,v))
            {
            }
            //all other subject to
            else
            {
            }
        }

    }
}

/**
 * Pick candidate motifs in different modes.
 * Invokes motif tree and the candidate motifs are stored in a
 * heap within the motif tree.
 */
void VNTRAnnotator::pick_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, uint32_t mode)
{
    if (mode==REFERENCE)
    {
        mt->detect_candidate_motifs(bcf_get_ref(v), strlen(bcf_get_ref(v)));
    }
    else if (mode==ALLELE_EXACT)
    {
        ReferenceRegion region = extract_regions_by_exact_alignment(h, v);
        mt->detect_candidate_motifs(region.ref);
    }
    else if (mode==ALLELE_FUZZY)
    {
        ReferenceRegion region = extract_regions_by_fuzzy_alignment(h, v);
        mt->detect_candidate_motifs(region.ref);
    }
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 *
 */
VNTR VNTRAnnotator::choose_best_motif(bcf_hdr_t* h, bcf1_t* v, MotifTree* mt, uint32_t mode)
{
    if (debug) std::cerr << "pcm size : " << mt->pcm.size() << "\n";
    if (debug) std::cerr << "*********************************************\n";

    if (debug) std::cerr << "============================================\n";

    if (mode==REFERENCE)
    {
        //choose candidate motif
        bool first = true;
        float cp = 0;
        float mscore = 0;
        uint32_t clen = 0;
        std::string ref(bcf_get_ref(v));
        if (ref.size()>256) ref = ref.substr(0,256);

        uint32_t no = 0;

        VNTR vntr;

        while (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            if (first)
            {
                std::string ru = choose_repeat_unit(ref, cm.motif);
                ahmm->set_model(ru.c_str());
                ahmm->align(ref.c_str(), qual.c_str());

                vntr.motif = cm.motif;
                vntr.motif_score = ahmm->get_motif_concordance();;
                vntr.ru = ru;

                cp = cm.score;
                clen = cm.len;
                first = false;
                mscore = ahmm->get_motif_concordance();

                if (debug)
                {
                    //ahmm->print_alignment();

                    printf("    selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", cm.motif.c_str(),
                                                                    cm.score,
                                                                    cm.fit,
                                                                    ahmm->get_motif_concordance(),
                                                                    (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                                    ahmm->get_exact_motif_count(),
                                                                    ahmm->get_motif_count());
                }
            }
            else
            {
                if (no<3 && (cp-cm.score<((float)1/cm.len)*cp || cm.score>0.4))
                {
                    //if score are not perfect, use AHMM
                    std::string ru = choose_repeat_unit(ref, cm.motif);
                    ahmm->set_model(ru.c_str());
                    ahmm->align(ref.c_str(), qual.c_str());

                    cp = cm.score;
                    clen = cm.len;

                    if (ahmm->get_motif_concordance()>mscore)
                    {
                        vntr.motif = cm.motif;
                        vntr.motif_score = ahmm->get_motif_concordance();;
                        vntr.ru = ru;
                    }

                    if (debug)
                    {
                        //ahmm->print_alignment();

                        printf("    selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score,
                                                                        mt->pcm.top().fit,
                                                                        ahmm->get_motif_concordance(),
                                                                        (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                                        ahmm->get_exact_motif_count(),
                                                                        ahmm->get_motif_count());
                    }
                }
                else
                {
                    if (debug)
                    {
                        printf("not selected: %10s %.2f \n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score);
                    }
                }
            }

            if (!mt->pcm.empty()) mt->pcm.pop();

            ++no;
        }

        return vntr;
    }
    else if (mode==ALLELE_EXACT)
    {
        if (debug) std::cerr << "ANNOTATING INDEL EXACTLY\n";

        VNTR vntr;

        if (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            vntr.motif = cm.motif;
            vntr.motif_score = cm.score;
        }
    }
    else if (mode==ALLELE_FUZZY)
    {
        //choose candidate motif
        bool first = true;
        float cp = 0;
        float mscore = 0;
        uint32_t clen = 0;
        std::string ref(bcf_get_ref(v));
        if (ref.size()>256) ref = ref.substr(0,256);

        uint32_t no = 0;

        VNTR vntr;

        while (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            if (first)
            {
                std::string ru = choose_repeat_unit(ref, cm.motif);
                ahmm->set_model(ru.c_str());
                ahmm->align(ref.c_str(), qual.c_str());

                vntr.motif = cm.motif;
                vntr.motif_score = ahmm->get_motif_concordance();;
                vntr.ru = ru;

                cp = cm.score;
                clen = cm.len;
                first = false;
                mscore = ahmm->get_motif_concordance();

                if (debug)
                {
                    ahmm->print_alignment();

                    printf("    selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", cm.motif.c_str(),
                                                                    cm.score,
                                                                    cm.fit,
                                                                    ahmm->get_motif_concordance(),
                                                                    (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                                    ahmm->get_exact_motif_count(),
                                                                    ahmm->get_motif_count());
                }
            }
            else
            {
                if (no<3 && (cp-cm.score<((float)1/cm.len)*cp || cm.score>0.4))
                {
                    //if score are not perfect, use AHMM
                    std::string ru = choose_repeat_unit(ref, cm.motif);
                    ahmm->set_model(ru.c_str());
                    ahmm->align(ref.c_str(), qual.c_str());

                    cp = cm.score;
                    clen = cm.len;

                    if (ahmm->get_motif_concordance()>mscore)
                    {
                        vntr.motif = cm.motif;
                        vntr.motif_score = ahmm->get_motif_concordance();;
                        vntr.ru = ru;
                    }

                    if (debug)
                    {
                        ahmm->print_alignment();

                        printf("    selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score,
                                                                        mt->pcm.top().fit,
                                                                        ahmm->get_motif_concordance(),
                                                                        (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                                        ahmm->get_exact_motif_count(),
                                                                        ahmm->get_motif_count());
                    }
                }
                else
                {
                    if (debug)
                    {
                        printf("not selected: %10s %.2f \n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score);
                    }
                }
            }

            if (!mt->pcm.empty()) mt->pcm.pop();

            ++no;
        }

        return vntr;
    }

    if (debug) std::cerr << "============================================\n";

    VNTR vntr;
    return vntr;
}

///**
// * Constructor.
// */
//void VNTRAnnotator::initialize_factors(int32_t max_len)
//{
//    if (factors)
//    {
//        for (size_t i=1; i<=this->max_len; ++i)
//        {
//            delete factors[i];
//        }
//        free(factors);
//    }
//
//    this->max_len = max_len;
//
//    factors = new int32_t*[max_len+1];
//    for (size_t i=1; i<=max_len; ++i)
//    {
//        factors[i] = new int32_t[max_len];
//        int32_t count = 0;
//
//        for (size_t j=1; j<=i; ++j)
//        {
//            if ((i%j)==0)
//            {
//                factors[i][count++] = j;
//            }
//        }
//    }
//};

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
std::string VNTRAnnotator::choose_repeat_unit(std::string& ref, std::string& motif)
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
 * Checks if a variant is a homopolymer.
 */
bool VNTRAnnotator::is_homopolymer(bcf_hdr_t* h, bcf1_t* v)
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
 * Infer flanks  motif discovery.
 *
 * returns
 * a. motif concordance
 * b. purity concordance
 * c. left flank
 * d. right flank
 */
void VNTRAnnotator::infer_flanks(bcf_hdr_t* h, bcf1_t* v, std::string& motif)
{
    //fit into flanks

    //left alignment

    //right alignment


    //        int32_t no_candidate_motifs;
//
//
//        const char* chrom = bcf_get_chrom(h, v);
//
//
//
//        int32_t pos1 = bcf_get_pos1(v);
//
//
//
//
//        std::vector<CandidateMotif> candidate_motifs;
//        pick_candidate_motifs(h, v, candidate_motifs);
//
//        std::string eru = candidate_motifs[0].motif;
//        std::string emotif = get_motif(eru);
//
//        std::cerr << "exact motif : " << emotif << "\n";
//        std::cerr << "exact ru    : " << eru << "\n";
//
//        variant.eru = eru;
//        variant.emotif = emotif;
//
//        ///////////////
//        //INEXACT STRs
//        //////////////
//
//        bool debug = false;
//
////        std::cerr << "\t";
////        bcf_print(h, v);
//
//        if (debug)
//        {
//            std::cerr << "//////////////////////\n";
//            std::cerr << "//FUZZY LEFT ALIGNMENT\n";
//            std::cerr << "//////////////////////\n";
//        }
//
//        int32_t lflank_len, rflank_len, lref_len, rref_len;
//        char* lflank, *ru, *rflank, *lrefseq, *rrefseq;
//        int32_t ru_len;
//        std::string qual(2048, 'K');
//
//        int32_t flank_len = 10;
//        int32_t check_len = 100;
//
//        //left flank
//        lflank = faidx_fetch_seq(fai, chrom, pos1-flank_len, pos1-1, &lflank_len);
//        lrefseq = faidx_fetch_seq(fai, chrom, pos1-flank_len, pos1+check_len, &lref_len);
//
//        int32_t penalty = emotif.size()>6? 5 : emotif.size();
//        penalty = 10;
//        //std::cerr << "PENALTY : " << penalty << "\n";
//
//        if (pos1==164284)
//        {
//            //std::cerr << "MOTIF AND LFLANK  " << emotif << " " << lflank << "\n";
//        }
//        lfhmm->set_model(lflank, eru.c_str());
//        lfhmm->set_mismatch_penalty(penalty);
//        lfhmm->set_delta(0.0000000001);
//        lfhmm->initialize_T();
//        lfhmm->align(lrefseq, qual.c_str());
//        if (debug) lfhmm->print_alignment();

//    std::cerr << "model: " << "(" << lflank_start[MODEL] << "~" << lflank_end[MODEL] << ") "
//                          << "[" << motif_start[MODEL] << "~" << motif_end[MODEL] << "]\n";
//    std::cerr << "read : " << "(" << lflank_start[READ] << "~" << lflank_end[READ] << ") "
//                          << "[" << motif_start[READ] << "~" << motif_end[READ] << "]"
//                          << "[" << rflank_start[READ] << "~" << rflank_end[READ] << "]\n";

//        int32_t lflank_start[2], lflank_end[2], motif_start[2], motif_end[2], rflank_start[2], rflank_end[2];
//        int32_t motif_count, exact_motif_count, motif_m, motif_xid;

//
//        if (debug) std::cerr << "\n#################################\n\n";
//
//        if (debug)
//        {
//            std::cerr << "///////////////////////\n";
//            std::cerr << "//FUZZY RIGHT ALIGNMENT\n";
//            std::cerr << "///////////////////////\n";
//        }
//
//
//        int32_t motif_end1 = pos1+lfhmm->get_motif_read_epos1()-flank_len;
//
//        if (pos1==164284)
//        {
////            rfhmm->set_debug(true);
//        }

        //right flank
//        rflank = faidx_fetch_seq(fai, chrom, motif_end1, motif_end1+flank_len-1, &rflank_len);
//        rrefseq = faidx_fetch_seq(fai, chrom, motif_end1-check_len, motif_end1+flank_len-1, &rref_len);
//        if (pos1==164284)
//        {
//        //    std::cerr << "MOTIF AND RFLANK  " << emotif << " " << rflank << "\n";
//        }
//        //std::cerr << "MOTIF AND RFLANK  " << emotif << " " << rflank << "\n";
//
//        rfhmm->set_model(eru.c_str(), rflank);
//        rfhmm->set_mismatch_penalty(penalty);
//        rfhmm->set_delta(0.0000000001);
//        rfhmm->initialize_T();
//        rfhmm->align(rrefseq, qual.c_str());
//        if (debug) rfhmm->print_alignment();
//
//        if (debug)
//        {
//            std::cerr << "pos1 " << pos1 << "\n";
//            std::cerr << "motif end1 " << motif_end1 << "\n";
//            std::cerr << "read spos1 " << rfhmm->get_motif_read_spos1() << "\n";
//            std::cerr << "read epos1 " << rfhmm->get_motif_read_epos1() << "\n";
//        }
//
//        int32_t motif_beg1 = motif_end1-check_len + rfhmm->get_motif_read_spos1();
//
//        variant.iregion.beg1 = motif_beg1;
//        variant.iregion.end1 = motif_end1;
//        if (pos1==164284)
//        {
//         //   exit(1);
//        }
//
//        if (lflank_len) free(lflank);
//        if (lref_len) free(lrefseq);
//        if (rflank_len) free(rflank);
//        if (rref_len) free(rrefseq);


}

/**
 * Extract reference sequence region for motif discovery.
 */
ReferenceRegion VNTRAnnotator::extract_regions_by_exact_alignment(bcf_hdr_t* h, bcf1_t* v)
{
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
            std::cerr << "indel fragment : " << (ref.size()<alt.size()? alt : ref) << "\n";
            std::cerr << "               : " << ref << ":" << alt << "\n";
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

    ReferenceRegion region = ReferenceRegion(min_beg1, max_end1, seq);

    if (seq_len) free(seq);

    return region;
}

/**
 * Pick candidate motifs.
 * candidate_motifs contain motifs and a measure of confidence
 */
//void VNTRAnnotator::pick_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, std::vector<CandidateMotif>& candidate_motifs)
//{
//    const char* chrom = bcf_get_chrom(h, v);
//
//    bcf_print(h,v);
//
//    int32_t min_beg1 = bcf_get_pos1(v);
//    int32_t max_end1 = min_beg1;
//
//    //merge candidate search region
//    for (size_t i=1; i<bcf_get_n_allele(v); ++i)
//    {
//        std::string ref(bcf_get_alt(v, 0));
//        std::string alt(bcf_get_alt(v, i));
//        int32_t pos1 = bcf_get_pos1(v);
//
//        /////////////////////////////////////////
//        //EXACT STRs in inserted/deleted sequence
//        /////////////////////////////////////////
////        std::cerr << "before trimming : " << pos1 << " " << ref << "\n";
////        std::cerr << "                : " << pos1 << " " << alt << "\n";
//
//        trim(pos1, ref, alt);
////        std::cerr << "after trimming : " << pos1 << " " << ref << "\n";
////        std::cerr << "               : " << pos1 << " " << alt << "\n";
//
//        std::cerr << "indel fragment : " << (ref.size()<alt.size()? alt : ref) << "\n";
//        std::cerr << "               : " << ref << ":" << alt << "\n";
//
//
//        int32_t end1 = pos1 + ref.size() - 1;
//        right_align(chrom, end1, ref, alt);
////        std::cerr << "after right alignment : " << end1 << " " << ref << "\n";
////        std::cerr << "                      : " << end1 << " " << alt << "\n";
//
//        //do this later as you want to pefrom the alignment shortly.
//        int32_t beg1 = end1 - ref.size() + 1;
//        left_align(chrom, beg1, ref, alt);
////        std::cerr << "after left alignment : " << beg1 << " " << ref << "\n";
////        std::cerr << "                     : " << beg1 << " " << alt << "\n";
//
//        min_beg1 = beg1<min_beg1 ? beg1 : min_beg1;
//        max_end1 = end1>max_end1 ? end1 : max_end1;
//
//        int32_t seq_len;
//        char* seq;
//        seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);
//        std::cerr << "EXACT REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
//        std::cerr << "             " << seq << "\n";
//        if (seq_len) free(seq);
//    }
//
//    int32_t seq_len;
//    char* seq;
//    seq = faidx_fetch_seq(fai, chrom, min_beg1-1, max_end1-1, &seq_len);
//
//    std::cerr << "EXACT REGION " << min_beg1 << "-" << max_end1 << " (" << max_end1-min_beg1+1 <<") " << "\n";
//    std::cerr << "             " << seq << "\n";
//
//    //detect motif
//    //mt->set_sequence(seq);
//  //  mt->detect_candidate_motifs(candidate_motifs, 6);
//
//    std::string sequence(seq);
//
//    int32_t bases[4] = {0,0,0,0};
//    for (size_t i=0; i<sequence.size(); ++i)
//    {
//        if (seq[i]=='A')
//        {
//            ++bases[0];
//        }
//        else if (seq[i]=='C')
//        {
//            ++bases[1];
//        }
//        else if (seq[i]=='G')
//        {
//            ++bases[2];
//        }
//        else if (seq[i]=='T')
//        {
//            ++bases[3];
//        }
//    }
//
//    std::cerr << "             A " << ((float)bases[0]/sequence.size()) << "\n";
//    std::cerr << "             C " << ((float)bases[1]/sequence.size()) << "\n";
//    std::cerr << "             G " << ((float)bases[2]/sequence.size()) << "\n";
//    std::cerr << "             T " << ((float)bases[3]/sequence.size()) << "\n";
//
//    if (seq_len>2)
//    {
//        sequence = sequence.substr(1, sequence.size()-2);
//    }
//
//    if (seq_len) free(seq);
//
//
//
//    scan_exact_motif(sequence);
//
//}

/**
 * This is a quick scan for a motif that is exactly repeated.
 */
std::string VNTRAnnotator::scan_exact_motif(std::string& sequence)
{
    size_t i = 0;
    size_t len = sequence.size();
    size_t sub_motif_len;
    const char* seq = sequence.c_str();
    size_t d = 0;
    size_t max_sub_motif_len = len/2;
    for (sub_motif_len = 1; sub_motif_len<=max_sub_motif_len; ++sub_motif_len)
    {
        size_t n_sub_motif = len/sub_motif_len;

        bool exact = true;
        size_t concordant = 0;
        size_t c1 = 0;

        for (size_t j=0; j<n_sub_motif; ++j)
        {
            if ((strncmp(&seq[0], &seq[j*sub_motif_len], sub_motif_len)))
            {
                exact = false;

                for (size_t k=0; k<sub_motif_len; ++k)
                {
                    if (seq[j*sub_motif_len+k]==seq[k])
                    {
                       ++concordant;
                    }
                }
            }
            else
            {
                ++c1;
                concordant += sub_motif_len;
            }
        }

        if (n_sub_motif*sub_motif_len<len)
        {
            if (strncmp(&seq[0], &seq[n_sub_motif*sub_motif_len], len-n_sub_motif*sub_motif_len))
            {
                exact = false;
            }
        }

        if (exact)
        {
            return sequence.substr(0, sub_motif_len);
        }
    }

    return "";
}


/**
 * Pick shortest consensus motif.
 */
std::string VNTRAnnotator::pick_consensus_motif(std::string& ref)
{
    //maybe not necessary later



    return "";
}

/**
 * Trim alleles.
 */
void VNTRAnnotator::trim(int32_t& pos1, std::string& ref, std::string& alt)
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

/**
 * Left align alleles.
 */
void VNTRAnnotator::left_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt)
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
void VNTRAnnotator::right_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt)
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
ReferenceRegion VNTRAnnotator::extract_regions_by_fuzzy_alignment(bcf_hdr_t* h, bcf1_t* v)
{
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

    ReferenceRegion region = ReferenceRegion(min_beg1, max_end1, seq);

    if (seq_len) free(seq);

    return region;
}

/**
 * Fuzzy left align alleles allowing for mismatches and indels defined by penalty.
 */
uint32_t VNTRAnnotator::fuzzy_left_align(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
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
 */
uint32_t VNTRAnnotator::fuzzy_right_align(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty)
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
 * Detect allele lower bound extent.
 */
void VNTRAnnotator::detect_lower_bound_allele_extent(const char* chrom, int32_t& pos1, std::vector<std::string>& alleles, int32_t& start1, int32_t& end1)
{
    if (alleles.size()==2)
    {
        std::string ref = alleles[0];
        std::string alt = alleles[1];

        trim(pos1, ref, alt);
    }
    else
    {
    }
}

/**
 * Detect candidate flanks given a motif fit.
 * Update model atttributes.
 */
void VNTRAnnotator::search_flanks(const char* chrom, int32_t start1, char* motif)
{
    //given chromosome position and motif
    //attempt to

    ////////////////////
    //detect right flank
    ////////////////////

    //lfhmm
    int32_t lflank_len, ref_len;
    char* lflank, *ru, *rflank, *ref;
    int32_t ru_len;
    std::string qual(2048, 'K');

    //make left flank
    lflank = faidx_fetch_seq(fai, chrom, start1-10, start1-1, &lflank_len);
    ref = faidx_fetch_seq(fai, chrom, start1-10, start1+100, &ref_len);
    ++ru;

    lfhmm->set_model(lflank, motif);
    lfhmm->set_mismatch_penalty(5);
    lfhmm->align(ref, qual.c_str());
    lfhmm->print_alignment();

    //get genome position of right flank

    //extract rflank
    //pick flank
    lflank = faidx_fetch_uc_seq(fai, chrom, start1-10, start1-1, &lflank_len);

    //rfhmm
    rfhmm->set_model(motif, rflank);
    rfhmm->set_mismatch_penalty(5);
    rfhmm->align(ref, qual.c_str());
    rfhmm->print_alignment();


    //////////////////////////////
    //annotate STR characteristics
    //////////////////////////////

}

/**
 * Extracts the shortest repeat unit in a sequence.
 */
char* VNTRAnnotator::get_shortest_repeat_motif(char* allele, int32_t len)
{
    std::cerr << "get shortest repeatmotif " << allele << " : " << len << "\n";

    size_t i = 0;
    size_t sub_motif_len;
    while ((sub_motif_len=factors[len][i])!=len)
    {
        std::cerr << "sub motif len : " << sub_motif_len << " " <<  i << "\n";

        bool exact = true;

        size_t n_sub_motif = len/sub_motif_len;

        for (size_t i=0; i<sub_motif_len; ++i)
        {
            char b;
            for (size_t j=0; j<n_sub_motif; ++j)
            {
                if (j)
                {
                    if (b != allele[j*sub_motif_len+i])
                    {
                        exact = false;
                        break;
                    }
                }
                else
                {
                    b = allele[j*sub_motif_len+i];
                }
            }

            if (!exact) break;
        }

        if (exact) break;
        ++i;
    }

    char *motif = allele+len-sub_motif_len;

    return motif;
};

/**
 * Gets motif of a repeat unit.
 */
std::string VNTRAnnotator::get_motif(std::string& ru)
{
    std::string motif = "";
    for (size_t i=0; i<ru.size(); ++i)
    {
        std::string phase = shift_phase(ru, i);
        std::string rc = reverse_complement(phase);
        motif = phase < rc ? phase : rc;
    }

    return motif;
}

/**
 * Reverse complement a sequence.
 */
std::string VNTRAnnotator::reverse_complement(std::string& seq)
{
    std::string rc = "";

    for (size_t i=seq.size()-1; i>0; --i)
    {
        char b = seq.at(i);

        switch (b)
        {
            case 'A':
                rc.append(1, 'A');
                break;
            case 'C':
                rc.append(1, 'C');
                break;
            case 'G':
                rc.append(1, 'G');
                break;
            case 'T':
                rc.append(1, 'T');
                break;
        }
    }

    return rc;
}

/**
 * Shifts a sequence to the right by i bases.
 */
std::string VNTRAnnotator::shift_phase(std::string& seq, size_t i)
{
    i = i<seq.size() ? i : i%seq.size();
    std::string shifted = seq.substr(i, seq.size()-i);
    shifted.append(seq, 0, i);

    return shifted;
}
