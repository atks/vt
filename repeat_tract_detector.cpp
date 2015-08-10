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

#include "repeat_tract_detector.h"

/**
 * Constructor.
 */
RepeatTractDetector::RepeatTractDetector(std::string& ref_fasta_file, bool debug)
{
    vm = new VariantManip(ref_fasta_file.c_str());

    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }

    cre = new CandidateRegionExtractor(ref_fasta_file, debug);

    max_mlen = 10;
    mt = new MotifTree(max_mlen, debug);

    this->debug = debug;
    qual.assign(256, 'K');

    ahmm = new AHMM();
    rfhmm = new RFHMM();
    lfhmm = new LFHMM();
};

/**
 * Destructor.
 */
RepeatTractDetector::~RepeatTractDetector()
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
 * Infer flanks  motif discovery.
 *
 * returns
 * a. motif concordance
 * b. purity concordance
 * c. left flank
 * d. right flank
 */
void RepeatTractDetector::infer_flanks(bcf_hdr_t* h, bcf1_t* v, std::string& motif)
{
    //fit into flanks

    //left alignment

    //right alignment

//        int32_t no_candidate_motifs;
//
//        const char* chrom = bcf_get_chrom(h, v);
//
//        int32_t pos1 = bcf_get_pos1(v);

//        std::vector<CandidateMotif> candidate_motifs;
//        pick_candidate_motifs(h, v, candidate_motifs);
//
//        std::string eru = candidate_motifs[0].motif;
//        std::string emotif = get_motif(eru);
//
//        std::cerr << "exact motif : " << emotif << "\n";
//        std::cerr << "exact ru    : " << eru << "\n";
//
//        vntr.eru = eru;
//        vntr.emotif = emotif;
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
//        vntr.iregion.beg1 = motif_beg1;
//        vntr.iregion.end1 = motif_end1;
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
 * Pick candidate motifs.
 * candidate_motifs contain motifs and a measure of confidence
 */
//void RepeatTractDetector::pick_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, std::vector<CandidateMotif>& candidate_motifs)
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
 * Detect candidate flanks given a motif fit.
 * Update model atttributes.
 */
void RepeatTractDetector::search_flanks(const char* chrom, int32_t start1, char* motif)
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