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

#ifndef VNTR_ANNOTATOR_H
#define VNTR_ANNOTATOR_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include "hts_utils.h"
#include "htslib/kstring.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "rfhmm.h"
#include "lfhmm.h"
#include "ahmm.h"
#include "variant_manip.h"
#include "program.h"
#include "motif_tree.h"
#include "vntr.h"
#include "candidate_region_extractor.h"

//definition of STRs
#define LAI_2003_STR       0
#define KELKAR_2008_STR    1
#define FONDON_2012_STR    2
#define ANANDA_2013_STR    3
#define WILLEMS_2014_STR   4
#define TAN_KANG_2015_VNTR 5


//forms of alignment
#define REFERENCE                                0
#define EXACT_LEFT_RIGHT_ALIGNMENT               1
#define FUZZY_LEFT_RIGHT_ALIGNMENT               2
#define FUZZY_LEFT_RIGHT_ALIGNMENT_WITH_PENALTY  3

//forms of choosing a motif
#define PICK_BEST_MOTIF             0

#define ALLELE_EXACT  1
#define ALLELE_FUZZY  2

#define CLIP_ENDS 0
#define CLIP_1L2R 1

/**
 * Class for determining basic traits of an indel
 * motifs, flanks and VNTR type statistics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
class VNTRAnnotator
{
    public:
//##INFO=<ID=VT_RU,Number=1,Type=String,Description=\"Repeat unit in a STR or Homopolymer\">");
//##INFO=<ID=VT_RL,Number=1,Type=Integer,Description=\"Repeat Length\">");
//##INFO=<ID=VT_LFLANK,Number=1,Type=String,Description=\"Right Flank Sequence\">");
//##INFO=<ID=VT_RFLANK,Number=1,Type=String,Description=\"Left Flank Sequence\">");
//##INFO=<ID=VT_LFLANKPOS,Number=2,Type=Integer,Description=\"Positions of left flank\">");
//##INFO=<ID=VT_RFLANKPOS,Number=2,Type=Integer,Description=\"Positions of right flank\">");
//##INFO=<ID=VT_MOTIF_DISCORDANCE,Number=1,Type=String,Description=\"Descriptive Discordance for each reference repeat unit.\">");
//##INFO=<ID=VT_STR_CONCORDANCE,Number=1,Type=Float,Description=\"Overall discordance of RUs.\">");

    uint32_t max_mlen; //maximum length for motif search in the fast tree.

    //model
    char* motif;
    int32_t motif_len;
    int32_t ref_len;
    char* lflank;
    char* rflank;
    bool exact;
    int32_t* motif_concordance;
    float* motif_completeness;
    float concordance;
    std::vector<CandidateMotif> candidate_motifs;
    bool debug;
    int32_t max_len;

    AHMM* ahmm;
    std::string qual;

    std::string MOTIF;
    std::string RU;
    std::string RL;
    std::string REF;
    std::string REFPOS;
    std::string SCORE;
    std::string TR;

    ///////
    //tools
    ///////
    VariantManip *vm;
    faidx_t* fai;
    CandidateRegionExtractor* cre;
    MotifTree* mt;

    //for retrieving sequences
    int8_t* seq;

    //factors[n][index], for determining what sub repeat units to examine
    int32_t** factors;

    /**
     * Constructor.
     */
    VNTRAnnotator(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~VNTRAnnotator();

    /**
     * Annotates VNTR characteristics.
     * @mode
     *   e - determine by exact alignment
     *   f - determine by fuzzy alignment
     *   p - determine by penalized fuzzy alignment
     *   h - using HMMs     
     *   x - integrated models     
     */
    void annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant, std::string mode);

    /**
     * Pick candidate motifs.
     * candidate_motifs contain motifs and a measure of confidence
     */
    void pick_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr);

    /**
     * Chooses a phase of the motif that is appropriate for the alignment
     */
    void choose_best_motif(bcf_hdr_t* h, bcf1_t* v, MotifTree* mt, VNTR& vntr, uint32_t mode);

    /**
     * Pick shortest motif.
     */
    std::string pick_motif(std::string& sequence);

    /**
     * This is a quick scan for a motif that is exactly repeated.
     */
    std::string scan_exact_motif(std::string& sequence);

    /**
     * Detect repeat region.
     */
    void detect_repeat_region(bcf_hdr_t* h, bcf1_t *v, Variant& variant, uint32_t mode);

    /**
     * Chooses a phase of the motif that is appropriate for the alignment
     */
    std::string choose_repeat_unit(std::string& ref, std::string& motif);

    /**
     * Trim alleles.
     */
    void trim(int32_t& pos1, std::string& ref, std::string& alt);

    /**
     * Checks if a vntr is a homopolymer.
     */
    bool is_homopolymer(bcf_hdr_t* h, bcf1_t* v);

    /**
     * Detect allele lower bound extent.
     */
    void detect_lower_bound_allele_extent(const char* chrom, int32_t& pos1, std::vector<std::string>& alleles, int32_t& start1, int32_t& end1);

    /**
     * Extracts the shortest repeat unit in a sequence.
     */
    char* get_shortest_repeat_motif(char* allele, int32_t len);

    /**
     * Gets motif of a repeat unit.
     */
    std::string get_motif(std::string& ru);

    /**
     * Reverse complement a sequence.
     */
    std::string reverse_complement(std::string& seq);

    /**
     * Shifts a sequence to the right by i bases.
     */
    std::string shift_phase(std::string& seq, size_t i);

    /**
     * Prints vntr information.
     */
    void print();

    /**
     * Returns true if is to be classified as a VNTR
     */
    bool is_vntr(Variant& variant, int32_t mode);

};

#endif