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

#ifndef CANDIDATE_REGION_EXTRACTOR_H
#define CANDIDATE_REGION_EXTRACTOR_H

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

//forms of alignment
#define REFERENCE                                0
#define EXACT_LEFT_RIGHT_ALIGNMENT               1
#define FUZZY_LEFT_RIGHT_ALIGNMENT               2
#define FUZZY_LEFT_RIGHT_ALIGNMENT_WITH_PENALTY  3

//forms of choosing a motif
#define PICK_BEST_MOTIF             0

//
#define ALLELE_EXACT  1
#define ALLELE_FUZZY  2

#define CLIP_ENDS 0
#define CLIP_1L2R 1

/**
 * Class for determining basic traits of an indel
 * motifs, flanks and VNTR type statistics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
class CandidateRegionExtractor
{
    public:

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

    ///////
    //tools
    ///////
    VariantManip *vm;
    faidx_t* fai;
    MotifTree* mt;
    RFHMM* rfhmm;
    LFHMM* lfhmm;

    //for retrieving sequences
    int8_t* seq;

    //factors[n][index], for determining what sub repeat units to examine
    int32_t** factors;

    /**
     * Constructor.
     */
    CandidateRegionExtractor(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~CandidateRegionExtractor();

    /**
     * Annotates VNTR characteristics.
     * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
     */
    void annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant, std::string mode);

    /**
     * Pick candidate region.
     *
     * @mode - REFERENCE     use refence field
     *       - ALLELE_EXACT  by exact alignment
     *       - ALLELE_FUZZY  by fuzzy alignment
     */
    void pick_candidate_region(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr, uint32_t mode);

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
     * Extract region to for motif discovery.
     */
    void extract_regions_by_exact_alignment(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr);

    /**
     * Left align alleles.
     */
    void left_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt);

    /**
     * Right align alleles.
     */
    void right_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt);

    /**
     * Extract reference sequence region for motif discovery in a fuzzy fashion.
     */
    void extract_regions_by_fuzzy_alignment(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr);

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
    uint32_t fuzzy_left_align(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty);

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
    uint32_t fuzzy_right_align(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty);

    /**
     * Extract reference sequence region for motif discovery in a fuzzy fashion.
     */
    void extract_regions_by_fuzzy_alignment_with_penalty(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr);

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
    uint32_t fuzzy_left_align_with_penalty(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty);

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
    uint32_t fuzzy_right_align_with_penalty(const char* chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty);

};

#endif