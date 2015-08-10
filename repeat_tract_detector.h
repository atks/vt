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

#ifndef REPEAT_TRACT_DETECTOR_H
#define REPEAT_TRACT_DETECTOR_H

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
class RepeatTractDetector
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
    CandidateRegionExtractor* cre;
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
    RepeatTractDetector(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~RepeatTractDetector();

    /**
     * Infer flanks  motif discovery.
     *
     * returns
     * a. motif concordance
     * b. purity concordance
     * c. left flank
     * d. right flank
     */
    void infer_flanks(bcf_hdr_t* h, bcf1_t* v, std::string& motif);

    /**
     * Detect repeat region.
     */
    void detect_repeat_region(bcf_hdr_t* h, bcf1_t *v, Variant& variant, uint32_t mode);

    /**
     * Detect candidate flanks given a motif fit.
     * Update model atttributes.
     */
    void search_flanks(const char* chrom, int32_t start1, char* motif);

};

#endif