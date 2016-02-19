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

#include "candidate_region_extractor.h"
#include "candidate_motif_picker.h"
#include "flank_detector.h"

//definition of STRs
#define LAI_2003_STR       1
#define KELKAR_2008_STR    2
#define FONDON_2012_STR    3
#define ANANDA_2013_STR    4
#define WILLEMS_2014_STR   5
#define TAN_KANG_2015_VNTR 6

//VNTR annotation modes
#define REPEAT_TRACT_FEATURES  0 // compute composition, entropy, score, trf_score, 

/**
 * Class for determining basic traits of an indel
 * motifs, flanks and VNTR type statistics.
 * MOTIF
 * RU
 * RL
 * FLANKS
 * FZ_FLANKS
 * FZ_CONCORDANCE
 */
class VNTRAnnotator
{
    public:

    /////////
    //options
    /////////
    bool debug;
    int32_t max_len;

    ///////
    //tools
    ///////
    VariantManip *vm;
    faidx_t* fai;
    CandidateRegionExtractor* cre;
    CandidateMotifPicker* cmp;
    FlankDetector* fd;

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
    void annotate(Variant& variant, int32_t amode);

    /**
     * Returns true if is to be classified as a VNTR
     */
    bool is_vntr(Variant& variant, int32_t mode, std::string& method);
};

#endif