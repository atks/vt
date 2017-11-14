/* The MIT License

   Copyright (c) 2017 Adrian Tan <atks@umich.edu>

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

#ifndef INDEL_ANNOTATOR_H
#define INDEL_ANNOTATOR_H

#include "candidate_region_extractor.h"
#include "candidate_motif_picker.h"
#include "flank_detector.h"

#define START                          0
#define EXACT_LEFT_AND_RIGHT_ALIGNMENT 1
#define ANNOTATE_EXACT                 1
#define PICK_CANDIDATE_REPEAT_UNITS    1
#define FUZZY_LEFT_AND_RIGHT_ALIGNMENT 2
#define EVALUATION                     3
#define CHOOSE_FINAL_ANNOTATION        4
#define ANNOTATE_FUZZY                 7
#define END                            8

/**
 * Class for determining basic traits of an indel
 * motifs, flanks and VNTR type statistics.
 * MOTIF
 * RU
 * RL
 * FLANKS
 * FZ_FLANKS
 * FZ_CONCORDANCE
 *
 * This is intended to replace VNTRAnnotator by implementing the
 * annotation process as a finite state machine.
 *
1. START

2. EXACT_LEFT_AND_RIGHT_ALIGNMENT

3. ANNOTATE_EXACT

4. PICK_CANDIDATE_REPEAT_UNITS

5. FUZZY_LEFT_AND_RIGHT_ALIGNMENT

6. EVALUATION
   if fail requirements, go back to PICK_CANDIDATE_REPEAT_UNITS
        a. fuzzy should be inclusive of exact
        b. when fuzzy does not overlap exact
           - retry with another repeat unit
      redo till no acceptable solution
   last choice is to use EX_RU

7. CHOOSE_FINAL_ANNOTATION

8. END

ALGOTRACE tags
ALGO_TRIES - number of iterations
ALGO_EXACT - final answer is exact==fuzzy
ALGO_RU_TRIES - repeat unit tries

other useful tags
EXNEFZRU - fuzzy repeat unit not the same as exact repat unit

requirements: a. fuzzy should be inclusive of exact
              b. when fuzzy does not overlap exact
                  - retry with another repeat unit

#other useful annotations = relationship between EX_RU and FZ_RU

#consider consolidation program as a form of evaluation.

 */
class IndelAnnotator
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

    //finite state machine
    int32_t state;


    /**
     * Constructor.
     */
    IndelAnnotator(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~IndelAnnotator();

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
};

#endif