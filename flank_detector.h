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

#ifndef FLANK_DETECTOR_H
#define FLANK_DETECTOR_H

#include "hts_utils.h"
#include "utils.h"
#include "rfhmm.h"
#include "lfhmm.h"
#include "ahmm.h"
#include "variant_manip.h"
#include "reference_sequence.h"
#include "vntr.h"

//modes for flank detection
#define CLIP_ENDS 0   //literally just clip the flanking bases
#define FRAHMM    1   //raHMM alignment
#define POLISH_ENDS 2 //for ensuring that the motif occurs on the ends of the repeat tract

/**
 * Class for determining flanks of an Indel.
 */
class FlankDetector
{
    public:

    bool debug;

    ////////
    //raHMMs
    ////////
    AHMM* ahmm;
    LFHMM* lfhmm;
    RFHMM* rfhmm;
    std::string qual;

    /////////////////////
    //polishing variables
    /////////////////////
    int32_t min_beg0;
    int32_t max_end0;
    std::string polished_repeat_tract;

    /////////////////
    //flank variables
    /////////////////
    int32_t beg1;
    int32_t end1;
    std::string repeat_tract;

    ////////////////////////
    //purity score variables
    ////////////////////////
    std::string ru;
    float score;
    int32_t no_perfect_ru;
    int32_t no_ru;
    float ref;
    int32_t rl;
    int32_t ll;
    int32_t trf_score;

    ///////////////////////////////////
    //composition and entropy variables
    ///////////////////////////////////
    int32_t comp[4];
    float entropy;
    float kl_divergence;
    float entropy2;
    float kl_divergence2;

    ///////
    //tools
    ///////
    ReferenceSequence *rs;

    /**
     * Constructor.
     */
    FlankDetector(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~FlankDetector();

//    /**
//     * Computes purity score of a sequence with respect to a motif.
//     */
//    void compute_purity_score(Variant& variant, int32_t amode);
//
//    /**
//     * Computes purity score of a sequence with respect to a motif.
//     */
//    void detect_flanks(std::string& repeat_tract, std::string& motif);
//        

    /**
     * Detect repeat region.
     */
    void detect_flanks(Variant& variant, uint32_t mode);

    /**
     * Shifts a string.
     */
    std::string shift_str(std::string& seq, uint32_t i);

    /**
     * Score string.
     */
    int32_t compute_score(int32_t start, int32_t len, std::string& a, std::string& b);
        
    /**
     * Chooses a phase of the motif that is appropriate for the alignment from the 5 prime end.
     * This differs from choose_exact_repeat_unit() where the motif is returned
     * if not suitable repeat unit is found.
     */
    std::string choose_5prime_repeat_unit(std::string& seq, std::string& motif);

    /**
     * Chooses a phase of the motif that is appropriate for the alignment from the 5 prime end.
     * If no exact match is available, the best possible match is returned. 
     */
    std::string choose_fuzzy_5prime_repeat_unit(std::string& seq, std::string& motif);
            
    /**
     * Chooses a phase of the motif that is appropriate for the alignment from the 3 prime end.
     * This differs from choose_exact_repeat_unit() where the motif is returned
     * if not suitable repeat unit is found.
     */
    std::string choose_3prime_repeat_unit(std::string& seq, std::string& motif);

    /**
     * Chooses a phase of the motif that is appropriate for the alignment from the 3 prime end.
     * If no exact match is available, the best possible match is returned. 
     */
    std::string choose_fuzzy_3prime_repeat_unit(std::string& seq, std::string& motif);
    
    /**
     * Chooses a phase of the motif that is appropriate for the alignment.
     * This returns the empty string if the motif does not have an exact
     * match in all its phases.
     */
    std::string choose_exact_repeat_unit(std::string& seq, std::string& motif);

    /**
     * Polish repeat tract ends.
     */
    void polish_repeat_tract_ends(Variant& variant);

    /**
     * Polish repeat tract ends.
     */
    void polish_repeat_tract_ends(std::string& repeat_tract, std::string& motif, bool debug=false);

    /**
     * Computes purity score of a sequence with respect to a motif.
     */
    void compute_purity_score(Variant& variant, int32_t amode);

    /**
     * Computes purity score of a sequence with respect to a motif.
     */
    void compute_purity_score(std::string& repeat_tract, std::string& motif);
    
    /**
    * Computes composition and entropy of repeat tract.
    */
    void compute_composition_and_entropy(Variant& variant, int32_t amode);

    /**
     * Computes composition and entropy of repeat tract.
     */
    void compute_composition_and_entropy(std::string& repeat_tract);
};

#endif