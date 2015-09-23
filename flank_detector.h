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
#define LAI_2003_STR       1
#define KELKAR_2008_STR    2
#define FONDON_2012_STR    3
#define ANANDA_2013_STR    4
#define WILLEMS_2014_STR   5
#define TAN_KANG_2015_VNTR 6


//forms of alignment
#define REFERENCE                                0
#define EXACT_LEFT_RIGHT_ALIGNMENT               1
#define FUZZY_LEFT_RIGHT_ALIGNMENT               2
#define FUZZY_LEFT_RIGHT_ALIGNMENT_WITH_PENALTY  3

#define CLIP_ENDS 0
#define CLIP_1L2R 1
#define FRAHMM    2

/**
 * Class for determining basic traits of an indel
 * motifs, flanks and VNTR type statistics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
class FlankDetector
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
    bool debug;
    int32_t max_len;

    ////////
    //raHMMs
    ////////

    ///////
    //tools
    ///////
    faidx_t* fai;
    
    std::string qual;
    AHMM* ahmm;
    LFHMM* lfhmm;
    RFHMM* rfhmm;    
    

    //for retrieving sequences
    int8_t* seq;

    /**
     * Constructor.
     */
    FlankDetector(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~FlankDetector();

    /**
     * Detect repeat region.
     */
    void detect_flanks(bcf_hdr_t* h, bcf1_t *v, Variant& variant, uint32_t mode);

    /**
     * Shifts a string.
     */
    std::string shift_str(std::string& seq, uint32_t i);
        
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
};

#endif