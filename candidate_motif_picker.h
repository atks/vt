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

#ifndef CANDIDATE_MOTIF_PICKER_H
#define CANDIDATE_MOTIF_PICKER_H

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
#include "flank_detector.h"

//forms of choosing a motif
#define PICK_BEST_MOTIF      0

/**
 * Class for determining basic traits of an indel
 * motifs, flanks and VNTR type statistics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
class CandidateMotifPicker
{
    public:

    //model
    bcf1_t* v;    
    uint32_t max_mlen; //maximum length for motif search in the fast tree.

    bool debug;
    int32_t max_len;

    ///////
    //tools
    ///////
    MotifTree* mt;

    /**
     * Constructor.
     */
    CandidateMotifPicker(bool debug=false);

    /**
     * Destructor.
     */
    ~CandidateMotifPicker();

    /**
     * Initialize motif tree and generate a pool candidate motifs.
     */
    void generate_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, Variant& variant);

    /**
     * Initialize candidate motif from VCF record.
     */
    void set_motif_from_info_field(Variant& variant);
    
    /**
     * Choose the next best motif.
     */
    bool next_motif(bcf_hdr_t* h, bcf1_t* v, Variant& variant);
    
    /**
     * Checks if motif is in indel fragment.
     */
    bool is_in_indel_fragment(std::string motif);
 
    /**
     * Chooses a phase of the motif that is appropriate for the alignment
     */
    std::string choose_repeat_unit(std::string& ref, std::string& motif);
};

#endif