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

#include "rfhmm.h"
#include "lfhmm.h"
#include "ahmm.h"
#include "variant.h"
#include "motif_tree.h"

//modes for picking motifs
#define NO_REQUIREMENT                 0
#define CHECK_MOTIF_PRESENCE_IN_ALLELE 1

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
    void generate_candidate_motifs(Variant& variant);
    
    /**
     * Generate candidate motids from a repeat tract.
     *
     * 1. assigns the repeat_tract to the fuzzy_repeat_tract field of the variants VNTR.
     * 2. examines this sequence for candidate motifs that are stored in the priority queue in the MotifTree class.
     * 3. the candidate motifs are then to be accessed via next_motif()
     */
    void generate_candidate_motifs(char* repeat_tract, Variant& variant);

    /**
     * Initialize candidate motif from VCF record.
     */
    void set_motif_from_info_field(Variant& variant);
  
    /**
     * Updates the motif of an indel allele.
     *
     * Returns true if the motif is from a simple indel.
     * Returns false if the motif is ambiguous.
     */
    void update_exact_repeat_unit(Variant& variant);
          
    /**
     * Choose the next best motif.
     */
    bool next_motif(Variant& variant, int32_t mode=NO_REQUIREMENT);
    
    /**
     * Gets inserted or deleted allele of a biallelic indel.
     */
    bool get_indel(std::string ref, std::string alt, std::string& indel);
    
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