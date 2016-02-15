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

#ifndef VNTR_H
#define VNTR_H

#include "utils.h"

//VNTR update modes
#define FINAL  1 // update final attributes
#define EXACT  2 // update exact attributes
#define FUZZY  4 // update fuzzy attributes

/**
 * Class for representing a VNTR.
 *
 * Includes functions specific to tandem repeats.
 */
class VNTR
{
    public:

    //chromosome
    int32_t rid;  //rid, redundant data with Variant. todo: something about this.

    ///////
    //motif
    ///////
    std::string motif;         //motif
    std::string ru;            //repeat unit on the reference
    int32_t mlen;              //length of motif
    std::string basis;         //unique bases found in motif
    int32_t blen;              //basis length

    std::string definition_support; //either exact or fuzzy

    //////////////
    //repeat tract
    //////////////
    std::string repeat_tract;   //repeat tract
    int32_t comp[4];            //composition of bases in repeat tract
    float entropy ;             //sequence entropy of repeat tract
    float entropy2;            //dinucleotide sequence entropy of repeat tract
    int32_t beg1;               //beginning of repeat tract
    int32_t end1;               //end of repeat tract
    int32_t rl;                 //length of repeat tract in base pairs
    int32_t ll;                 //length of repeat tract (including longest alternate allele) in base pairs
    int32_t no_exact_ru;              //number exact repeat units from hmm
    int32_t total_no_ru;        //total no of repeat units from hmm
    float score;                //motif concordance from hmm
    int32_t trf_score;          //TRF score of exact repeat tract
    std::string lflank;         //left flank
    std::string rflank;         //right flank

    ////////////////////
    //exact repeat tract
    ////////////////////
    std::string exact_repeat_tract;   //repeat tract
    int32_t exact_comp[4];            //composition of bases in repeat tract
    float exact_entropy ;             //sequence entropy of repeat tract
    float exact_entropy2 ;            //dinucleotide sequence entropy of repeat tract
    int32_t exact_beg1;               //beginning of repeat tract
    int32_t exact_end1;               //end of repeat tract
    int32_t exact_rl;                 //length of repeat tract in base pairs
    int32_t exact_ll;                 //length of repeat tract (including longest alternate allele) in base pairs
    int32_t exact_no_exact_ru;        //number exact repeat units from hmm
    int32_t exact_total_no_ru;        //total no of repeat units from hmm
    float exact_score;                //motif concordance from hmm
    int32_t exact_trf_score;          //TRF score of exact repeat tract
    std::string exact_lflank;         //left flank
    std::string exact_rflank;         //right flank

    /////////////////
    //fuzzy alignment
    /////////////////
    std::string fuzzy_repeat_tract;   //repeat tract
    int32_t fuzzy_comp[4];            //composition of bases in repeat tract
    float fuzzy_entropy ;             //sequence entropy of repeat tract
    float fuzzy_entropy2;             //dinucleotide sequence entropy of repeat tract
    int32_t fuzzy_beg1;               //beginning of repeat tract
    int32_t fuzzy_end1;               //end of repeat tract
    int32_t fuzzy_rl;                 //length of repeat tract in base pairs
    int32_t fuzzy_ll;                 //length of repeat tract (including longest alternate allele) in base pairs
    int32_t fuzzy_no_exact_ru;        //number exact repeat units from hmm
    int32_t fuzzy_total_no_ru;        //total no of repeat units from hmm
    float fuzzy_score;                //motif concordance from hmm
    int32_t fuzzy_trf_score;          //TRF score of fuzzy repeat tract
    std::string fuzzy_lflank;         //left flank
    std::string fuzzy_rflank;         //right flank

    //types of repeat
    bool is_large_repeat_tract;
    bool is_interspersed_repeat_tract;

    /**
     * Constructor.
     */
    VNTR();

    /**
     * Clear object.
     */
    void clear();

    /**
     * Checks for equality.
     */
    bool equals(VNTR& vntr);

    /**
     * Return the string of unique bases in a motif.
     */
    static std::string get_basis(std::string& motif);

    /**
     * Return the string of unique bases in a motif.
     */
    static std::string get_basis(char* motif, uint32_t n);

    /**
     * Shifts a string.
     */
    static std::string shift_str(std::string& seq, uint32_t i);

    /**
     * Print object.
     */
    void print();
};
#endif