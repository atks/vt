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

#ifndef STR_H
#define STR_H

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

KHASH_MAP_INIT_STR(mdict, int32_t);

/**
 * STR Motifs
 */
class STRMotif
{
    public:
    
    //factors[n][index]
    int32_t** factors;
    
    khash_t(mdict) *motifs;
    
    /**
     * Constructor.
     */
    STRMotif()
    {
        motifs = kh_init(mdict);
        
        factors = new int32_t*[31];
//        factors[0] = {0};
//        factors[1] = {1};
//        factors[2] = {1,2};
//        factors[3] = {1,3};
//        factors[4] = {1,2,4};
//        factors[5] = {1,5};
//        factors[6] = {1,2,3,6};
//        factors[7] = {1,7};
//        factors[8] = {1,2,4,8};
//        factors[9] = {1,3,9};
//        factors[10] = {1,2,5,10};
//        factors[11] = {1,11};
//        factors[12] = {1,2,3,4,6,12};
//        factors[13] = {1,13};
//        factors[14] = {1,2,7,14};
//        factors[15] = {1,3,5,15};
//        factors[16] = {1,2,4,8,16};
//        factors[17] = {1,17};
//        factors[18] = {1,2,6,9,18};
//        factors[19] = {1,19};
//        factors[20] = {1,2,4,5,10,20};
//        factors[21] = {1,3,7,21};
//        factors[22] = {1,2,11,22};
//        factors[23] = {1,23};
//        factors[24] = {1,2,3,4,6,8,12,24};
//        factors[25] = {1,5,25};
//        factors[26] = {1,2,13,26};
//        factors[27] = {1,27};
//        factors[28] = {1,2,4,7,14,28};
//        factors[29] = {1,29};
//        factors[30] = {1,2,3,5,6,10,15,30};
    };

    /**
     * Destructor.
     */
    ~STRMotif();

    /**
     * Takes in a set of alleles and suggests a repeat unit.
     */
    char** infer_motif(char** alleles, int32_t n_allele);

    /**
     * Extracts the shortest repeat unit in a sequence.
     */
    char* get_shortest_repeat_motif(char* allele, int32_t len);

    /**
     * Prints variant information.
     */
    void print();
    
};

#endif