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
    
    //factors[n][index], for determining what sub repeat units to examine
    int32_t** factors;
    
    khash_t(mdict) *motifs;
    
    /**
     * Constructor.
     */
    STRMotif()
    {
        motifs = kh_init(mdict);
        
        //update factors
        factors = new int32_t*[32];
        for (size_t i=1; i<=32; ++i)
        {
            factors[i] = new int32_t[32];
            int32_t count = 0;
            for (size_t j=1; j<=i/2; ++j)
            {                
                if ((i%j)==0)
                {
                    factors[i][count++] = j;
                }    
            }   
        }        
    };

    /**
     * Destructor.
     */
    ~STRMotif();

    /**
     * Suggests a set of repeat motif candidates in a set of alleles.
     */
    char** suggest_motifs(char** alleles, int32_t n_allele, int32_t &no_candidate_motifs);

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