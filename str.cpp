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

#include "str.h"

/**
 * Takes in a set of alleles and suggest a set of candidate repeat units.
 */
char** STRMotif::infer_motif(char** alleles, int32_t n_allele)
{
    char* ref = alleles[0];
    int32_t ref_len = strlen(ref);
    char *alt, *motif;
    int32_t alt_len;

    //grab all candidate alleles
    for (size_t i=1; i<n_allele; ++i)
    {
        alt = alleles[i];
        alt_len = strlen(alleles[i]);

        //get length difference
        int32_t dlen = alt_len-ref_len;
        //extract fragment
        if (dlen>0)
        {
            motif = &alt[alt_len-dlen];
        }
        else
        {
            motif = &ref[ref_len+dlen];
        }

        int32_t len = abs(dlen);
        size_t j = 0;
        char* m = get_shortest_repeat_motif(alleles[i], len);

        int32_t ret;
        khiter_t k;
        if (kh_get(mdict, motifs, m)==kh_end(motifs))
        {
            k = kh_put(mdict, motifs, m, &ret);
            kh_value(motifs, k) = 1;
        }
        else
        {
            kh_value(motifs, k) += 1;
        }
    }
    
    int32_t no_candidate_motifs = kh_size(motifs);
    
    char** candidate_motifs = (char**) malloc(no_candidate_motifs*sizeof(char*)); 
    
    khiter_t k;
    int32_t i = 0;
    for (k=kh_begin(motifs); k!=kh_end(motifs); ++k)
    {
        if (kh_exist(motifs, k))
        {
            candidate_motifs[i] = (char*) kh_key(motifs, k);
        }
    }
    kh_clear(mdict, motifs);
}

/**
 * Extracts the shortest repeat unit in a sequence.
 */
char* STRMotif::get_shortest_repeat_motif(char* allele, int32_t len)
{
    size_t i = 0;
    size_t sub_motif_len;
    while ((sub_motif_len=factors[len][i])!=len)
    {
        bool exact = true;

        size_t n_sub_motif = len/sub_motif_len;

        for (size_t i=0; i<sub_motif_len; ++i)
        {
            char b;            
            for (size_t j=0; j<n_sub_motif; ++j)
            {
                if (j)
                {
                    if (b != allele[j*sub_motif_len+i])
                    {
                        exact = false;
                        break;
                    }
                }
                else
                {
                    b = allele[j*sub_motif_len+i];
                }     
            }        
            
            if (!exact) break;
        }    
        
        if (exact) break;
    }

    char *motif = allele+len-sub_motif_len;
    
    return motif;
};