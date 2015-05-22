/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#include <cstdint>
#include <iostream>
#include "htslib/kstring.h"
#include <vector>

#ifndef ALLELE_H
#define ALLELE_H

#define VT_REF      0   //dlen==0 && diff==0
#define VT_SNP      1   //min(rlen,alen)==1 && diff==1
#define VT_MNP      2   //min(rlen,alen)==diff  || all allele lengths are the same and greater than 1.
#define VT_INDEL    4   //diff!=0 && (rlen==1 || alen==1)
#define VT_CLUMPED  8   //all others sequence explicit
#define VT_VNTR     16  //variable number of tandem repeats
#define VT_SV       32  //structural variant

/**
 * Allele.
 *
 * Allele can be described to be a SNP or INDEL or etc. with respect to the reference
 * type - give the variant type
 */
class Allele
{
    public:

    int32_t type;  //allele type
    int32_t diff;  //number of difference bases when bases are compared
    int32_t alen;  //length(alt)
    int32_t dlen;  //length(alt)-length(ref)
    int32_t mlen;  //min shared length
    int32_t ts;    //no. of transitions
    int32_t tv;    //no. of tranversions (mlen-ts)
    int32_t ins;   //no. of insertions
    int32_t del;   //no. of deletions
    std::string sv_type; //hierarchical descriptor for the imprecise allele type

    /**
     * Constructor for explict sequence variant.
     */
    Allele(int32_t type, int32_t diff, int32_t alen, int32_t dlen, int32_t mlen, int32_t ts, int32_t tv);

    /**
     * Constructor for VNTR.
     */
    Allele(int32_t type);

    /**
     * Constructor for SV.
     */
    Allele(int32_t type, std::string& sv_type);

    /**
     * Destructor.
     */
    ~Allele();

    /**
     * Special dictionary for some reserve types.
     * CN\d+ be CNV and VN\d+ be VNTR
     */
    std::string reduce_sv_type(std::string& sv_type);

    /**
     * Clear variables.
     */
    void clear();

    /**
     * Print allele.
     */
    void print();
};

#endif
