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
#include <cstring>
#include <iostream>
#include <vector>
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "allele.h"
#include "hts_utils.h"

#ifndef VARIANT_H
#define VARIANT_H

class Region
{
    public:

    int32_t beg1;
    int32_t end1;
};

/**
 * Variant.
 */
class Variant
{
    public:

    int32_t type;       //aggegrated type from the alleles
    
    //sum from all the alleles
    int32_t ts;         //no. of transitions
    int32_t tv;         //no. of tranversions (mlen-ts)
    int32_t ins;        //no. of insertions
    int32_t del;        //no. of deletions

    std::string motif;   //motif of VNTR
    std::string ru;      //repeat unit of VNTR
    int32_t mlen;        //length of motif
    float rcn;           //reference copy number
    int32_t rlen;        //reference length of repeat tract in bases
        
    std::string emotif;
    std::string eru;
    Region eregion;

    std::string imotif;
    std::string iru;
    Region iregion;

    std::vector<Allele> alleles;

    /**
     * Constructor.
     */
    Variant();

    /**
     * Destructor.
     */
    ~Variant();

    /**
     * Prints variant information.
     */
    void print();

    /**
     * Returns true if variant contains an allele that is potentially frame shifting.
     */
    bool exists_frame_shift();

    /**
     * Converts VTYPE to string.
     */
    std::string vtype2string(int32_t VTYPE);

    /**
     * Clears variant information.
     */
    void clear();
};

#endif
