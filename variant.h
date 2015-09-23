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
#include "vntr.h"
#include "hts_utils.h"

#ifndef VARIANT_H
#define VARIANT_H

/**
 * Variant.
 */
class Variant
{
    public:

    //aggegrated type from the alleles
    int32_t type;

    //location information
    std::string chrom;
    uint32_t rid;
    uint32_t pos1;
    uint32_t end1;

    //linked VCF file
    bcf1_t* v;

    //contains alleles
    std::vector<Allele> alleles;

    //sum from all the alleles
    int32_t ts;         //no. of transitions
    int32_t tv;         //no. of tranversions (mlen-ts)
    int32_t ins;        //no. of insertions
    int32_t del;        //no. of deletions

    //describes VNTR
    VNTR vntr;

    /**
     * Constructor.
     */
    Variant(bcf1_t* v);

    /**
     * Constructor.
     */
    Variant();

    /**
     * Destructor.
     */
    ~Variant();

    /**
     * Clears variant information.
     */
    void clear();

    /**
     * Prints variant information.
     */
    void print();

    /**
     * Gets a string representation of the underlying VNTR.
     */
    void get_vntr_string(kstring_t* s);

    /**
     * Gets a fuzzy string representation of the underlying VNTR.
     */
    void get_fuzzy_vntr_string(kstring_t* s);

    /**
     * Converts VTYPE to string.
     */
    static std::string vtype2string(int32_t VTYPE);
};

#endif
