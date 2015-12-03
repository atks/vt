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
 * This represents a Variant and is augmented on top of VCF's record to handle the notion of variants as defined by us.
 */
class Variant
{
    public:

    //aggegrated type from the alleles
    int32_t type;

    //location information
    std::string chrom;
    uint32_t rid;
    uint32_t pos1; //position of first reference base in VCF record
    uint32_t beg1; //for detecting overlaps, for indels and VNTRs, this will reflect lend1 
    uint32_t end1; //for detecting overlaps, for indels and VNTRs, this will reflect rbeg1

    //linked VCF record
    bcf_hdr_t* h;
    bcf1_t* v;
    
    //associated VCF records for merging into a multiallelic
    std::vector<bcf1_t*> vs;
    std::vector<bcf1_t*> snp_vs;
    std::vector<bcf1_t*> indel_vs;
    std::vector<bcf1_t*> vntr_vs;
    
    //contains alleles
    std::vector<Allele> alleles;
        
    //sum from all the alleles
    int32_t ts;         //no. of transitions
    int32_t tv;         //no. of tranversions (mlen-ts)
    int32_t ins;        //no. of insertions
    int32_t del;        //no. of deletions

    //overlapping statistics
    //for normal variants
    // - the number of other variants overlapping with this 
    //for candidate multiallelic/complex variant
    // - the number of variants considered for merging into a multiallelic/complex variant
    int32_t no_overlapping_snps;
    int32_t no_overlapping_indels;
    int32_t no_overlapping_vntrs;
    
    //describes VNTR
    VNTR vntr;

    /**
     * Constructor.
     */
    Variant(bcf_hdr_t* h, bcf1_t* v);

    /**
     * Constructor.
     */
    Variant(Variant* v1, Variant* v2);

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
     * Classifies variants based on observed alleles in vcf record.
     */
    int32_t classify(bcf_hdr_t *h, bcf1_t *v);

    /**
     * Updates VNTR related information from INFO fields.
     */
    void update_vntr_from_info_fields(bcf_hdr_t *h, bcf1_t *v);

    /**
     * Prints variant information.
     */
    void print();
    
    /**
     * Gets a string representation of the underlying VNTR by exact alignment.
     */
    void get_vntr_string(kstring_t* s);

    /**
     * Gets a string representation of the underlying VNTR by fuzzy alignment.
     */
    void get_fuzzy_vntr_string(kstring_t* s);

    /**
     * Converts VTYPE to string.
     */
    static std::string vtype2string(int32_t VTYPE);
};

#endif
