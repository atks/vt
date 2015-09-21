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

#ifndef GENOTYPING_RECORD_H
#define GENOTYPING_RECORD_H

#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "bcf_ordered_writer.h"
#include "variant.h"

/**
 * A generic record that holds information for genotyping a
 * variant across multiple samples.
 *
 * Maintains read information and allows for additional reads
 * till VCF record can be printed out.
 */
class GenotypingRecord
{
    public:
    bcf_hdr_t *h;
    bcf1_t *v;
    int32_t rid;
    int32_t pos1; //position of variant
    //[beg1,end1] is the required overlapping of the variant against the aligned read necessary to make a genotype call.
    //for SNPs, beg1=end1=pos1
    //
    //for Indels, this refers to the flanking positions
    //   insertion 
    //        if T/TG - beg1=pos1, end1=pos1+1
    //        if T/GT - beg1=pos1-1, end1=pos1
    //   deletion  
    //        if TG/T - beg1=pos1, end1=pos1+length(REF)
    //        if TG/G - beg1=pos1-1, end1=pos1+length(REF)-1
    int32_t beg1;
    int32_t end1;
    
    int32_t vtype;

    //indel specific record
    int32_t dlen;
    uint32_t len;
    std::string indel;

    //vntr specific record
    std::string motif;

    //vntr specific record
    std::vector<float> counts;

    //for records that observe at least one alternate observation
    std::vector<uint32_t> quals;  //for SNPs, store BQ, for Indels, store AQ
    std::vector<uint32_t> map_quals;
    std::string strands;
    std::vector<int32_t> alleles;
    std::vector<uint32_t> cycles;
    std::vector<uint32_t> no_mismatches;

    //for records that only have reference observation
    uint32_t no_nonref;
    std::vector<uint32_t> allele_depth_fwd;
    std::vector<uint32_t> allele_depth_rev;
    uint32_t depth, depth_fwd, depth_rev;
    uint32_t base_qualities_sum;

    /**
     * Constructor.
     * @v - VCF record.
     */
    GenotypingRecord(bcf_hdr_t*h, bcf1_t *v, int32_t vtype);

    /**
     * Clears this record.
     */
    void clear();

    /**
     * Destructor.
     */
    ~GenotypingRecord();
};

#endif