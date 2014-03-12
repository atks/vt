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

#ifndef LHMM_GENOTYPING_RECORD_H
#define LHMM_GENOTYPING_RECORD_H

#include <iostream>  
#include <fstream>   
#include <iomanip>   
#include <string>
#include <regex.h>
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "htslib/faidx.h"
#include "utils.h"
#include "lhmm.h"
#include "program.h"
#include "hts_utils.h"
#include "bam_ordered_reader.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "variant_manip.h"
#include "genotyping_record.h"

#include "lhmm_genotyping_record.h"

typedef struct
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

/**
 * Maintains read information and allows for additional reads
 * till VCF record can be printed out.
 */
class LHMMGenotypingRecord : public GenotypingRecord
{
    public:
    bcf1_t *v;
    int32_t vtype;
    uint32_t read_no;

    khash_t(rdict) *reads;
    khiter_t k;
    int32_t ret;

    //probes for probing indels
    char* ref_probe;
    char* alt_probe;
    int32_t plen;

    kstring_t rqs;
    kstring_t aqs;

    LHMMGenotypingRecord(faidx_t *fai=NULL);

    LHMMGenotypingRecord(bcf_hdr_t *h, bcf1_t *v, faidx_t *fai=NULL);
    
    ~LHMMGenotypingRecord();
    
    /**
     * Initializes a candidate variant for genotyping.
     */
    bool initialize(bcf_hdr_t *h, bcf1_t *v, faidx_t *fai);
    
    /**
     * Initializes a candidate VCF record. Returns false if failure.
     */
    bool set(bcf1_t *v);
   
    /**
     * Prints records
     */
    void genotype(bam1_t *b);
    
    void genotype_indel(bam1_t* s);
    
    /**
     * Prints records
     */
    void print(BCFOrderedWriter *odw);
    
    /**
     * Prints records
     */
    void clear();
    
    /**
     * Generates a probing haplotype with flanks around the variant of interest.
     * Flanks are equal length
     */
    void generate_probes(const char* chrom,
                            int32_t pos1, 
                            uint32_t probeDiff, //
                            std::vector<std::string>& alleles, //store alleles
                            std::vector<std::string>& probes, //store probes
                            uint32_t min_flank_length,
                            int32_t& preambleLength); //store preamble length
    
    /**
     * Iteratively called function for generating a haplotype.
     */
    void generate_probes(const char* chrom,
                            int32_t pos1,
                            uint32_t flankLength,
                            uint32_t& currentDiff,
                            uint32_t& length,
                            uint32_t gald,
                            std::vector<uint32_t>& diff,
                            std::vector<std::string>& alleles,
                            std::vector<std::string>& probes);
};

#endif