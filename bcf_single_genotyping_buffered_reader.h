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

#ifndef BCF_SINGLE_GENOTYPING_BUFFERED_READER_H
#define BCF_SINGLE_GENOTYPING_BUFFERED_READER_H

#include "hts_utils.h"
#include "utils.h"
#include "genotyping_record.h"
#include "snp_genotyping_record.h"
#include "indel_genotyping_record.h"
#include "vntr_genotyping_record.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "variant.h"
#include "variant_manip.h"
#include "log_tool.h"
#include "augmented_bam_record.h"

/**
 * Wrapper for BCFOrderedReader.
 *
 * VCF records are wrapped in GenotyingRecord and are
 * maintained in a buffer.
 *
 */
class BCFSingleGenotypingBufferedReader
{
    public:

    ///////
    //i/o//
    ///////
    BCFOrderedReader* odr;
    BCFOrderedWriter* odw;

    //////////////////
    //buffer related//
    //////////////////
    std::list<GenotypingRecord*> buffer;
    std::string chrom;
    AugmentedBAMRecord as;

    Variant variant;

    ///////////
    //options//
    ///////////
    bool output_annotations;

    /////////
    //stats//
    /////////
    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;
    uint32_t no_vntrs_genotyped;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    LogTool lt;
    faidx_t *fai;

    /**
     * Constructor.
     */
    BCFSingleGenotypingBufferedReader(std::string input_vcf_file, std::vector<GenomeInterval>& intervals, std::string output_vcf_file);

    /**
     * Collects sufficient statistics from read for variants to be genotyped.
     */
    void process_read(bam_hdr_t *h, bam1_t *s);

    /**
     * Compute SNP genotype likelihoods in PHRED scale.
     */
    void compute_snp_pl(std::vector<int32_t>& alleles, 
                        std::vector<uint32_t>& quals, 
                        uint32_t ploidy, uint32_t no_alleles, 
                        std::vector<uint32_t>& pls,
                        float& pl_offset);

    /**
     * Compute Indel genotype likelihoods in PHRED scale.
     */
    void compute_indel_pl(std::vector<int32_t>& alleles, 
                          std::vector<uint32_t>& quals, 
                          uint32_t ploidy, uint32_t no_alleles, 
                          std::vector<uint32_t>& pls);

    /**
     * Compute Indel allele likelihoods in PHRED scale.
     */
    void compute_indel_al(char lflank_state[], uint32_t lflank_qual[], 
                          char rflank_state[], uint32_t rflank_qual[],
                          std::vector<std::string>& alleles,
                          std::string& obs_indel,
                          std::vector<uint32_t>& aqs,
                          std::vector<int32_t>& als,
                          std::string& dls);
                        
    /**
     * Collects sufficient statistics from read for variants to be genotyped.
     */
    void collect_sufficient_statistics(GenotypingRecord *g,  AugmentedBAMRecord& as);

    /**
     * Flush records.
     */
    void flush(bam_hdr_t *h, bam1_t *s, bool flush_all=false);

    /**
     * Genotype variant and print to odw.
     */
    void genotype_and_print(GenotypingRecord* g);
    
    /**
     * Create appropriate genotyping record.
     */
    GenotypingRecord* create_genotyping_record(bcf_hdr_t* h, bcf1_t* v, uint32_t ploidy, Variant& variant);
};

#endif