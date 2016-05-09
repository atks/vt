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

#ifndef BCF_GENOTYPING_BUFFERED_READER_H
#define BCF_GENOTYPING_BUFFERED_READER_H

#include <string>
#include "htslib/kseq.h"
#include "htslib/vcf.h"
#include "hts_utils.h"
#include <list>
#include "program.h"
#include "genotyping_record.h"
#include "bcf_ordered_reader.h"
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
class BCFGenotypingBufferedReader
{
    public:

    ///////
    //i/o//
    ///////
    BCFOrderedReader* odr;

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

    /**
     * Constructor.
     */
    BCFGenotypingBufferedReader(std::string filename, std::vector<GenomeInterval>& intervals);

    /**
     * Collects sufficient statistics from read for variants to be genotyped.
     */
    void process_read(bam_hdr_t *h, bam1_t *s);

    /**
     * Compute SNP genotype likelihoods in PHRED scale.
     */
    void compute_snp_pl(std::vector<int32_t>& alleles, std::vector<uint32_t>& quals, uint32_t ploidy,  std::vector<uint32_t>& pls);

    /**
     * Compute Indel genotype likelihoods in PHRED scale.
     */
    void compute_indel_pl(std::vector<int32_t>& alleles, std::vector<uint32_t>& quals, uint32_t ploidy,  std::vector<uint32_t>& pls);

    /**
     * Collects sufficient statistics from read for variants to be genotyped.
     */
    void collect_sufficient_statistics(GenotypingRecord *g,  AugmentedBAMRecord& as);

    /**
     * Flush records.
     */
    void flush(BCFOrderedWriter* odw, bam_hdr_t *h, bam1_t *s, bool flush_all=false);

    /**
     * Genotype variant and print to odw.
     */
    void genotype_and_print(BCFOrderedWriter* odw, GenotypingRecord* g);

};

#endif
