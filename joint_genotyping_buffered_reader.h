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

#ifndef JOINT_GENOTYPING_BUFFERED_READER_H
#define JOINT_GENOTYPING_BUFFERED_READER_H

#include <string>
#include <vector>
#include "htslib/kseq.h"
#include "htslib/vcf.h"
#include "hts_utils.h"
#include "program.h"
#include "joint_genotyping_record.h"
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
class JointGenotypingBufferedReader
{
    public:

    ///////
    //i/o//
    ///////
    BCFOrderedReader* odr; // anchor VCF

    //////////////////
    //buffer related//
    //////////////////
    std::vector<JointGenotypingRecord*> gRecords;
    std::string chrom;
    AugmentedBAMRecord as;
    Variant variant;

    int32_t lastFirst;
    //int32_t currentSampleIndex;

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

    std::vector<std::string> sample_names;
    std::vector<double> sample_contams;

    /**
     * Constructor.
     */
    JointGenotypingBufferedReader(std::string in_vcf_filename, std::vector<GenomeInterval>& intervals, std::string out_vcf_filename, int32_t nsamples);


    void set_sample(int32_t sampleIndex, const char* sampleName, double contam);
    void flush_sample(int32_t sampleIndex);
    
    /**
     * Collects sufficient statistics from read for variants to be genotyped.
     */
    int32_t process_read(bam_hdr_t *h, bam1_t *s, int32_t sampleIndex);

    inline int32_t numVariants() { return (int32_t)gRecords.size(); }
    bcf1_t* flush_variant(int32_t variantIndex, bcf_hdr_t* hdr);
    void write_header(BCFOrderedWriter* odw);
};

#endif
