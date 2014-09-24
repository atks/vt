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

#ifndef BAM_VARIANT_EXTRACTOR_H
#define BAM_VARIANT_EXTRACTOR_H

#include <vector>
#include <map>
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/faidx.h"
#include "program.h"
#include "hts_utils.h"
#include "variant_manip.h"
#include "utils.h"
#include "variant_buffer.h"

/**
 * Class for mining candidate variants.
 *
 * Processes an align read and records the variants observed against the reference
 */
class BAMVariantExtractor
{
    public:

    VariantBuffer *vb;
    
    std::string chrom;
    uint32_t start, end;
    uint32_t min_empty_buffer_size;
    uint32_t start_genome_pos0;
    uint32_t max_used_buffer_size_threshold;
    uint32_t max_indel_length;
    uint32_t baseq_cutoff;
    uint32_t evidence_allele_count_cutoff;
    double fractional_evidence_allele_count_cutoff;
    faidx_t *fai;
    int32_t vtype;
    kstring_t s;
    kstring_t alleles;
    kstring_t read_seq;
    kstring_t qual;
    kstring_t cigar;
    bool debug;
    bcf1_t* v;
        
    /**
     * Constructor
     * baseq_cutoff - q value cutoff to select candidate SNPs
     */
    BAMVariantExtractor(int32_t vtype,
                  size_t evidence_allele_count_cutoff,
                  double fractional_evidence_allele_count_cutoff,
                  size_t baseq_cutoff,
                  std::string& ref_fasta_file);    

    /**
     * Transfer read into a buffer for processing later
     */
    void process_read(bam_hdr_t *h, bam1_t *s);

    /**
     * Processes buffer to pick up variants
     */
    void extract_candidate_variants();

    /**
     * Empty variant buffer records that are completed.
     */
    bool flush_variant_buffer();

    /**
     * Processes buffer to pick up variant
     */
    bool next_variant(bcf1_t* v);
   
    /**
     * Checks if a variant is normalized.
     */
    bool is_biallelic_normalized(std::string& ref, std::string& alt);
    
    /**
     * Normalize a biallelic variant.
     */
    void normalize_biallelic(size_t pos0, std::string& ref, std::string& alt);

    private:
};


#endif