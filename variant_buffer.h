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

#ifndef VARIANT_BUFFER_H
#define VARIANT_BUFFER_H

#include <vector>
#include <map>
#include "hts_utils.h"
#include "variant_manip.h"
#include "utils.h"

/**
 * Class for mining candidate variants.
 *
 * Processes an align read and records the variants observed against the reference
 */
class VariantBuffer
{
    public:
    /**
     * Constructor
     * baseq_cutoff - q value cutoff to select candidate SNPs
     */
    VariantBuffer(uint32_t vtype,
                  uint32_t evidence_allele_count_cutoff,
                  double fractional_evidence_allele_count_cutoff,
                  uint32_t baseq_cutoff,
                  faidx_t *fai);

    private:

    uint32_t buffer_size;
    std::vector<std::vector<char> > X; // contains read bases that differ from the genome
    std::vector<std::vector<std::string> > I; //contains inserted bases
    std::vector<std::vector<std::string> > D; //contains reference bases that are deleted
    std::vector<int32_t> N; // number of evidences observed here - combination of X, I and D
    std::vector<char> REF; 
    std::vector<char> ANCHOR;
    std::vector<std::string> ALT;
    char* chrom;

    //key control variables for circular buffer
    uint32_t start, end;
    uint32_t empty_buffer_space;
    uint32_t min_empty_buffer_size;
    uint32_t start_genome_pos0;
    uint32_t max_used_buffer_size_threshold;
    uint32_t max_indel_length;
    uint32_t baseq_cutoff;
    uint32_t evidence_allele_count_cutoff;
    double fractional_evidence_allele_count_cutoff;
    faidx_t *fai;
    uint32_t vtype;
    kstring_t s;
    kstring_t alleles;
    kstring_t read_seq;
    kstring_t qual;
    kstring_t cigar;

    bool debug;

    /**
     * Processes buffer to pick up variants
     * Empty buffer to recover space.
     * @chrom  - remove variants on chrom
     * @pos1   - remove variants up to pos1
     * @flush  - remove all variants
     */
    void extract_candidate_variants(const char* chrom, uint32_t pos1, bool flush=false);

    /**
     * Checks if buffer is empty
     */
    bool is_empty();

    /**
     *Increments buffer index i by 1.
     */
    void add(uint32_t& i);

    /**
     * Increments buffer index i by j.
     */
    uint32_t add(uint32_t i, uint32_t j);

    /**
     * Decrements buffer index i by j.
     */
    uint32_t minus(uint32_t& i, uint32_t j);

    /**
     * Decrements buffer index i by 1.
     */
    void minus(uint32_t& i);

    /**
     * Returns the difference between 2 buffer positions
     */
    uint32_t diff(uint32_t i, uint32_t j);

    /**
     * Gets the position in the buffer that corresponds to
     * the genome position indicated by pos.
     */
    uint32_t get_cur_pos0(uint32_t genome_pos0);

    /**
     * Print buffer contents for debugging purpose
     */
    void printBuffer();
};


#endif