/* The MIT License

   Copyright (c) Adrian Tan <atks@umich.edu>

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

#ifndef READ_FILTER_H
#define READ_FILTER_H

#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "hts_utils.h"
#include "utils.h"
#include "augmented_bam_record.h"

/**
 * For detecting overlapping reads.
 */
typedef struct
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

/**
 * Filter for reads.
 */
class ReadFilter
{
    public:

    //variables for keeping track of chromosome
    std::string chrom; //current chromosome
    int32_t tid; // current sequence id in bam

    //read filters
    uint32_t read_mapq_cutoff;
    uint16_t read_exclude_flag;
    bool ignore_overlapping_read;
    khash_t(rdict) *reads;

    /////////
    //stats//
    /////////
    uint32_t no_reads;
    uint32_t no_overlapping_reads;
    uint32_t no_passed_reads;
    uint32_t no_exclude_flag_reads;
    uint32_t no_low_mapq_reads;
    uint32_t no_unaligned_cigars;
    uint32_t no_malformed_del_cigars;
    uint32_t no_malformed_ins_cigars;
    uint32_t no_salvageable_ins_cigars;

    //////////////////
    //buffer related//
    //////////////////
    AugmentedBAMRecord as;

    ///////////
    //options//
    ///////////
    bool output_annotations;

    /////////
    //stats//
    /////////
    std::vector<std::string> sample_names;
    std::vector<double> sample_contams;
    
    /**
     * Constructor.
     */
    ReadFilter(uint32_t read_mapq_cutoff, uint16_t read_exclude_flag, bool ignore_overlapping_read);

    /**
     * Destructor.
     */
    ~ReadFilter()
    {
        kh_destroy(rdict, reads);
    }
    
    /**
     * Filter reads.
     *
     * Returns true if read is failed.
     */
    bool filter_read(bam_hdr_t* h, bam1_t *s);

    /**
     * Clear reads from hash.
     */
    void clear_reads();

    /**
     * Print BAM for debugging purposes.
     */
    void bam_print_key_values(bam_hdr_t *h, bam1_t *s);

    private:
};

#endif
