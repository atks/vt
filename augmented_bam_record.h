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

#ifndef AUGMENTED_CIGAR_H
#define AUGMENTED_CIGAR_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include "htslib/sam.h"
#include "hts_utils.h"

/**
 * The augmented BAM record adds functionalities to process the
 * cigar and MD5 tag in an integrated fashion.
 *
 * Augmented Cigar.
 *
 * The augmented cigar is basically a cigar that converts a normal
 * cigar to include mismatches.  An additional structure is added
 * to point to the releveant sequences.
 *
 * The goal is to simply allow a single iteration of the augmented
 * cigar for genotyping purposes.  In addition, this structure is
 * amenable to left alignment.
 *
 * Contains
 * 1. seq from bam1_t
 * 2. qual from bam1_t
 * 3. cigar from bam1_t
 * 4. MD from bam1_t (or reconstructed)
 * 5. aux_cigar that includes mismatches
 * 6. aux_seq that points to
 *    a. base substitution in MD
 *    b. inserted sequence in seq (I)
 *    c. deleted sequence in MD (D)
 *
 * For ease of left alignment of indels, and extracting a SNP
 */
class AugmentedBAMRecord
{
    public:
    bam_hdr_t *h;
    bam1_t *s;
    int32_t beg1;
    int32_t end1;
    
    //new augmented cigar with X's
    std::vector<uint32_t> aug_cigar;

    //points to mismatch, deleted and inserted sequences
    std::vector<std::string> aug_ref;

    //points to mismatch, deleted and inserted sequences
    std::vector<std::string> aug_alt;

    //statistics
    uint32_t no_mismatches;

    /**
     * Constructor.
     */
    AugmentedBAMRecord();

    /**
     * Constructor.
     */
    AugmentedBAMRecord(bam_hdr_t* h, bam1_t* s);

    /**
     * Initialize.
     */
    void initialize(bam_hdr_t* h, bam1_t* s);

    /**
     * Left align indels in an augmented cigar.
     *
     * returns
     * 1 - if left alignment was performed
     * 2 - if left alignment was not possible
     * 3 - if left alignment is possible beyond the extent of the alignment
     */
    bool left_align();

    /**
     * Right align indels in an augmented cigar.
     *
     * returns
     * 1 - if left alignment was performed
     * 2 - if left alignment was not possible
     * 3 - if left alignment is possible beyond the extent of the alignment
     */
    bool right_align();

    /**
     * Clear.
     */
    void clear();

    /**
     * Prints alignment of record.
     */
    void print();
};

#endif