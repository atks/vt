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

#ifndef CIRCULAR_BUFFER_H
#define CIRCULAR_BUFFER_H

#include <map>
#include <vector>
#include "utils.h"
#include "hts_utils.h"
#include "variant_manip.h"

/**
 * A Reference Sequence object wrapping htslib's faidx.
 * This allows for buffered reading of seqeunces.
 */
class CircularBuffer
{
    public:
    uint32_t buffer_size;
    uint32_t buffer_size_mask;
    uint32_t window_size;
    std::vector<char> P;

    int32_t tid;
    uint32_t beg0, end0; // index of P for the start and end position in 0 base coordinates


    //  beg0 points to the beginning of the pileup
    //  end0 points to the position after the end of the pileup
    //
    //      end0
    //       |
    //       V
    //  +--+--+--+--+--+--+--+--+--+--+--+--+
    //  |  |  |  |  |  |  |  |  |  |  |  |  |
    //  +--+--+--+--+--+--+--+--+--+--+--+--+
    //       ^
    //       |
    //      beg0
    //
    //  when beg0==end0, pileup is empty
    //  gbeg1 - is the coordinate of the genome position that beg0 represents.
    //
    //  invariance:  pileup is always continuous
    //
    //  how do we check if
    //


    //where should we check if a read is from another chromosome?

    //genome position at the beginning of the pileup
    uint32_t gbeg1;

    int32_t debug;

    faidx_t *fai;

    public:

    /**
     * Get a base.
     */
    char fetch_base(std::string& chrom, uint32_t& pos1);  
   
    /**
     * Fetches sequence.
     *
     * Sequence is buffered in a circular buffer.
     */
    void fetch_seq(std::string& chrom, uint32_t start0, uint32_t end0);

    /**
     * Constructor.
     *
     * @k - size of pileup is 2^k
     */
    CircularBuffer(uint32_t k=10, uint32_t window_size=256);

    /**
     * Overloads subscript operator for accessing pileup positions.
     */
    char& operator[] (const int32_t i);

    /**
     * Returns the maximum size of the pileup.
     */
    uint32_t max_size();

    /**
     * Returns the size of the pileup.
     */
    uint32_t size();

    /**
     * Checks if buffer is empty
     */
    bool is_empty();

    /**
     * Set reference fasta file.
     */
    void set_reference(std::string& ref_fasta_file);

    /**
     * Set debug.
     */
    void set_debug(int32_t debug);

    /**
     * Sets tid.
     */
    void set_tid(uint32_t tid);

    /**
     * Gets tid.
     */
    uint32_t get_tid();

    /**
     * Sets chrom.
     */
    void set_chrom(std::string& chrom);

    /**
     * Gets chrom.
     */
    std::string get_chrom();

    /**
     * Gets window_size.
     */
    uint32_t get_window_size();

    /**
     * Sets gbeg1.
     */
    void set_gbeg1(uint32_t gbeg1);

    /**
     * Gets gbeg1.
     */
    uint32_t get_gbeg1();

    /**
     * Gets gend1.
     */
    uint32_t get_gend1();

    /**
     * Sets beg0.
     */
    void set_beg0(uint32_t beg0);

    /**
     * Sets end0.
     */
    void set_end0(uint32_t end0);

    /**
     * Gets the index of the first element.
     */
    uint32_t begin();

    /**
     * Gets the index of the last element.
     */
    uint32_t end();

    /**
     * Increments i by 1 circularly.
     */
    uint32_t inc(uint32_t i);

    /**
     * Increments beg0 by 1.
     */
    void inc_beg0();

    /**
     * Increments end0 by 1.
     */
    void inc_end0();

    /**
     * Increments index i by j cyclically.
     */
    uint32_t inc(uint32_t i, uint32_t j);

    /**
     * Converts gpos1 to index in P.
     */
    uint32_t g2i(uint32_t gpos1);

    /**
     * Returns the difference between 2 buffer positions
     */
    uint32_t diff(uint32_t i, uint32_t j);

    /**
     * Print pileup state.
     */
    void print_state();

    /**
     * Print pileup state extent.
     */
    void print_state_extent();
};

#endif