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

#include "hts_utils.h"
#include <map>
#include "utils.h"
#include "variant_manip.h"
#include <vector>

/**
 * A circular buffer for containing SNPs and Indels.
 */
class VariantBuffer
{
    public:

    /**
     * Constructor
     */
    VariantBuffer(size_t buffer_size=1000);

    size_t buffer_size;
    std::vector<char> R;                      // reference bases
    std::vector<std::vector<char> > X;        // contains read bases that differ from the reference
    std::vector<std::vector<std::string> > I; // contains inserted bases at the left anchor position including anchor
    std::vector<std::vector<std::string> > J; // contains inserted bases at the right anchor position including anchor
    std::vector<std::vector<std::string> > D; // contains deleted bases at the anchor position including anchor
    std::vector<int32_t> N;                   // number of alternative evidences observed here - sum of R, X, D and I

    size_t start, end;                        // non empty location in the buffer
    size_t empty_buffer_space;                // remaining buffer space left
    size_t start_genome_pos0;

    //not necessary?
    size_t min_empty_buffer_size;
    size_t max_used_buffer_size_threshold;
    size_t max_indel_length;

    bool debug;

    /**
     * Inserts a reference base at pos0 into the buffer.
     */
    void insertR(size_t pos0, char r);

    /**
     * Inserts an alternate base at pos0 into the buffer.
     */
    void insertX(size_t pos0, char x);

    /**
     * Inserts a deletion at pos0 into the buffer.
     */
    void insertD(size_t pos0, std::string& ref, std::string& alt);

    /**
     * Inserts an insertion base at pos0 into the buffer.
     */
    void insertI(size_t pos0, std::string& ref, std::string& alt);

    private:

    /**
     * Checks if buffer is empty
     */
    bool is_empty();

    /**
     *Increments buffer index i by 1.
     */
    void add(size_t& i);

    /**
     * Increments buffer index i by j.
     */
    size_t add(size_t i, size_t j);

    /**
     * Decrements buffer index i by j.
     */
    size_t minus(size_t& i, size_t j);

    /**
     * Decrements buffer index i by 1.
     */
    void minus(size_t& i);

    /**
     * Returns the difference between 2 buffer positions
     */
    size_t diff(size_t i, size_t j);

    /**
     * Gets the position in the buffer that corresponds to
     * the genome position indicated by pos.
     */
    size_t get_cur_pos0(size_t genome_pos0);

    /**
     * Print buffer contents for debugging purpose
     */
    void printBuffer();
};

#endif