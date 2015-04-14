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

#ifndef MOTIF_MAP_H
#define MOTIF_MAP_H

#include <cstdint>
#include <vector>
#include <iostream>
#include "motif_map.h"

/**
 * Motif Map for the mapping functions to and from index.
 *
 * Function for canonical form too.
 */
class MotifMap
{
    public:
    std::vector<uint32_t> len_count;

    /**
     * Constructor.
     */
    MotifMap(uint32_t size);

    /**
     * Destructor.
     */
    ~MotifMap();

    /**
     * Get canonical representation.
     */
    uint32_t canonical(uint32_t motif);

    /**
     * Converts index to sequence.
     */
    uint32_t index2seq(uint32_t index);

    /**
     * Converts sequence to index.
     */
    uint32_t seq2index(uint32_t index);

    private:
};

#undef A
#undef C
#undef G
#undef T
#undef N

#endif