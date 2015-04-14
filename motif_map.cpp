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

#include "motif_map.h"

#define A 1
#define C 2
#define G 4
#define T 8
#define N 15

/**
 * Constructor.
 */
MotifMap::MotifMap(uint32_t size)
{
    if (size>8)
    {
        fprintf(stderr, "[%s:%d %s] size > 8  not supported : %d\n", __FILE__, __LINE__, __FUNCTION__, size);
        exit(1);
    }    
    
    uint32_t max_size = 0;
    len_count.resize(size,0);
    for (uint32_t i=1; i<size; ++i)
    {
        len_count[i] = len_count[i-1] + (1<<(2*i));
        std::cerr << "\t1<< " << (2*i) << " : " << (1<<(2*i)) << " " << len_count[i] << "\n";
        max_size += (1<<(2*i));
    }
};

/**
 * Destructor.
 */
MotifMap::~MotifMap()
{
};

/**
 * Get canonical representation.
 */
uint32_t MotifMap::canonical(uint32_t motif)
{
    uint32_t cmotif = motif;
    uint32_t smotif = motif;
    std::cerr << "\t" << 0 << ") " << smotif << " - ";
    std::cerr << "\n";
    for (uint32_t i=1; i<8; ++i)
    {
        std::cerr << "\t" << i << ") " << smotif << " - ";
        std::cerr << "\n";
        cmotif = smotif<cmotif ? smotif : cmotif;
    }

    return cmotif;
}

/**
 * Converts index to sequence.
 */
uint32_t MotifMap::index2seq(uint32_t index)
{
    return 0;
}

/**
 * Converts index to sequence.
 */
uint32_t MotifMap::seq2index(uint32_t index)
{
    return 0;
}

#undef A
#undef C
#undef G
#undef T
#undef N
