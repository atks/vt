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

#ifndef MOTIF_SUFFIX_TREE_H
#define MOTIF_SUFFIX_TREE_H

#include <cstdint>
#include <cstring>
#include <vector>
#include <list>
#include <iostream>
#include <queue>
#include "candidate_motif.h"
//#include <arpa/net.h>

#define A 1
#define C 2
#define G 4
#define T 8
#define N 15

#define shift1(m) (((0x0FFFFFFF&(m))<<4) | ((0xF0000000&(m))>>28))
#define shift2(m) ((0x00FFFFFF&(m)<<8) | (0xFF000000&(m)>>24))
#define shift3(m) ((0x000FFFFF&(m)<<12) | (0xFFF00000&(m)>>20))
#define shift4(m) ((0x0000FFFF&(m)<<16) | (0xFFFF0000&(m)>>16))
#define shift5(m) ((0x00000FFF&(m)<<20) | (0xFFFFF000&(m)>>12))
#define shift6(m) ((0x000000FF&(m)<<24) | (0xFFFFFF00&(m)>>8))
#define shift7(m) ((0x0000000F&(m)<<28) | (0xFFFFFFF0&(m)>>4))

#define seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

#define index2base(i) ("OACXGXXXTXXXXXXN"[(i)])

/**
 * Motif Suffix Tree for selecting candidate motifs.
 */
class MotifSuffixTree
{
    public:
    uint64_t* tree;

    /**
     * Constructor.
     */
    MotifSuffixTree();

    /**
     * Destructor.
     */
    ~MotifSuffixTree();

    /**
     * Clear the suffix tree.
     */
    void clear();

    /**
     * Construct suffix tree based on sequence.
     */
    void set_sequence(char* sequence);

    /**
     * Construct suffix tree based on sequence up to max_motif_len.
     */
    void set_sequence(char* sequence, int32_t max_motif_len);

    /**
     * Gets candidate motifs up to max_motif_len.
     */
    void get_candidate_motifs(std::vector<CandidateMotif>& candidate_motifs);

    /**
     * Get canonical representation.
     */
    uint32_t canonical(uint32_t motif);

    private:

    /**
     * Adds a suffix of sequence from start to end.
     */
    void add_suffix(char* sequence, int32_t start, int32_t end);

    /**
     * Converts base to index.
     */
    int32_t base2index(char base);

    /**
     * Print sequence.
     */
    void print(uint32_t seq);

};

#undef A
#undef C
#undef G
#undef T
#undef N

#endif