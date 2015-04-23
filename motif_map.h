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

#define get_seqi(s, i) (((s)>>((15-(i))<<1)) & 3)
#define set_seqi(s, i, b) (((b)<<((15-(i))<<1)) | (s&~(((3)<<((15-(i))<<1)))))
//shift a sequence s of length l
#define shift(s,l)  ((((s) >> 30) << ((16-(l))<<1)) | ((s) << 2))


/**
 * Motif Map for the mapping functions to and from index.
 *
 * Function for canonical form too.
 */
class MotifMap
{
    public:
    std::vector<uint32_t> len_count;
    uint32_t max_len;  
    uint32_t max_index;    
    
    /**
     * Constructor.
     */
    MotifMap(uint32_t max_len);

    /**
     * Destructor.
     */
    ~MotifMap();

    /**
     * Get canonical representation.
     */
    uint32_t canonical(uint32_t motif, uint32_t len);

    /**
     * Checks if a string is aperiodic.
     */
    bool is_aperiodic(uint32_t motif, uint32_t len);

    /**
     * Converts index to sequence.
     */
    uint32_t index2seq(uint32_t index);

    /**
     * Converts sequence to index.
     */
    uint32_t seq2index(uint32_t seq, uint32_t len);
    
    /**
     * Prints sequence.
     */
    void print_seq(uint32_t seq, uint32_t len);
    
    /**
     * Converts sequence to string.
     */
    std::string seq2str(uint32_t seq, uint32_t len);

    private:
};

#endif