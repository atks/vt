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

#define A 0
#define C 1
#define G 2
#define T 3

/**
 * Constructor.
 */
MotifMap::MotifMap(uint32_t max_len)
{
    if (max_len>16)
    {
        fprintf(stderr, "[%s:%d %s] max_len > 16  not supported : %d\n", __FILE__, __LINE__, __FUNCTION__, max_len);
        exit(1);
    }

    this->max_len = max_len;

    uint32_t max_index = 0;
    len_count.resize(max_len+1,0);
    std::cerr << "\tmax_len: " << max_len << "\n";

    for (uint32_t i=1; i<=max_len; ++i)
    {
        len_count[i] = len_count[i-1] + (1<<(2*i));
        std::cerr << "\t" << i << " : " << (1<<(2*i)) << " " << len_count[i] << "\n";
        max_index += (1<<(2*i));
    }

    std::cerr << "\tmax index: " << (max_index-1) << "\n\n";

    //perform mapping.
    for (uint32_t len=1; len<=max_len; ++len)
    {
        std::cerr << "\tlen: " << len << "\n";
        uint32_t max_index = len_count[len-1] + (1<<(2*len));

        for (uint32_t index=len_count[len-1]; index<max_index; ++index)
        {
            std::cerr << "\t" << index << "\n";
            index2seq(index);
        }

        if (len==3) exit(1);    
    }
}

/**
 * Destructor.
 */
MotifMap::~MotifMap()
{
}

/**
 * Get canonical representation.
 */
uint32_t MotifMap::canonical(uint32_t motif)
{
    uint32_t cmotif = motif;
    uint32_t smotif = motif;
    std::cerr << "\t" << 0 << ") " << smotif << " - ";
    std::cerr << "\n";
    for (uint32_t i=1; i<max_len; ++i)
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
    uint32_t s = 0;
    for (uint32_t len=1; len<max_len; ++len)
    {
        if (index<len_count[len])
        {
            //std::cerr << "len: " << (len) << "\n";
            index -= len_count[len-1];
            for (int32_t i=len; i>0; --i)
            {
                uint32_t j = ((index >> ((i-1)<<1)) & 3);
//                std::cerr << "\tindex: " << index << "\n";
//                std::cerr << "\ti: " << i << "\n";
//                std::cerr << "\ts: " << s << "\n";
//                std::cerr << "\tj: " << (i-1) << "\n";
//                std::cerr << "\tbase: " << "ACGT"[j] << "\n";
                
                s = set_seqi(s,(len-i),j);
            }
            std::cerr << "seq : ";
            print_seq(s, len);

            return 0;
        }
    }

    return 0;
}

/**
 * Converts sequence to index.
 */
uint32_t MotifMap::seq2index(uint32_t seq, uint32_t len)
{
//    for (uint32_t len=1; len<max_len; ++len)
//    {
//        if (index<len_count[len])
//        {
//            //std::cerr << "len: " << (len) << "\n";
//            index -= len_count[len-1];
//            for (int32_t i=len; i>0; --i)
//            {
//                uint32_t j = ((index >> ((i-1)<<1)) & 3);
////                std::cerr << "s: " << j << " ";
////                std::cerr << "j: " << j << " ";
//                std::cerr << "ACGT"[j] << ""; 
//            }
//            std::cerr << "\n";
//
//            return 0;
//        }
//    }

    return 0;
}

/**
 * Prints sequence.
 */
void MotifMap::print_seq(uint32_t seq, uint32_t len)
{
    for (int32_t i=0; i<len; ++i)
    {
        std::cerr << "ACGT"[get_seqi(seq, i)] << "";
    }
    std::cerr << "\n";
};


#undef A
#undef C
#undef G
#undef T
