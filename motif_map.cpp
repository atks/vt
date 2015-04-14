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

#include "motif_map.h"

#define A 1
#define C 2
#define G 4
#define T 8
#define N 15

/**
 * Constructor.
 */
MotifMap::MotifMap()
{
    uint32_t size = (1<<2) + (1<<4) + (1<<6) + (1<<8) + (1<<10) + (1<<12) + (1<<14) + (1<<16);

    //len_count = (uint32_t *) malloc(sizeof(uint32_t)*9);

    len_count.resize(8,0);
    for (uint32_t i=1; i<8; ++i)
    {
        len_count[i] = len_count[i-1] + (1<<(2*i));
        std::cerr << "\t1<< " << (2*i) << " : " << (1<<(2*i)) << " " << len_count[i] << "\n";
    }

    std::cerr << "size : " << size << "\n";
    std::cerr << "\t1<<2 : " << (1<<2) << "\n";
    std::cerr << "\t1<<4 : " << (1<<4) << "\n";
    std::cerr << "\t1<<6 : " << (1<<6) << "\n";
    std::cerr << "\t1<<8 : " << (1<<8) << "\n";
    std::cerr << "\t1<<10 : " << (1<<10) << "\n";
    std::cerr << "\t1<<12 : " << (1<<12) << "\n";
    std::cerr << "\t1<<14 : " << (1<<14) << "\n";
    std::cerr << "\t1<<16 : " << (1<<16) << "\n";
    std::cerr << "\tmax int16_t : " << 0xFFFF << "\n";

   
    uint64_t value = 0;
    //you only want to loop through multiples of 2
//    for (uint32_t i=0; i<size; ++i)
//    {
//        uint32_t c = canonical(i);
//        tree[i] = ((uint64_t) c)<<32;
//        std::cerr << i << ":" << c << ":" << tree[i] << "\n";
//        if (i==1) exit(1);
//    }

    //enumerate size

    //map index to sequence


    //map sequence to index - 4 based

    uint32_t next_len_index;
    uint32_t clen = 1;
    uint32_t seq = 0;

exit(1);

    for (uint32_t i=0; i<size; ++i)
    {
        //node tnode = tree[i];
        if (i==next_len_index)
        {
            ++clen;
        }

        //for the len - extract base
        uint32_t c = canonical(seq);


    }

};

/**
 * Destructor.
 */
MotifMap::~MotifMap()
{
};

/**
 * Construct suffix tree based on sequence.
 */
void MotifMap::set_sequence(char* sequence)
{
    //translate sequence to binary form
    uint32_t len = strlen(sequence);



};

/**
 * Construct suffix tree based on sequence up to max_motif_len.
 */
void MotifMap::set_sequence(char* sequence, int32_t max_motif_len)
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
uint32_t MotifMap::index2sequence(uint32_t index)
{
    return 0;
}

/**
 * Gets index of child.
 */
uint32_t MotifMap::get_first_child(uint32_t index)
{
    return 0;
}

/**
 * Adds a suffix of sequence from start to end.
 */
void MotifMap::add_suffix(char* sequence, int32_t start, int32_t end)
{

};

/**
 * Converts base to index.
 */
int32_t MotifMap::base2index(char base)
{
    switch (base)
    {
        case 'A':
            return A;
            break;
        case 'C':
            return C;
            break;
        case 'G':
            return G;
            break;
        case 'T':
            return T;
            break;
        default:
            return N;
    }
};

#undef A
#undef C
#undef G
#undef T
#undef N
