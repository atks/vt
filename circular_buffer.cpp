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

#include "circular_buffer.h"

/**
 * Constructor.
 *
 * @k - size of pileup is 2^k
 */
CircularBuffer::CircularBuffer(uint32_t k, uint32_t window_size)
{
    //Buffer size is a power of 2^k.
    buffer_size = 1 << k;
    //this provides a cheaper way to do modulo operations for a circular array.
    buffer_size_mask = (0xFFFFFFFF >> (32-k));
    this->window_size = window_size;

    P.resize(buffer_size);

    tid = -1;
    beg0 = end0 = 0;
    gbeg1 = 0;

    debug = 0;
};

/**
 * Overloads subscript operator for accessing pileup positions.
 */
char& CircularBuffer::operator[] (const int32_t i)
{
    return P[i];
}

/**
 * Returns the maximum size of the pileup.
 */
uint32_t CircularBuffer::max_size()
{
    return buffer_size - 1;
}

/**
 * Returns the size of the pileup.
 */
uint32_t CircularBuffer::size()
{
    return (end0>=beg0 ? end0-beg0 : buffer_size-(beg0-end0));
}

/**
 * Checks if buffer is empty.
 */
bool CircularBuffer::is_empty()
{
    return beg0==end0;
};

/**
 * Set reference fasta file.
 */
void CircularBuffer::set_reference(std::string& ref_fasta_file)
{
    if (ref_fasta_file!="")
    {
        fai = fai_load(ref_fasta_file.c_str());
        if (fai==NULL)
        {
            fprintf(stderr, "[%s:%d %s] cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
            exit(1);
        }
    }
};

/**
 * Set debug.
 */
void CircularBuffer::set_debug(int32_t debug)
{
    this->debug = debug;
};

/**
 * Sets tid.
 */
void CircularBuffer::set_tid(uint32_t tid)
{
    this->tid = tid;
}

/**
 * Gets tid.
 */
uint32_t CircularBuffer::get_tid()
{
    return this->tid;
}

/**
 * Sets chrom.
 */
void CircularBuffer::set_chrom(std::string& chrom)
{
    this->chrom = chrom;
}

/**
 * Gets chrom.
 */
std::string CircularBuffer::get_chrom()
{
    return chrom;
}

/**
 * Gets window_size.
 */
uint32_t CircularBuffer::get_window_size()
{
    return window_size;
}

/**
 * Converts gpos1 to index in P.
 * If P is empty, initialize first position as gpos1.
 */
uint32_t CircularBuffer::g2i(uint32_t gpos1)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        return beg0;
    }
    else
    {
        if (gpos1<gbeg1)
        {
            fprintf(stderr, "[%s:%d %s] buffer out of extent: gpos1 %d < gbeg1 %d\n", __FILE__, __LINE__, __FUNCTION__, gpos1, gbeg1);
            abort();
        }
        return (beg0 + (gpos1-gbeg1)) & buffer_size_mask;
    }
}

/**
 * Increments i by 1 circularly.
 */
uint32_t CircularBuffer::inc(uint32_t i)
{
    return (i+1) & buffer_size_mask;
};

/**
 * Sets gbeg1.
 */
void CircularBuffer::set_gbeg1(uint32_t gbeg1)
{
    this->gbeg1 = gbeg1;
}

/**
 * Gets gbeg1.
 */
uint32_t CircularBuffer::get_gbeg1()
{
    return gbeg1;
}

/**
 * Gets gend1.
 */
uint32_t CircularBuffer::get_gend1()
{
    if (is_empty())
    {
        return 0;
    }
    else
    {
        return gbeg1 + diff(end0, beg0) - 1;
    }
}

/**
 * Sets beg0.
 */
void CircularBuffer::set_beg0(uint32_t beg0)
{
    this->beg0 = beg0;
}

/**
 * Sets end0.
 */
void CircularBuffer::set_end0(uint32_t end0)
{
    this->end0 = end0;
}

/**
 * Gets the index of the first element.
 */
uint32_t CircularBuffer::begin()
{
    return beg0;
}

/**
 * Gets the index of the last element.
 */
uint32_t CircularBuffer::end()
{
    return end0;
}

/**
 * Returns the difference between 2 buffer positions
 */
uint32_t CircularBuffer::diff(uint32_t i, uint32_t j)
{
    return (i>=j ? i-j : buffer_size-(j-i));
};

/**
 * Increments beg0 by 1.
 */
void CircularBuffer::inc_beg0()
{
    beg0 = (beg0+1) & buffer_size_mask;
};

/**
 * Increments end0 by 1.
 */
void CircularBuffer::inc_end0()
{
    end0 = (end0+1) & buffer_size_mask;
};

/**
 * Increments index i by j cyclically.
 */
uint32_t CircularBuffer::inc(uint32_t i, uint32_t j)
{
    return (i+j) & buffer_size_mask;
};

/**
 * Get a base.
 */
char CircularBuffer::fetch_base(std::string& chrom, uint32_t& pos1)
{
    int ref_len = 0;
    char *refseq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos1-1, pos1-1, &ref_len);
    if (!refseq)
    {
        fprintf(stderr, "[%s:%d %s] failure to extrac base from fasta file: %s:%d: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos1-1);
        exit(1);
    }
    char base = refseq[0];
    free(refseq);

    return base;
};

///**
// * Get a sequence.  User have to free the char* returned.
// */
//char* CircularBuffer::get_sequence(std::string& chrom, uint32_t pos1, uint32_t len)
//{
//    int ref_len = 0;
//    char* seq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos1-1, pos1+len-2, &ref_len);
//    if (!seq || ref_len!=len)
//    {
//        fprintf(stderr, "[%s:%d %s] failure to extract sequence from fasta file: %s:%d: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos1-1);
//        exit(1);
//    }
//
//    return seq;
//};

/**
 * Print pileup state.
 */
void CircularBuffer::print_state()
{
    std::cerr << "******************" << "\n";
    std::cerr << "gindex   : " << gbeg1 << "-" << get_gend1() << "\n";
    std::cerr << "index   : " << beg0 << "-" << end0 << " (" << size() << ")\n";
    std::cerr << "******************" << "\n";
    uint32_t k = 0;
    for (uint32_t i=beg0; i!=end0; i=inc(i))
    {
        //P[i].print(gbeg1+k);
        ++k;
    }
    std::cerr << "******************" << "\n";
}

/**
 * Print pileup state extent.
 */
void CircularBuffer::print_state_extent()
{
    std::cerr << "******************" << "\n";
    std::cerr << "gindex   : " << gbeg1 << "-" << get_gend1() << "\n";
    std::cerr << "index   : " << beg0 << "-" << end0 << " (" << size() << ")\n";
    std::cerr << "******************" << "\n";
}