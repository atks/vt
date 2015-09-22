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

#include "reference_sequence.h"

/**
 * Constructor.
 *
 * @k - size of buffered sequence is 2^k
 */
ReferenceSequence::ReferenceSequence(uint32_t k, uint32_t window_size)
{
    //Buffer size is a power of 2^k.
    buffer_size = 1 << k;
    //this provides a cheaper way to do modulo operations for a circular array.
    buffer_size_mask = (0xFFFFFFFF >> (32-k));
    this->window_size = window_size;

    seq.resize(buffer_size);

    beg0 = end0 = 0;
    gbeg1 = 0;

    debug = 0;
};

/**
 * Get a base.
 */
char ReferenceSequence::fetch_base(std::string& chrom, uint32_t& pos1)
{
    
    //check buffer and retrieve base if it is in it.
    if (this->chrom == chrom && pos1>=gbeg1 && pos1<=gbeg1+(end0-beg0))
    {
        
    }
    else
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

    }
    
    return 'N';
}

/**
 * Fetches sequence chrom:start1-end1
 * 
 * Retrieved sequence is in seq with the length of n.
 */
void ReferenceSequence::fetch_seq(std::string& chrom, uint32_t start1, uint32_t end1, char* seq, int32_t n)
{
}

/**
 * Overloads subscript operator for accessing buffered sequence positions.
 */
char& ReferenceSequence::operator[] (const int32_t i)
{
    return seq[i];
}

/**
 * Returns the maximum size of the buffered sequence.
 */
uint32_t ReferenceSequence::max_size()
{
    return buffer_size - 1;
}

/**
 * Returns the size of the buffered sequence.
 */
uint32_t ReferenceSequence::size()
{
    return (end0>=beg0 ? end0-beg0 : buffer_size-(beg0-end0));
}

/**
 * Checks if buffer is empty.
 */
bool ReferenceSequence::is_empty()
{
    return beg0==end0;
};

/**
 * Set reference fasta file.
 */
void ReferenceSequence::set_reference(std::string& ref_fasta_file)
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
void ReferenceSequence::set_debug(int32_t debug)
{
    this->debug = debug;
};

/**
 * Sets tid.
 */
void ReferenceSequence::set_tid(uint32_t tid)
{
    this->tid = tid;
}

/**
 * Gets tid.
 */
uint32_t ReferenceSequence::get_tid()
{
    return this->tid;
}

/**
 * Sets chrom.
 */
void ReferenceSequence::set_chrom(std::string& chrom)
{
    this->chrom = chrom;
}

/**
 * Gets chrom.
 */
std::string ReferenceSequence::get_chrom()
{
    return chrom;
}

/**
 * Gets window_size.
 */
uint32_t ReferenceSequence::get_window_size()
{
    return window_size;
}

/**
 * Converts gpos1 to index in seq.
 * If P is empty, initialize first position as gpos1.
 */
uint32_t ReferenceSequence::g2i(uint32_t gpos1)
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

