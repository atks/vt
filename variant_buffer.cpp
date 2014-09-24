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

#include "variant_buffer.h"

VariantBuffer::VariantBuffer(size_t buffer_size)
{
    R.resize(buffer_size);
    X.resize(buffer_size);
    I.resize(buffer_size);
    D.resize(buffer_size);
    N.resize(buffer_size, 0);
    
    start = 0;
    end = 0;
    
    empty_buffer_space = buffer_size;

    start_genome_pos0 = 0;    
};


/**
 * Inserts a reference base at pos0 into the buffer.
 */
void VariantBuffer::insertR(size_t pos0, char r)
{
    if (pos0)
    {
    }
};

/**
 * Inserts an alternate base at pos0 into the buffer.
 */
void VariantBuffer::insertX(size_t pos0, char x)
{
};

/**
 * Inserts a deletion at pos0 into the buffer.
 */       
void VariantBuffer::insertD(size_t pos0, std::string& ref, std::string& alt)
{
};

/**
 * Inserts an insertion base at pos0 into the buffer.
 */
void VariantBuffer::insertI(size_t pos0, std::string& ref, std::string& alt)
{
};

/**
 * Checks if buffer is empty.
 */
bool VariantBuffer::is_empty()
{
    return start==end;
};

/**
 *Increments buffer index i by 1.
 */
void VariantBuffer::add(size_t& i)
{
    if (i>=buffer_size)
    {
        std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
        exit(1);
    }

    size_t temp = (i+1)%buffer_size;
    i = end==i ? (end=temp) : temp;
};

/**
 * Increments buffer index i by j.
 */
size_t VariantBuffer::add(size_t i, size_t j)
{
    if (i>=buffer_size)
    {
        std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
        exit(1);
    }

    return (i+j)%buffer_size;
};

/**
 * Decrements buffer index i by j.
 */
size_t VariantBuffer::minus(size_t& i, size_t j)
{
    if (i>=buffer_size)
    {
        std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
        exit(1);
    }

    return (i>=j ? i-j : buffer_size-(j-i));
};

/**
 * Decrements buffer index i by 1.
 */
void VariantBuffer::minus(size_t& i)
{
    if (i>=buffer_size)
    {
        std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
        exit(1);
    }

    i = (i>=1 ? i-1 : buffer_size-1);
};

/**
 * Returns the difference between 2 buffer positions
 */
size_t VariantBuffer::diff(size_t i, size_t j)
{
    return (i>=j ? i-j : buffer_size-(j-i));
};

/**
 * Gets the position in the buffer that corresponds to
 * the genome position indicated by pos.
 */
size_t VariantBuffer::get_cur_pos0(size_t genome_pos0)
{
    //when buffer is empty
    if (is_empty())
    {
        start_genome_pos0 = genome_pos0;
        return start;
    }
    else
    {
        if (genome_pos0-start_genome_pos0>buffer_size)
        {
            std::cerr << "overflow buffer\n" ;
            //should allow for unbuffering here

        }
        return (start + (genome_pos0-start_genome_pos0))%buffer_size;
    }
};

/**
 * Print buffer contents for debugging purpose
 */
void VariantBuffer::printBuffer()
{
    std::cout << "PRINT BUFFER" << "\n";
    std::cout << "usedBufferSize: " << diff(end,start) << "\n";
    size_t cur_pos0 = start;
    size_t genome_pos = start_genome_pos0;

    while (cur_pos0!=end)
    {
        std::cerr << genome_pos << "\t" << cur_pos0 << "\t" << R[cur_pos0] << "\t";

        for (size_t j=0; j<I[cur_pos0].size(); ++j)
        {
            std::cerr << I[cur_pos0][j] << ",";
        }

        for (size_t j=0; j<D[cur_pos0].size(); ++j)
        {
            std::cerr << D[cur_pos0][j] << ",";
        }  

        std::cerr << "\t" <<  N[cur_pos0] << "\n";

        add(cur_pos0);
        ++genome_pos;
    }
};
