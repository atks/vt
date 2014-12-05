/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#include "genome_interval.h"

/**
 * Constructs a Genome Interval.
 */
GenomeInterval::GenomeInterval(std::string& seq, int32_t start1, int32_t end1)
{
    set(seq, start1, end1);
};

/**
 * Constructs a Genome Interval from a string representation.
 *
 * @interval    string representation of an interval.
 *
 * e.g X:2000-4000   position 2000 to 4000 on chromosome X
 *     Y:2000        position 2000 on chromosome Y
 *     Y             the entirety of chromosome Y
 */
GenomeInterval::GenomeInterval(std::string interval)
{
    set(interval);
};

/**
 * Sets an interval.
 */
void GenomeInterval::set(std::string& seq, int32_t start1, int32_t end1)
{
    this->seq = seq;
    this->start1 = start1;
    this->end1 = end1;
};

/**
 * Sets an interval.
 */
void GenomeInterval::set(std::string interval)
{
    std::vector<std::string> v;
    split(v, ":-", interval);

    if (v.size()==1)
    {
        seq = v[0];
        start1 = 1;
        end1 = (1<<29) - 1;
    }
    else if (v.size()==2)
    {
        seq = v[0];
        if (!str2int32(v[1], start1))
        {
            fprintf(stderr, "[%s:%d %s] Invalid genomic interval: %s\n", __FILE__,__LINE__,__FUNCTION__, interval.c_str());
            exit(1);
        }

        end1 = start1;
    }
    else if (v.size()==3)
    {
        seq = v[0];
        if (!str2int32(v[1], start1) || !str2int32(v[2], end1))
        {
            fprintf(stderr, "[%s:%d %s] Invalid genomic interval: %s\n", __FILE__,__LINE__,__FUNCTION__, interval.c_str());
            exit(1);
        }
    }
    else
    {
        fprintf(stderr, "[%s:%d %s] Invalid genomic interval: %s\n", __FILE__,__LINE__,__FUNCTION__, interval.c_str());
        exit(1);
    }
}

/**
 * Converts genome interval into the entire chromosome.
 */
void GenomeInterval::chromosomify()
{
    start1 = 1;
    end1 = (1<<29) - 1;
}

/**
 * Returns a string representation of this Genome Interval.
 */
std::string GenomeInterval::to_string()
{
    kstring_t s = {0,0,0};
    kputs(seq.c_str(), &s);
    if (start1!=1 || end1!=((1<<29)-1))
    {
        kputc(':', &s);
        kputw(start1, &s);
        kputc('-', &s);
        kputw(end1, &s);
    }
    std::string interval(s.s);
    if (s.m) free(s.s);
    return interval;
};

/**
 * Returns a string representation of this Genome Interval.
 */
void GenomeInterval::to_string(kstring_t *interval)
{
    interval->l = 0;
    kputs(seq.c_str(), interval);
    if (start1!=1 || end1!=((1<<29)-1))
    {
        kputc(':', interval);
        kputw(start1, interval);
        kputc('-', interval);
        kputw(end1, interval);
    }
};

/**
 * Checks if this interval overlap.
 */
bool GenomeInterval::overlaps_with(std::string& chrom, int32_t start1, int32_t end1)
{
    return (seq==chrom && this->start1<=end1 && this->end1>=start1);
};