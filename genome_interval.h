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

#ifndef GENOME_INTERVAL_H
#define GENOME_INTERVAL_H

#include <queue>
#include <sstream>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include "hts_utils.h"
#include "utils.h"

class GenomeInterval
{
    public:
    std::string seq;
    int32_t start1; // 1 based coordinate
    int32_t end1;   // 1 based coordinate

    /**
     * Constructs a Genome Interval.
     */
    GenomeInterval() {};

    /**
     * Constructs a Genome Interval.
     */
    GenomeInterval(std::string& seq, int32_t start1, int32_t end1);

    /**
     * Constructs a Genome Interval from a string representation.
     *
     * @interval    string representation of an interval.
     *
     * e.g X:2000-4000   position 2000 to 4000 on chromosome X
     *     Y             the entirety of chromosome Y
     */
    GenomeInterval(std::string interval);

    /**
     * Sets an interval.
     */
    void set(std::string& seq, int32_t start1, int32_t end1);

    /**
     * Sets an interval.
     */
    void set(std::string interval);

    /**
     * Converts genome interval into the entire chromosome.
     */
    void chromosomify();

    /**
     * Returns a string representation of this Genome Interval.
     */
    std::string to_string();

    /**
     * Returns a string representation of this Genome Interval.
     */
    void to_string(kstring_t *interval);

    /**
     * Checks if this interval overlap.
     */
    bool overlaps_with(std::string& chrom, int32_t start1, int32_t end1);
};

#endif