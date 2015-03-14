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

#ifndef BINOMIAL_DIST_H
#define BINOMIAL_DIST_H

#include <assert.h>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <float.h>
#include <iostream>
#include <map>
#include <queue>
#include <vector>
#include "Rmath/Rmath.h"
#include "binomial_distribution.h"

/**
 * Class implementing log space arithmetic.
 */
class BinomialDistribution
{
    private:
    float p;
    std::vector<std::vector<float> > pvalues;

    public:
    /**
     * Constructor.
     */
    BinomialDistribution();

    /**
     * Constructor.
     */
    BinomialDistribution(float p);

    /**
     * Sets p.
     */
    void set_p(float p);

    /**
     * Get P(X<=x).
     */
    float get_pvalue(uint32_t x, uint32_t n);
    
    /**
     * Get pvalues size.
     */
    uint32_t get_pvalue_size();
};

#endif