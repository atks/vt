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

#ifndef LOG_TOOL_H
#define LOG_TOOL_H

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

#define LOGZERO -DBL_MAX

/**
 * Class implementing log space arithmetic.
 */
class LogTool
{
    private:
    std::vector<double> PL;
    std::vector<double> PL_one_minus_p;
    std::vector<double> LOG10_VARP;
    std::vector<double> LOG10FACT;

    public:
    LogTool () {};

    /**
     * Convert -10log(p) to p.
     */
    double pl2prob(uint32_t PL);
    
    /**
     * Convert -10log(p) to -10log(1-p).
     */
    double pl2pl_one_minus_p(uint32_t pl);

    /**
     * Convert -10log(p) to log10(sqrt(e(1-e))).
     */
    double pl2log10_varp(uint32_t PL);

    /**
     * Convert probabilities to PHRED score.
     */
    uint32_t prob2pl(double x);

    /**
     * Compute log(x)
     */
    double log10(double x);

    /**
     * Compute log(xy)
     */
    double log10prod(double x, double y);

    /**
     * Compute log(x+y)
     */
    double log10sum(double x, double y);

    /**
     * Compute log10 factorial x.
     */
    double log10fact(uint32_t x);

    /**
     * Compute log10 nCr
     */
    double log10choose(uint32_t n, uint32_t r);

    /**
     * Round a value
     */
    double round(double x);
};

#endif