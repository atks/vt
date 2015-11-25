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

#include <log_tool.h>

/**
 * Round a value
 */
double LogTool::round(double x)
{
    return (x > 0.0) ? std::floor(x + 0.5) : std::ceil(x - 0.5);
};

/**
 * Convert -10log(p) to p.
 */
double LogTool::pl2prob(uint32_t pl)
{
    if (pl>=PL.size())
    {
        //cap
        if (pl > 3236)
        {
            pl = 3236;
        }

        for (uint32_t i=PL.size(); i<=pl; ++i)
        {
            double p = std::pow(10, -((double) i)/10.0);
            PL.push_back(p);
            PL_one_minus_p.push_back(-10*std::log10(1-p));
        }
    }

    return PL[pl];
}

/**
 * Convert -10log(p) to -10log(1-p).
 */
double LogTool::pl2pl_one_minus_p(uint32_t pl)
{
    if (pl>=PL.size())
    {
        //cap
        if (pl > 3236)
        {
            pl = 3236;
        }

        for (uint32_t i=PL.size(); i<=pl; ++i)
        {
            double p = std::pow(10, -((double) i)/10.0);
            PL.push_back(p);
            PL_one_minus_p.push_back(-10*std::log10(1-p));
        }
    }

    return PL_one_minus_p[pl];
}

/**
 * Convert -10log(p) to log10(sqrt(p(1-p))).
 */
double LogTool::pl2log10_varp(uint32_t pl)
{
    if (pl>=LOG10_VARP.size())
    {
        if (pl > 3236)
        {
            pl = 3236;
        }

        for (size_t i=LOG10_VARP.size(); i<=pl; ++i)
        {
            double e = std::pow(10, -((double) i)/10.0);
            double v = log10(e*(1-e))/2;
            LOG10_VARP.push_back(v);
        }
    }

    return LOG10_VARP[pl];
}

/**
 * Convert probabilities to PHRED score.
 */
uint32_t LogTool::prob2pl(double x)
{
    if (x>1 || x<0)
    {
        std::cerr << "[e] x is not a probability\n";
        exit(1);
    }

    return (uint32_t) (round(-10*std::log10(x)));
}

/**
 * Compute log(x)
 */
double LogTool::log10(double x)
{
    return x==0? LOGZERO : std::log10(x);
}

/**
 * Compute log(xy)
 */
double LogTool::log10prod(double x, double y)
{
    if (x==LOGZERO || y==LOGZERO) return LOGZERO;

    return x + y;
}

/**
 * Compute log(x+y).
 */
double LogTool::log10sum(double x, double y)
{
    if (x==LOGZERO) return y;
    if (y==LOGZERO) return x;

    if (x<y)
    {
        x = y-x;
        y -= x;
        x += y;
    }
    else if (x==y)
    {
        return log10(2) + x;
    }

    return x + std::log10(1+pow(10,y-x));
}

/**
 * Compute log10 factorial x.
 */
double LogTool::log10fact(uint32_t x)
{
    if(LOG10FACT.size()==0)
    {
        LOG10FACT.push_back(0);
    }

    if (x>=LOG10FACT.size())
    {
        for (uint32_t i = LOG10FACT.size(); i<=x; ++i)
        {
            LOG10FACT.push_back(LOG10FACT[i-1] + std::log10(i));
        }
    }

    return LOG10FACT[x];
}

/**
 * Compute log10 nCr.
 */
double LogTool::log10choose(uint32_t n, uint32_t r)
{
    assert(r<=n);
    return log10fact(n) - log10fact(r) - log10fact(n-r);
}
