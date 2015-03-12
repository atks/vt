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

#include <binomial_distribution.h>


/**
 * Constructor.
 */
BinomialDistribution::BinomialDistribution()
{
    this->p = -1;
};

/**
 * Constructor.
 */
BinomialDistribution::BinomialDistribution(float p)
{
    set_p(p);
};

/**
 * Sets p.
 */
void BinomialDistribution::set_p(float p)
{
    if (this->p!=p)
    {
        this->p = p;
        pvalues.clear();
    }
}

/**
 * Get P(X<=x) where X~Binomial(0.5, n)
 */
float BinomialDistribution::get_pvalue(uint32_t x, uint32_t n)
{
    if (x<=n)
    {
        uint32_t current_size = pvalues.size();
        if (current_size>n)
        {
            return pvalues[n][x];
        }
        else
        {
            pvalues.resize(n+1);
            for (uint32_t i=current_size; i<=n; ++i)
            {
                for (uint32_t j=0; j<=i; ++j)
                {
                    //std::cerr << "x=" << j << " n=" << i << " p=" << p << " = " << pbinom(j, i, p, 1, 0) << "\n";
                    
                    //double pbinom(double x, double n, double p, int lower_tail, int log_p)
                    pvalues[i].push_back(pbinom(j, i, p, 1, 0));
                    
                    //double dbinom(double x, double n, double p, int give_log)
                    //pvalues[i].push_back(dbinom(j, i, p, 1));
                    
                }
            }

            return pvalues[n][x];
        }
    }

    return 0;
}
