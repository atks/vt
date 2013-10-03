/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>

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
 * Convert PLs to probabilities.
 */
double LogTool::pl2prob(uint32_t PL)
{
	if (PL>=PLs.size())
	{
		if (PL > 3236)
			PL = 3236;
					
		for (uint32_t i=PLs.size(); i<=PL; ++i)
		{
			PLs.push_back(pow(10, -((double) i)/10.0));
		}
	}
	
	return PLs[PL];
}

/**
 * Compute log(x)
 */	
double LogTool::elog10(double x)
{
	return x==0? LOGZERO : log10(x); 
}

/**
 * Compute log(xy)
 */
double LogTool::elog10prod(double x, double y)
{
	if (x==LOGZERO || y==LOGZERO) return LOGZERO;
	
	return x + y;
}

/**
 * Compute log(x+y)
 */
double LogTool::elog10sum(double x, double y)
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
	
	return x + log10(1+pow(10,y-x));
}
