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

#include <iostream>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <map>
#include <queue>
#include <time.h>
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "normalize.h"
#include "merge_duplicate_variants.h"

void print_time(double t)
{
    if (t<60)
    {
        std::clog << "Time elapsed: " << t << "s\n\n";
   	}
   	else if (t<60*60) //less than an hour
   	{
   	    std::clog << "Time elapsed: " << ((int32_t)(t/60)) << "m " << fmod(t,60) << "s\n\n";
   	}
   	else if (t<60*60*24) //less than a day
   	{
   	    std::clog << "Time elapsed: " << ((int32_t)(t/(60*60))) << "h ";
   	    t = fmod(t,60*60); //remaining minutes
   	    std::clog <<                     ((int32_t)(t/(60))) << "m " << fmod(t,60) << "s\n\n";
   	}
   	else if (t<60*60*24*365) //less than a year
   	{
   	    std::clog << "Time elapsed: " << ((int32_t)(t/(60*60*24))) << "d ";
   	    t = fmod(t,60*60*24); //remaining hours
   	    std::clog <<                     ((int32_t)(t/(60*60))) << "h ";
   	    t = fmod(t,60*60); //remaining minutes
   	    std::clog <<                      ((int32_t)(t/60)) << "m " << fmod(t,60) << "s\n\n";
   	}
};

void help()
{
	std::clog << "Help page on http://statgen.sph.umich.edu/wiki/Vt\n\n";

	std::clog << "normalize                 normalize indels\n";
	std::clog << "merge_duplicate_variants  merge duplicate variants\n";
	std::clog << "\n";
}

int main(int argc, char ** argv)
{
    time_t t0;
	std::time(&t0);

    std::clog << "\n=======\n";
	std::clog << "vt v0.5\n";
	std::clog << "=======\n\n";

    if (argc==1)
    {
        help();
		exit(0);
    }

    std::string cmd(argv[1]);
    
    if (argc>1 && cmd=="normalize")
	{
	    normalize(argc-1, ++argv);
	}
	else if (argc>1 && cmd=="merge_duplicate_variants")
	{
	    merge_duplicate_variants(argc-1, ++argv);
	}	
	else
    {
        std::clog << "Command not found: " << argv[1] << "\n\n";
        help();
        exit(0);
    }

	time_t t1;
	std::time(&t1);
    
    print_time(difftime(t1,t0));

    return 0;
}
