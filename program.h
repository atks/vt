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

#ifndef PROGRAM_H
#define PROGRAM_H

#include <vector>
#include <map>
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/join.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "htslib/hts.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"

class VTOutput : public TCLAP::StdOutput
{
	public:

    void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e);
    
    void usage(TCLAP::CmdLineInterface& c);        
};

/**
 * Provides an interface 
 *
 * @intervals       - intervals stored in this vector
 * @interval_list   - file containing intervals
 * @interval_string - comma delimited intervals in a string 
 */ 
class Program
{
    public:

    std::string version;
    
    /**
     * Process arguments.
     */
    Program(){};

	/**
     * Parse intervals. Processes the interval list first followed by the interval string. Duplicates are dropped.
     *
     * @intervals       - intervals stored in this vector
     * @interval_list   - file containing intervals
     * @interval_string - comma delimited intervals in a string 
     */ 
	void parse_intervals(std::vector<std::string>& intervals, std::string interval_list, std::string interval_string);

    /**
     * Initialize I/O and shared objects.
     */
    void initialize(){};
    
    /**
     * Print options.
     */
    void print_options(){};
    
    /**
     * Print run stats.
     */
    void print_stats(){};
 	
    private:
};
    
#endif