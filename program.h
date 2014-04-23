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

#ifndef PROGRAM_H
#define PROGRAM_H

#include <vector>
#include <map>
#include <typeinfo>
#include "htslib/hts.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "hts_utils.h"
#include "utils.h"
#include "genome_interval.h"

class VTOutput : public TCLAP::StdOutput
{
	public:

    void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e);

    void usage(TCLAP::CmdLineInterface& c);
};

/**
 * Provides an interface for programs in vt.
 *
 *
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
     * Parse multiple files from command line unlabeled arguments or -L denoted file list.  If both are defined, the files are merged.
     *
     * @files          - file names are stored in this vector
     * @argument_files - vector of input files
     * @file_list      - file names stored in a file
     *
     */
    void parse_files(std::vector<std::string>& files, const std::vector<std::string>& arg_files, std::string file_list);

	/**
     * Parse intervals. Processes the interval list first followed by the interval string. Duplicates are dropped.
     *
     * @intervals       - intervals stored in this vector
     * @interval_list   - file containing intervals
     * @interval_string - comma delimited intervals in a string
     */
	void parse_intervals(std::vector<GenomeInterval>& intervals, std::string interval_list, std::string interval_string);

    /**
     * Parse samples. Processes the sample list. Duplicates are dropped.
     *
     * @nsamples     - number of unique samples found in list
     * @sample_list  - file containing sample names
     */
    char** read_sample_list(int32_t& nsamples, std::string sample_list);

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

 	/**
     * Print reference FASTA file option.
     */
    void print_ref_op(const char* option_line, std::string ref_fasta_file);
    
    /**
     * Print string option, hide if not present.
     */
    void print_str_op(const char* option_line, std::string str_value);

    /**
     * Print input files.
     */
    void print_ifiles(const char* option_line, std::vector<std::string>& files);

  	/**
     * Print intervals option.
     */
    void print_int_op(const char* option_line, std::vector<GenomeInterval>& intervals);

    private:
};

#endif