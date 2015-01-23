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

#include "seq.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::vector<GenomeInterval> intervals;
    std::string ref_fasta_file;
    bool print;

    ///////
    //i/o//
    ///////

    /////////
    //stats//
    /////////

    /////////
    //tools//
    /////////
    faidx_t *fai;
    bool reference_present;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "extract sequences from a fasta file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_quiet("q", "q", "quiet [false]", cmd, false);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_fasta_file = arg_ref_fasta_file.getValue();
            print = !arg_quiet.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////

        ////////////////////////
        //stats initialization//
        ////////////////////////

        ////////////////////////
        //tools initialization//
        ////////////////////////
        if (ref_fasta_file!="")
        {
            fai = fai_load(ref_fasta_file.c_str());
            if (fai==NULL)
            {
                fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
                exit(1);
            }
            reference_present = (fai!=NULL);
        }
    }

    void seq()
    {
        if (print) std::clog << "\n";
        for (size_t i=0; i<intervals.size();++i)
        {
            char* ref = 0;
            int32_t ref_len = 0;
            ref = faidx_fetch_uc_seq(fai, intervals[i].seq.c_str(), intervals[i].start1-1, intervals[i].end1-1, &ref_len);

            if (ref_len!=0)
            {
                if (!print)
                {
                    std::cout << ref;
                    free(ref);
                }
                else
                {
                    std::clog << "interval " << intervals[i].to_string() << "\n";
                    std::clog << "sequence " << ref << "\n";
                    free(ref);
                }
            }
            else
            {
                if (!print)
                {
                    std::cout << "X";
                }
                else
                {
                    std::clog << "interval " << intervals[i].to_string() << "\n";
                    std::clog << "sequence : cannot be read\n";
                }
                
            }
        }
        if (print) std::clog << "\n";
    };

    void print_options()
    {
        if (!print) return;
        
        std::clog << "seq v" << version << "\n";
        std::clog << "\n";
        std::clog << "options: [r] reference FASTA file  " << ref_fasta_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
    };

    ~Igor() {};

    private:
};

}

bool seq(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.seq();
    igor.print_stats();
    return igor.print;
};
