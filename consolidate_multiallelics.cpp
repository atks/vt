/* The MIT License

   Copyright (c) 2016 Adrian Tan <atks@umich.edu>

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

#include "consolidate_multiallelics.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string ref_fasta_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    int32_t window_overlap;
    bool debug;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    ///////
    //i/o//
    ///////
    MultiallelicsConsolidator* mc;

    /////////
    //stats//
    /////////
    int32_t no_vntrs_annotated;

    ////////////////
    //common tools//
    ////////////////
    VariantManip* vm;
    ReferenceSequence* rs;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "consolidates multiallelic variants from biallelic variants represented with their sequence explicitly.  Requires END annotation.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_window_overlap("w", "w", "window overlap for merging variants [0]", false, 0, "integer", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            window_overlap = arg_window_overlap.getValue();
            fexp = arg_fexp.getValue();
            debug = arg_debug.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
            abort();
        }
    };

    ~Igor()
    {
    };

    void initialize()
    {
        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        //////////////////////
        //i/o initialization//
        //////////////////////
        mc = new MultiallelicsConsolidator(input_vcf_file, intervals, output_vcf_file, fexp, ref_fasta_file, window_overlap);

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void print_options()
    {
        std::clog << "consolidate_multiallelics v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)                    " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file                      " << output_vcf_file << "\n";
        print_boo_op("         [d] debug                                ", debug);
        print_num_op("         [w] window overlap for merging variants  ", window_overlap);
        print_ref_op("         [r] ref FASTA file                       ", ref_fasta_file);
        print_int_op("         [i] intervals                            ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants              " << mc->no_variants << "\n";
        std::cerr << "       no. of multiallelics added   " << mc->no_added_multiallelics << "\n";
        std::clog << "\n";
    }

    void consolidate_multiallelics()
    {
        mc->process();
    };

    private:
};
}

void consolidate_multiallelics(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.consolidate_multiallelics();
    igor.print_stats();
};
