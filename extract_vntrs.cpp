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

#include "extract_vntrs.h"

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
    bool debug;
    std::string vntr_classification;
    int32_t vntr_classification_code;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    ///////
    //i/o//
    ///////
    VNTRExtractor* ve;

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
            std::string desc = "extracts vntrs from indels annotated with annotate_indels.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_vntr_classification("c", "c", "classification schemas of tandem repeat [tan_kang2015]\n"
                 "              1 : lai2003      \n"
                 "              2 : kelkar2008   \n"
                 "              3 : fondon2012   \n"
                 "              4 : ananda2013   \n"
                 "              5 : willems2014  \n"
                 "              6 : tan_kang2015 \n"
                 "              7 : exact_vntr   \n"
                 "              8 : fuzzy_vntr   \n",
                 false, "tan_kang2015", "string", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            vntr_classification = arg_vntr_classification.getValue();
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
        ///////////
        //options//
        ///////////
        if (vntr_classification == "lai2003")
        {
            vntr_classification_code = LAI2003;
        }
        else if (vntr_classification == "kelkar2008")
        {
            vntr_classification_code = KELKAR2008;
        }
        else if (vntr_classification == "fondon2012")
        {
            vntr_classification_code = FONDON2012;
        }
        else if (vntr_classification == "ananda2013")
        {
            vntr_classification_code = ANANDA2013;
        }
        else if (vntr_classification == "willems2014")
        {
            vntr_classification_code = WILLEMS2014;
        }
        else if (vntr_classification == "tan_kang2015")
        {
            vntr_classification_code = TAN_KANG2015;
        }
        else if (vntr_classification == "exact_vntr")
        {
            vntr_classification_code = EXACT_VNTR;
        }
        else if (vntr_classification == "fuzzy_vntr")
        {
            vntr_classification_code = FUZZY_VNTR;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] VNTR classification not recognized: %s\n", __FILE__,__LINE__,__FUNCTION__, vntr_classification.c_str());
            exit(1);
        }

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        //////////////////////
        //i/o initialization//
        //////////////////////
         ve = new VNTRExtractor(input_vcf_file, intervals, output_vcf_file, fexp,  vntr_classification_code, ref_fasta_file);

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void print_options()
    {
        std::clog << "extract_vntrs v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file          " << output_vcf_file << "\n";
        print_boo_op("         [d] debug                    ", debug);
        print_str_op("         [c] VNTR classification      ", vntr_classification);
        print_ref_op("         [r] ref FASTA file           ", ref_fasta_file);
        print_int_op("         [i] intervals                ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of VNTRs added   " << no_vntrs_annotated << "\n";
        std::clog << "\n";
    }

    void extract_vntrs()
    {
        ve->process();
    };

    private:
};
}

void extract_vntrs(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.extract_vntrs();
    igor.print_stats();
};
