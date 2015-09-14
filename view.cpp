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

#include "view.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::vector<std::string> samples;
    std::string variant;
    uint32_t sort_window_size;
    bool print_header;
    bool print_header_only;
    bool print_sites_only;
    bool print;
    int32_t no_subset_samples;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    uint32_t no_samples;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Views a VCF or BCF or VCF.GZ file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::SwitchArg arg_print_header("h", "h", "omit header, this option is honored only for STDOUT [false]", cmd, false);
            TCLAP::SwitchArg arg_print_header_only("H", "H", "print header only, this option is honored only for STDOUT [false]", cmd, false);
            TCLAP::SwitchArg arg_print_sites_only("s", "s", "print site information only without genotypes [false]", cmd, false);
            TCLAP::ValueArg<uint32_t> arg_sort_window_size("w", "w", "local sorting window size [0]", false, 0, "int", cmd);
            //TCLAP::ValueArg<std::string> arg_sample_list("s", "s", "file containing list of sample []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF/VCF.GZ/BCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            print_header = arg_print_header.getValue();
            print_header_only = arg_print_header_only.getValue();
            no_subset_samples = arg_print_sites_only.getValue() ? 0 : -1;
            print = arg_print.getValue();
            sort_window_size = arg_sort_window_size.getValue();
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
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, sort_window_size);
        if (no_subset_samples==-1)
        {
            odw->link_hdr(odr->hdr);
        }
        else if (no_subset_samples==0)
        {
            odw->link_hdr(bcf_hdr_subset(odr->hdr, 0, 0, 0));
        }

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

//        exit(1);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_samples = 0;

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");
    }

    void view()
    {
        if (print_header_only && output_vcf_file == "-")
        {
            odw->write_hdr();
            odr->close();
            odw->close();
            return;
        }

        if (print_header || output_vcf_file != "-") odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odr->hdr;
        Variant variant;
        while (odr->read(v))
        {
//            bcf_print(h,v);
            
            if (filter_exists)
            {   
                vm->classify_variant(h, v, variant);
                if (!filter.apply(h, v, &variant, false))
                {
                    continue;
                }
            }

            if (no_subset_samples==0)
            {
                bcf_subset(odw->hdr, v, 0, 0);
                //maybe add some additional adhoc fixing for BCF files that do not have a complete header.
            }
            odw->write(v);
            if (sort_window_size)
            {
                v = odw->get_bcf1_from_pool();
            }
            ++no_variants;
        }

        odw->close();
        odr->close();
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "view v" << version << "\n\n";

        std::clog << "options:     input VCF file              " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file             " << output_vcf_file << "\n";
        std::clog << "         [w] sort window size            " << sort_window_size << "\n";
        std::clog << "         [h] print header                " << (print_header ? "yes" : "no") << "\n";
        std::clog << "         [H] print header only           " << (print_header_only ? "yes" : "no") << "\n";
        std::clog << "         [s] print site information only " << (print_sites_only ? "yes" : "no") << "\n";
        std::clog << "         [p] print options and stats     " << (print ? "yes" : "no") << "\n";
        print_str_op("         [f] filter                      ", fexp);
        print_int_op("         [i] intervals                   ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::clog << "stats: no. variants  : " << no_variants << "\n";
        std::clog << "       no. samples   : " << no_samples << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool view(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.view();
    igor.print_stats();
    return igor.print;
};
