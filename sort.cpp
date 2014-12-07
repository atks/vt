/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#include "sort.h"

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
    std::string sort_mode;
    bool print;

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
            std::string desc = "Sorts a VCF or BCF or VCF.GZ file.\n";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::ValueArg<uint32_t> arg_sort_window_size("w", "w", "local sorting window size [0]", false, 0, "int", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF/VCF.GZ/BCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_sort_mode("m", "m", ""
                               "sorting modes [full]\n"
                 "              local : locally sort within a window.\n"
                 "              chrom : sort chromosomes based on order of contigs in header.\n"
                 "                      input must be an indexed vcf.gz\n"
                 "              full  : full sort with no assumptions",
                                 false, "full", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            sort_mode = arg_sort_mode.getValue();
            print = arg_print.getValue();
            sort_window_size = arg_sort_window_size.getValue();
            
            if (sort_mode=="local")
            {
                sort_window_size = sort_window_size ? sort_window_size : 1000;
            }
            else if (sort_mode=="chrom")
            {
                
            }
            else if (sort_mode=="full")
            {
                
            }
            else
            {
                fprintf(stderr, "[%s:%d %s] Sort modecan only be full, local or chrom.\n", __FILE__,__LINE__,__FUNCTION__);
            }
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
                
        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        
        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");
    }

    void sort()
    {
        if (sort_mode=="local")
        {
            odw = new BCFOrderedWriter(output_vcf_file, sort_window_size);
            odw->link_hdr(odr->hdr);
            odw->write_hdr();
            
            bcf1_t *v = odw->get_bcf1_from_pool();
            bcf_hdr_t *h = odr->hdr;
            Variant variant;
            
            while (odr->read(v))
            {
                if (filter_exists)
                {
                    vm->classify_variant(h, v, variant);
                    if (!filter.apply(h, v, &variant, false))
                    {
                        continue;
                    }
                }
    
                odw->write(v);
                v = odw->get_bcf1_from_pool();
                ++no_variants;
            
            }
            
            odw->close();
            odr->close();
        }
        else if (sort_mode=="chrom")
        {
            
        }
        else if (sort_mode=="full")
        {
            
        }
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "view v" << version << "\n\n";

        std::clog << "options:     input VCF file              " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file             " << output_vcf_file << "\n";
        std::clog << "         [w] sort window size            " << sort_window_size << "\n";
        std::clog << "         [m] sorting mode                " << sort_mode << "\n";
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
        std::clog << "       no. samples   : " <<  "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool sort(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.sort();
    igor.print_stats();
    return igor.print;
};
