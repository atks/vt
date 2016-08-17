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

#include "info2tab.h"

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
    std::string ref_fasta_file;
    std::vector<std::string> info_tags;
    bool debug;

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
    uint32_t new_no_variants;
    uint32_t no_biallelic;
    uint32_t no_multiallelic;
    uint32_t no_additional_biallelic;

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
            std::string desc = "extracts info fields to a tab delimited file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_info_tags("t", "t", "list of info tags to be removed []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::SwitchArg arg_smart("s", "s", "smart decomposition [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            fexp = arg_fexp.getValue();
            parse_string_list(info_tags, arg_info_tags.getValue());
            debug = arg_debug.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
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

        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);
        //bcf_hdr_append_info_with_backup_naming(odw->hdr, "OLD_MULTIALLELIC", "1", "String", "Original chr:pos:ref:alt encoding", false);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description=\"Original chr:pos:ref:alt encoding\">\n");
        odw->write_hdr();

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_biallelic = 0;
        no_multiallelic = 0;

        no_additional_biallelic = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip();
    }

    void info2tab()
    {
        bcf_hdr_t* h = odr->hdr;
        bcf1_t* v = bcf_init();
        Variant variant;

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            if (filter_exists)
            {
                vm->classify_variant(h, v, variant);
                if (!filter.apply(h, v, &variant, false))
                {
                    continue;
                }
            }


            int32_t ret = 0;
            for (uint32_t i=0; i<info_tags.size(); ++i)
            {
                ret += bcf_update_info(h, v, info_tags[i].c_str(), NULL, 0, 0);
            }

            ++no_variants;

            odw->write(v);
            v = odw->get_bcf1_from_pool();
        }

        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "info2tab v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output text file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. variants                 : " << no_variants << "\n";
        std::clog << "       no. biallelic variants       : " << no_biallelic << "\n";
        std::clog << "       no. multiallelic variants    : " << no_multiallelic << "\n";
        std::clog << "\n";
        std::clog << "       no. additional biallelics    : " << no_additional_biallelic << "\n";
        std::clog << "       total no. of biallelics      : " << no_additional_biallelic + no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void info2tab(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.info2tab();
    igor.print_stats();
};
