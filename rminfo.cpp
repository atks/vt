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

#include "rminfo.h"

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
    std::vector<std::string> info_tags;
    bool print;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    uint32_t no_variants_with_removed_info;

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
            std::string desc = "removes INFO tags from a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_info_tags("t", "t", "list of info tags to be removed []", true, "", "str", cmd);
            TCLAP::SwitchArg arg_quiet("q", "q", "do not print options and summary [false]", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_string_list(info_tags, arg_info_tags.getValue());
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
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);
        odw->write_hdr();

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_variants_with_removed_info = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void rminfo()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odr->hdr;

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            int32_t ret = 0;
            for (uint32_t i=0; i<info_tags.size(); ++i)
            {
                ret += bcf_update_info(h, v, info_tags[i].c_str(), NULL, 0, 0);
            }

            //todo: this is not correct, ret only returns non 0 upon an error.  
            if (!ret) ++no_variants_with_removed_info;

            ++no_variants;

            odw->write(v);
            v = odw->get_bcf1_from_pool();
        }

        odw->close();
        odr->close();
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "rminfo v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file                                  " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file                                 " << output_vcf_file << "\n";
        std::clog << "         [q] quiet                                           " << (!print ? "true" : "false") << "\n";
        print_int_op("         [i] intervals                                       ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::clog << "stats: total no. variants                   : " << no_variants << "\n";
        std::clog << "       total no. variants with removed info : " << no_variants_with_removed_info << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool rminfo(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.rminfo();
    igor.print_stats();

    return igor.print;
};
