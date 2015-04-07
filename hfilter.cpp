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

#include "hfilter.h"

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
    std::string regions_bed_file;
    std::string filter_tag;
    std::string filter_tag_desc;
    bool clear_filter;

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
    int32_t no_variants_filtered;
    int32_t no_variants;

    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "filters variants in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_filter_tag("t", "t", "filter tag []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_filter_tag_desc("d", "d", "filter tag description []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", true, "", "str", cmd);
            TCLAP::SwitchArg arg_clear_filter("x", "x", "clear filter [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            clear_filter = arg_clear_filter.getValue();
            filter_tag = arg_filter_tag.getValue();
            filter_tag_desc = arg_filter_tag_desc.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
            abort();
        }
    };

    ~Igor() {};

    void initialize()
    {
        //******************
        //i/o initialization
        //******************
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);

        std::string hrec = "##FILTER=<ID=" + filter_tag + ",Description=\"" + filter_tag_desc + "\">";
        bcf_hdr_append(odw->hdr, hrec.c_str());

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip();

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants_filtered = 0;
        no_variants = 0;
    }

    void print_options()
    {
        std::clog << "hfilter v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)     " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_str_op("         [f] filter expression     ", fexp);
        print_str_op("         [t] filter tag            ", filter_tag);
        print_str_op("         [d] filter description    ", filter_tag_desc);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants filtered     " << no_variants_filtered << "\n";
        std::cerr << "       total no. of variants        " << no_variants << "\n";std::clog << "\n";
    }

    void hfilter()
    {
        odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odr->hdr;
        Variant variant;
        int32_t filter_id = bcf_hdr_id2int(h, BCF_DT_ID, filter_tag.c_str());
        while (odr->read(v))
        {
            vm->classify_variant(h, v, variant);
            if (filter.apply(h, v, &variant, false))
            {
                if (clear_filter)
                {
                    bcf_update_filter(h, v, &filter_id, 1);
                } 
                else
                {
                    bcf_add_filter(h, v, filter_id);
                }
                
                ++no_variants_filtered;
            }
            
            ++no_variants;
            odw->write(v);
            v = odw->get_bcf1_from_pool();
        }

        odw->close();
    };

    private:
};
}

void hfilter(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.hfilter();
    igor.print_stats();
};