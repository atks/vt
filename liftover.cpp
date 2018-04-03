/* The MIT License

   Copyright (c) 2018 Adrian Tan <atks@umich.edu>

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

#include "liftover.h"

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
    std::string chain_file;


    int32_t compression_level;
    std::string streaming_selection_bed_file;
    uint32_t left_window;
    uint32_t right_window;
    std::vector<GenomeInterval> intervals;
    std::vector<std::string> samples;
    std::string variant;
    uint32_t sort_window_size;
    bool stream_selection;
    bool print_header;
    bool print_header_only;
    bool print_sites_only;
    bool print;
    bool debug;
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
            std::string desc = "Liftover a VCF.GZ or BCF file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_chain_file("c", "c", "chain file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output_vcf file [liftover.<chain_file>]", true, "", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            if (output_vcf_file=="")
            {
                //todo::update file name only incase it is a path.
                output_vcf_file = "liftover." + output_vcf_file;
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
        odw = new BCFOrderedWriter(output_vcf_file, sort_window_size, compression_level);
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
        filter.parse(fexp.c_str(), debug);
        filter_exists = fexp=="" ? false : true;

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

    void liftover()
    {
        bcf1_t *v = NULL;
        bcf_hdr_t *h = odr->hdr;
        Variant variant;


        while (odr->read(v))
        {
//            bcf_print(h,v);

            if (stream_selection)
            {
                bcf_unpack(v, BCF_UN_STR);
                std::string chrom = bcf_get_chrom(odr->hdr,v);
                int32_t start1 = bcf_get_pos1(v);
                int32_t end1 = bcf_get_end1(v);

                
            }

            if (filter_exists)
            {
                vm->classify_variant(h, v, variant);
                if (!filter.apply(h, v, &variant, debug))
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




    };

    void print_options()
    {
        if (!print) return;

        std::clog << "liftover v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "             output VCF file       " << output_vcf_file << "\n";
        std::clog << "             chain file            " << chain_file << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
    };

    ~Igor() {};

    private:
};
}

bool liftover(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.liftover();
    igor.print_stats();
    return igor.print;
};