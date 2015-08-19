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

#include "annotate_regions.h"

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
    std::string REGIONS_TAG;
    std::string REGIONS_TAG_DESC;
    uint32_t left_window; 
    uint32_t right_window;
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    int32_t no_variants_annotated;
    int32_t no_variants;

    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;
    OrderedRegionOverlapMatcher *orom_regions;
    OrderedRegionOverlapMatcher *orom_gencode_cds;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "annotates regions in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_regions_bed_file("b", "b", "regions BED file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_REGIONS_TAG("t", "t", "regions tag []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_REGIONS_TAG_DESC("d", "d", "regions tag description []", true, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_left_window("l", "l", "left window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_right_window("r", "r", "right window size for overlap []", false, 0, "int", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            regions_bed_file = arg_regions_bed_file.getValue();
            left_window = arg_left_window.getValue();
            right_window = arg_right_window.getValue();
            REGIONS_TAG = arg_REGIONS_TAG.getValue();
            REGIONS_TAG_DESC = arg_REGIONS_TAG_DESC.getValue();
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

        std::string hrec = "##INFO=<ID=" + REGIONS_TAG + ",Number=0,Type=Flag,Description=\"" + REGIONS_TAG_DESC + "\">";
        bcf_hdr_append(odw->hdr, hrec.c_str());

        ///////////////////////
        //tool initialization//
        ///////////////////////
        orom_regions = new OrderedRegionOverlapMatcher(regions_bed_file);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants_annotated = 0;
        no_variants = 0;
    }

    void print_options()
    {
        std::clog << "annotate_regions v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)     " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_str_op("         [m] regions bed file      ", regions_bed_file);
        print_num_op("         [l] left window           ", left_window);
        print_num_op("         [r] right window          ", right_window);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants annotated     " << no_variants_annotated << "\n";
        std::cerr << "       total no. of variants         " << no_variants << "\n";std::clog << "\n";
    }

    void annotate_regions()
    {
        odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        std::vector<Interval*> overlaps;
        Variant variant;
        kstring_t s = {0,0,0};
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            std::string chrom = bcf_get_chrom(odr->hdr,v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);

            if (orom_regions->overlaps_with(chrom, start1-left_window, end1+right_window))
            {
                bcf_update_info_flag(odr->hdr, v, REGIONS_TAG.c_str(), "", 1);
                ++no_variants_annotated;
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

void annotate_regions(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.annotate_regions();
    igor.print_stats();
};