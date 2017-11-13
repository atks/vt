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
    std::string regions_file;
    std::string REGIONS_TAG;
    std::string REGIONS_TAG_DESC;
    std::string REGIONS_LEFT_TAG;
    std::string REGIONS_LEFT_TAG_DESC;
    std::string REGIONS_RIGHT_TAG;
    std::string REGIONS_RIGHT_TAG_DESC;
    uint32_t left_window; 
    uint32_t right_window;
    bool use_bed;
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    //////////
    //filter//
    //////////
    std::vector<std::string> fexps;
    Filter filter;
    bool filter_exists;
    
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
    OrderedBCFOverlapMatcher *obom_regions;

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
            TCLAP::ValueArg<std::string> arg_regions_file("b", "b", "regions BED/BCF file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_REGIONS_TAG("t", "t", "regions tag []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_REGIONS_TAG_DESC("d", "d", "regions tag description []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_left_window("l", "l", "left window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_right_window("r", "r", "right window size for overlap []", false, 0, "int", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_filters(fexps, arg_fexp.getValue(), 2, false);
            regions_file = arg_regions_file.getValue();
            left_window = arg_left_window.getValue();
            right_window = arg_right_window.getValue();
            REGIONS_TAG = arg_REGIONS_TAG.getValue();
            REGIONS_TAG_DESC = arg_REGIONS_TAG_DESC.getValue();
            REGIONS_LEFT_TAG = REGIONS_TAG + "_LEFT";        
            REGIONS_RIGHT_TAG = REGIONS_TAG + "_RIGHT";        
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
        ////////////////////
        //i/o initialization
        ////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);

        std::string hrec = "##INFO=<ID=" + REGIONS_TAG + ",Number=0,Type=Flag,Description=\"" + REGIONS_TAG_DESC + "\">";
        bcf_hdr_append(odw->hdr, hrec.c_str());
        if (left_window)
        {
            std::string hrec = "##INFO=<ID=" + REGIONS_LEFT_TAG + ",Number=0,Type=Flag,Description=\"" + REGIONS_TAG_DESC + " (Left window)\">";
            bcf_hdr_append(odw->hdr, hrec.c_str());
        }
        if (right_window)
        {
            std::string hrec = "##INFO=<ID=" + REGIONS_RIGHT_TAG + ",Number=0,Type=Flag,Description=\"" + REGIONS_TAG_DESC + " (Right window)\">";
            bcf_hdr_append(odw->hdr, hrec.c_str());
        }
        
        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexps[0].c_str(), false);
        filter_exists = fexps[0]=="" ? false : true;
            
        ///////////////////////
        //tool initialization//
        ///////////////////////
        if (str_ends_with(regions_file, ".bed") || str_ends_with(regions_file, ".bed.gz"))
        {    
            use_bed = true;
            orom_regions = new OrderedRegionOverlapMatcher(regions_file);
        }
        else if (str_ends_with(regions_file, ".vcf") || str_ends_with(regions_file, ".vcf.gz") ||  str_ends_with(regions_file, ".bcf"))
        {    
            use_bed = false;
            obom_regions = new OrderedBCFOverlapMatcher(regions_file, intervals, fexps[1]);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Need to at least specify either a bed or bcf file\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
        
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
        std::clog << "options:     input VCF file(s)       " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file         " << output_vcf_file << "\n";
        print_str_op("         [f] filter 1                ", fexps[0]);
        print_str_op("             filter 2                ", fexps[1]);
        print_str_op("         [t] region INFO tag         ", REGIONS_TAG);    
        print_str_op("         [b] region INFO description ", REGIONS_TAG_DESC);
        print_str_op("         [m] regions file            ", regions_file);
        print_num_op("         [l] left window             ", left_window);
        print_num_op("         [r] right window            ", right_window);
        print_int_op("         [i] intervals               ", intervals);
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

        bcf_hdr_t *h = odr->hdr;
        bcf1_t *v = bcf_init1();
        std::vector<Interval*> overlaps;
        Variant variant;
        kstring_t s = {0,0,0};
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            if (filter_exists)
            {
                vm->classify_variant(h, v, variant);
                if (!filter.apply(h, v, &variant, false))
                {
                    continue;
                }
            }
           
            std::string chrom = bcf_get_chrom(odr->hdr,v);
            int32_t beg1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end1(v);

            if (use_bed)
            {
                if (orom_regions->overlaps_with(chrom, beg1-left_window, end1+right_window))
                {
                    if (left_window+right_window)
                    {    
                        std::vector<Interval>& regs = orom_regions->overlapping_regions;
                        for (int32_t i=0; i<regs.size(); ++i)
                        {
                            if (beg1>=regs[i].beg1-left_window && beg1<=regs[i].beg1+right_window)
                            {
                                bcf_update_info_flag(odr->hdr, v, REGIONS_LEFT_TAG.c_str(), "", 1);
                            }    
                            else if (end1>=regs[i].end1-right_window && end1<=regs[i].end1+right_window)
                            {
                                bcf_update_info_flag(odr->hdr, v, REGIONS_RIGHT_TAG.c_str(), "", 1);
                            }
                        }
                    }                    
                    
                    bcf_update_info_flag(odr->hdr, v, REGIONS_TAG.c_str(), "", 1);
                    ++no_variants_annotated;
                }
            }
            else
            {
                if (obom_regions->overlaps_with(chrom, beg1-left_window, end1+right_window))
                {
                    bcf_update_info_flag(odr->hdr, v, REGIONS_TAG.c_str(), "", 1);
                    ++no_variants_annotated;
                }
            }

            ++no_variants;
            odw->write(v);
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