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

#include "annotate_variants.h"

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
    std::string mdust_bed_file;
    std::string gencode_cds_bed_file;
    bool annotate_low_complexity_regions;
    bool annotate_coding;

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
    int32_t no_variants_annotated;

    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;
    OrderedRegionOverlapMatcher *orom_lcplx;
    OrderedRegionOverlapMatcher *orom_gencode_cds;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "annotates variants in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_mdust_bed_file("m", "m", "mdust Low Complexity Regions BED file []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_gencode_cds_bed_file("g", "g", "GENCODE Coding Regions BED file []", false, "", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            ref_fasta_file   = arg_ref_fasta_file.getValue();
            mdust_bed_file = arg_mdust_bed_file.getValue();
            annotate_low_complexity_regions = mdust_bed_file != "" ? true : false;
            gencode_cds_bed_file = arg_gencode_cds_bed_file.getValue();
            annotate_coding = gencode_cds_bed_file != "" ? true : false;
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
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant Type - SNP, MNP, INDEL, CLUMPED\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit in a STR or Homopolymer\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Repeat Length\">");

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip(ref_fasta_file);

        if (annotate_low_complexity_regions)
        {
            bcf_hdr_append(odw->hdr, "##INFO=<ID=LOW_COMPLEXITY,Number=0,Type=Flag,Description=\"INDEL is in low complexity region\">");
            orom_lcplx = new OrderedRegionOverlapMatcher(mdust_bed_file);
        }

        if (annotate_coding)
        {
            bcf_hdr_append(odw->hdr, "##INFO=<ID=GENCODE_FS,Number=0,Type=Flag,Description=\"Frameshift INDEL\">");
            bcf_hdr_append(odw->hdr, "##INFO=<ID=GENCODE_NFS,Number=0,Type=Flag,Description=\"Non Frameshift INDEL\">");
            orom_gencode_cds = new OrderedRegionOverlapMatcher(gencode_cds_bed_file);
        }

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants_annotated = 0;
    }

    void print_options()
    {
        std::clog << "annotate_indels v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)     " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_str_op("         [m] mdust bed file        ", mdust_bed_file);
        print_str_op("         [g] gencode cds bed file  ", gencode_cds_bed_file);
        print_ref_op("         [r] ref FASTA file        ", ref_fasta_file);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants annotated     " << no_variants_annotated << "\n";
        std::clog << "\n";
    }

    void annotate_variants()
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

            if (filter_exists)
            {
                if (!filter.apply(odr->hdr, v, &variant, false))
                {
                    continue;
                }
            }

            vm->vtype2string(vtype, &s);
            if (s.l)
            {
                bcf_update_info_string(odr->hdr, v, "VT", s.s);
            }

            std::string chrom = bcf_get_chrom(odr->hdr,v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);

            if (annotate_low_complexity_regions)
            {
                if (orom_lcplx->overlaps_with(chrom, start1, end1))
                {
                    bcf_update_info_flag(odr->hdr, v, "LOW_COMPLEXITY", "", 1);
                }
            }

            if (vtype==VT_SNP)
            {
                //synonymous and non synonymous annotation

            }
            else if (vtype&VT_INDEL)
            {
                if (annotate_coding)
                {
                    bool overlap = false;
                    if ((overlap = orom_gencode_cds->overlaps_with(chrom, start1, end1)))
                    {
                        if (abs(variant.alleles[0].dlen)%3!=0)
                        {
                            bcf_update_info_flag(odr->hdr, v, "GENCODE_FS", "", 1);
                        }
                        else
                        {
                            bcf_update_info_flag(odr->hdr, v, "GENCODE_NFS", "", 1);
                        }

                    }
                }
            }

            ++no_variants_annotated;
            odw->write(v);
            v = odw->get_bcf1_from_pool();
        }

        odw->close();
    };

    private:
};
}

void annotate_variants(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.annotate_variants();
    igor.print_stats();
};
