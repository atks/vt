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

#include "consolidate_vntrs.h"

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
    std::string ref_fasta_file;
    std::vector<GenomeInterval> intervals;
    bool debug;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    VNTRConsolidator *vc;

    Igor(int argc, char **argv)
    {
        version = "0.57";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Consolidates VNTRs by\n"
                 "              1. removing overlapping VNTRs and leaving behind the most complete VNTR\n"
                 "              2. resolves adjacent VNTRs\n"
                 "              3. marking VNTRs that are \n"
                 "              4. Adds INFO fields indicating the number of variants that overlap with this variant";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
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
        ////////////////////////
        //tools initialization//
        ////////////////////////
        vc = new VNTRConsolidator(input_vcf_file, intervals, output_vcf_file, ref_fasta_file);
        vm = new VariantManip();
    }

    void consolidate_vntrs()
    {
        bcf1_t *v = vc->odw->get_bcf1_from_pool();

        Variant* variant;
        while (vc->odr->read(v))
        {
            variant = new Variant(vc->odw->hdr, v);
            vc->flush_variant_buffer(variant);
            vc->insert_variant_record_into_buffer(variant);
            v = vc->odw->get_bcf1_from_pool();

            ++vc->no_total_variants;
        }

        vc->flush_variant_buffer();
        vc->close();
    };

    void print_options()
    {
        std::clog << "consolidate_vntrs v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals              ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "       VNTR Classification          \n";
        std::clog << "       Number of VNTRs      " << vc->no_vntrs << "\n";
        std::clog << "\n";
        std::clog << "       Number of isolated exact VNTRs      " << vc->no_isolated_exact_vntrs << "\n";
        std::clog << "                     perfect concordance   " << vc->no_perfect_concordance_isolated_exact_vntrs << "\n";
        std::clog << "                     imperfect concordance " << vc->no_imperfect_concordance_isolated_exact_vntrs << "\n";
        std::clog << "       Number of isolated inexact VNTRs    " << vc->no_isolated_inexact_vntrs << "\n";
        std::clog << "                     complete overlaps     " << vc->no_isolated_complete_overlap_vntrs << "\n";
        std::clog << "                     partial overlaps      " << vc->no_isolated_partial_overlap_vntrs << "\n";
        std::clog << "                     no overlaps           " << vc->no_isolated_no_overlap_vntrs << "\n";
        std::clog << "       Number of clustered exact VNTRs     " << vc->no_clustered_exact_vntrs << "\n";

        for (uint32_t i=0; i<vc->overlapping_vntr_hist.size(); ++i)
        {
            if (vc->overlapping_vntr_hist[i])
            {
                std::clog << "       " << i << "    "  << vc->overlapping_vntr_hist[i] << "\n";
            }
        }


        std::clog << "       Number of clustered inexact VNTRs   " << vc->no_clustered_inexact_vntrs << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void consolidate_vntrs(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.consolidate_vntrs();
    igor.print_stats();
};
