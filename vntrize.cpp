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

#include "vntrize.h"

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
    std::string ref_vntr_vcf_file;
    int32_t window_size;
    uint32_t left_window;
    uint32_t right_window;
    bool print;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    uint32_t no_variants_vntrized;

    //for allele counts [no indel alleles][no tandem repeat alleles]
    std::vector<std::vector<int32_t> > joint_allele_dist;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    OrderedBCFOverlapMatcher *orom_vntrs;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "updates variants with a VNTR it overlaps with in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_vntr_vcf_file("r", "r", "reference VNTR VCF file []",true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<int32_t> arg_window_size("w", "w", "window size for local sorting of variants [10000]", false, 10000, "integer", cmd);
            TCLAP::ValueArg<uint32_t> arg_left_window("l", "l", "left window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_right_window("t", "t", "right window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_quiet("q", "q", "do not print options and summary []", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            print = !arg_quiet.getValue();
            window_size = arg_window_size.getValue();
            left_window = arg_left_window.getValue();
            right_window = arg_right_window.getValue();
            ref_vntr_vcf_file = arg_ref_vntr_vcf_file.getValue();
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
        odw = new BCFOrderedWriter(output_vcf_file, window_size);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VNTR_OVERLAP_VARIANT,Number=.,Type=String,Description=\"Original chr:pos:ref:alt variant that overlaps with a VNTR\">\n");
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
        no_variants_vntrized = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        orom_vntrs = new OrderedBCFOverlapMatcher(ref_vntr_vcf_file, intervals);
    }

    /**
     * Updates joint allele distribution.
     */
    void update_joint_allele_dist(int32_t no_indel_alleles, int32_t no_tandem_repeat_alleles)
    {
//
//        std::cerr << "updating " << no_indel_alleles << " " << no_tandem_repeat_alleles << "\n";
//        std::cerr << joint_allele_dist.size() <<  " => ";
        if (joint_allele_dist.size()<=no_indel_alleles)
        {
            std::vector<int32_t> vec;
            for (uint32_t i = joint_allele_dist.size(); i<=no_indel_alleles; ++i)
            {
                joint_allele_dist.push_back(vec);
            }
        }

//        std::cerr << joint_allele_dist.size() << "\n";
//        std::cerr << joint_allele_dist[no_indel_alleles].size() <<  " => ";
        if (joint_allele_dist[no_indel_alleles].size()<=no_tandem_repeat_alleles)
        {
            for (uint32_t i = joint_allele_dist[no_indel_alleles].size(); i<=no_tandem_repeat_alleles; ++i)
            {
                joint_allele_dist[no_indel_alleles].push_back(0);
            }
        }
//        std::cerr << joint_allele_dist[no_indel_alleles].size() <<  "\n";

        ++joint_allele_dist[no_indel_alleles][no_tandem_repeat_alleles];
    }

    void vntrize()
    {
        uint32_t left_extended = 0;
        uint32_t left_trimmed = 0;
        uint32_t right_trimmed = 0;

        kstring_t old_alleles = {0,0,0};
        kstring_t new_alleles = {0,0,0};

        v = odw->get_bcf1_from_pool();
        bcf_hdr_t* h = odr->hdr;
        Variant variant;

        std::vector<bcf1_t *> overlap_vars;

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

            bcf_unpack(v, BCF_UN_INFO);

            std::string chrom = bcf_get_chrom(odr->hdr,v);
            int32_t rid = bcf_get_rid(v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end1(v);

            if (orom_vntrs->overlaps_with(rid, start1-left_window, end1+right_window, overlap_vars))
            {
                uint32_t no_indel_alleles = bcf_get_n_allele(v);
                uint32_t no_tandem_repeat_alleles = 0;

                for (size_t i=0; i<overlap_vars.size(); ++i)
                {
                    if (i) kputc(',', &old_alleles);
                    bcf_variant2string(orom_vntrs->odr->hdr, overlap_vars[i], &old_alleles);

                    no_tandem_repeat_alleles = bcf_get_n_allele(overlap_vars[i]) > no_tandem_repeat_alleles ? bcf_get_n_allele(overlap_vars[i]) : no_tandem_repeat_alleles;

//                    if (no_tandem_repeat_alleles == 1)
//                    {
//                        bcf_print(orom_vntrs->odr->hdr, overlap_vars[i]);
//                    }

//                    if (no_indel_alleles == 13)
//                    {
//                        bcf_print(h, v);
//                    }
                }

                update_joint_allele_dist(no_indel_alleles, no_tandem_repeat_alleles);

                bcf_variant2string(odw->hdr, v, &old_alleles);
                bcf_update_info_string(odw->hdr, v, "VNTR_OVERLAP_VARIANT", old_alleles.s);
                bcf_set_pos1(v, bcf_get_pos1(overlap_vars[0]));
                new_alleles.l = 0;
                char** allele = bcf_get_allele(overlap_vars[0]);
                int32_t n_allele = bcf_get_n_allele(overlap_vars[0]);
                for (size_t i=0; i<n_allele; ++i)
                {
                    if (i) kputc(',', &new_alleles);
                    kputs(allele[i], &new_alleles);
                }
                bcf_update_alleles_str(odw->hdr, v, new_alleles.s);

                ++no_variants_vntrized;
            }
            overlap_vars.clear();
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

        std::clog << "vntrize v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_num_op("         [l] left window           ", left_window);
        print_num_op("         [t] right window          ", right_window);
        std::clog << "         [w] sorting window size   " << window_size << "\n";
        std::clog << "         [r] reference VNTR file   " << ref_vntr_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::clog << "stats:\n";
        std::clog << "\n";
        std::clog << "       total no. variants vntrized              : " << no_variants_vntrized << "\n";
        std::clog << "       total no. variants observed              : " << no_variants << "\n";
        std::clog << "\n";


//        //determine max alleles
//        int32_t max_tr_allele_no = 0;
//        for (uint32_t i = 1; i<joint_allele_dist.size(); ++i)
//        {
//            if (joint_allele_dist[i].size()>max_tr_allele_no)
//            {
//                max_tr_allele_no = joint_allele_dist[i].size();
//            }
//        }
//        --max_tr_allele_no;
//
//        std::cerr << "Joint distribution of alleles\n";
//        for (uint32_t i = 1; i<= max_tr_allele_no; ++i) std::cerr << "\t" << i;
//        std::cerr << "\n";
//        for (uint32_t i = 1; i<joint_allele_dist.size(); ++i)
//        {
//            std::cerr << i ;
//            for (uint32_t j = 1; j<joint_allele_dist[i].size(); ++j)
//            {
//                std::cerr << "\t";
//                std::cerr << joint_allele_dist[i][j];
//            }
//            std::cerr << "\n";
//        }
    };

    ~Igor() {};

    private:
};

}

bool vntrize(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.vntrize();
    igor.print_stats();

    return igor.print;
};
