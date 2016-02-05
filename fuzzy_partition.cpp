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

#include "fuzzy_partition.h"

namespace
{
    
class FuzzyOverlapStats
{
    public:

    uint32_t a, ap, ab, bp, b;
    
    FuzzyOverlapStats()
    {
        a = 0;
        ap= 0;
        ab = 0;
        bp = 0;
        b = 0;
    };
};

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_vcf_files;
    std::vector<GenomeInterval> intervals;
    uint32_t left_window;
    uint32_t right_window;
    bool print;
    bool write_partition;

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
    
    /////////
    //stats//
    /////////
    FuzzyOverlapStats stats;

    //for allele counts [no indel alleles][no tandem repeat alleles]
    std::vector<std::vector<int32_t> > joint_allele_dist;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    OrderedBCFOverlapMatcher *orom_vcf;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "fuzzy partition 2 VCF files";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<uint32_t> arg_left_window("l", "l", "left window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_right_window("r", "r", "right window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_write_partition("w", "w", "write partitioned variants to file", cmd, false);
            TCLAP::SwitchArg arg_quiet("q", "q", "do not print options and summary []", cmd, false);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf><in2.vcf>", "2 input VCF files for comparison", true, "files", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print = !arg_quiet.getValue();
            left_window = arg_left_window.getValue();
            right_window = arg_right_window.getValue();
            fexp = arg_fexp.getValue();
            write_partition = arg_write_partition.getValue();
            input_vcf_files = arg_input_vcf_files.getValue();

            if (input_vcf_files.size()!=2)
            {
                std::cerr << "error: require 2 VCF files\n";
                exit(1);
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
        odr = new BCFOrderedReader(input_vcf_files[0], intervals);
        bcf_hdr_append(odr->hdr, "##INFO=<ID=OVERLAP,Number=.,Type=Integer,Description=\"Overlap Count\">\n");
      
        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;
                
        ////////////////////////
        //tools initialization//
        ////////////////////////
        orom_vcf = new OrderedBCFOverlapMatcher(input_vcf_files[1], intervals);
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
    
    void fuzzy_partition()
    {
        bcf1_t *v = bcf_init1();
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
            int32_t no_tr_alleles = bcf_get_n_allele(v);
            int32_t beg1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end1(v);

            if (orom_vcf->overlaps_with(chrom, beg1-left_window, end1+right_window, overlap_vars))
            {
                int32_t no_overlaps = overlap_vars.size();
                update_joint_allele_dist(no_tr_alleles, no_overlaps);
//                std::cerr << no_tr_alleles << " " << no_overlaps << "\n";
            }
            
//            std::cerr << "hello\n";
            
        }

        odr->close();
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "fuzzy partition v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF file a   " << input_vcf_files[0] << "\n";
        std::clog << "             input VCF file b   " << input_vcf_files[1] << "\n";
        print_num_op("         [l] left window           ", left_window);
        print_num_op("         [r] right window          ", right_window);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::clog << "stats:\n";
        fprintf(stderr, "    A-B  %10d \n", stats.a);
        fprintf(stderr, "    A-B~ %10d \n", stats.ap);
        fprintf(stderr, "    A&B  %10d \n", stats.ab);
        fprintf(stderr, "    B-A~ %10d \n", stats.bp);
        fprintf(stderr, "    B-A  %10d \n", stats.b);
        fprintf(stderr, "\n");
        
        
        //determine max alleles
        int32_t max_tr_allele_no = 0;
        for (uint32_t i = 1; i<joint_allele_dist.size(); ++i)
        {   
            if (joint_allele_dist[i].size()>max_tr_allele_no)
            {
                max_tr_allele_no = joint_allele_dist[i].size();
            }
        }
        --max_tr_allele_no;
        
        std::cerr << "Joint distribution of alleles\n";    
        for (uint32_t i = 1; i<= max_tr_allele_no; ++i) std::cerr << "\t" << i;
        std::cerr << "\n";
        for (uint32_t i = 1; i<joint_allele_dist.size(); ++i)
        {
            std::cerr << i ;
            for (uint32_t j = 1; j<joint_allele_dist[i].size(); ++j)
            {
                std::cerr << "\t";
                std::cerr << joint_allele_dist[i][j];
            }
            std::cerr << "\n";
        }
    };

    ~Igor() {};

    private:
};

}

void fuzzy_partition(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.fuzzy_partition();
    igor.print_stats();
};
