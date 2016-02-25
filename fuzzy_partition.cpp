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

    uint32_t a, ap, ab1, ab2, bp, b;

    FuzzyOverlapStats()
    {
        a = 0;
        ap = 0;
        ab1 = 0;
        ab2 = 0;
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
    std::string a_vcf_file;
    std::string b_vcf_file;

    //////////
    //filter//
    //////////
    std::vector<std::string> fexps;
    Filter filter;
    bool filter_exists;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    OrderedBCFOverlapMatcher *obom;
    BCFOrderedWriter *a_odw;
    BCFOrderedWriter *b_odw;

    /////////
    //stats//
    /////////
    FuzzyOverlapStats stats;
    int32_t no_variants;

    //for allele counts [no indel alleles][no tandem repeat alleles]
    std::vector<std::vector<int32_t> > joint_allele_dist;

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
            std::string desc = "fuzzy partition 2 VCF files";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<uint32_t> arg_left_window("l", "l", "left window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_right_window("r", "r", "right window size for overlap []", false, 0, "int", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_a_vcf_file("a", "a", "output file name for first VCF file [a.bcf]", false, "a.bcf", "str", cmd);
            TCLAP::ValueArg<std::string> arg_b_vcf_file("b", "b", "output file name for second VCF file [b.bcf]", false, "b.bcf", "str", cmd);
            TCLAP::SwitchArg arg_write_partition("w", "w", "write partitioned variants to file", cmd, false);
            TCLAP::SwitchArg arg_quiet("q", "q", "do not print options and summary []", cmd, false);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf><in2.vcf>", "2 input VCF files for comparison", true, "files", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print = !arg_quiet.getValue();
            left_window = arg_left_window.getValue();
            right_window = arg_right_window.getValue();
            parse_filters(fexps, arg_fexp.getValue(), 2, false);
            write_partition = arg_write_partition.getValue();
            a_vcf_file = arg_a_vcf_file.getValue();
            b_vcf_file = arg_b_vcf_file.getValue();
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
        obom = new OrderedBCFOverlapMatcher(input_vcf_files[1], intervals, fexps[1]);

        bcf_hdr_append_info_with_backup_naming(odr->hdr, "EXACT_OVERLAPS", "1", "Integer", "Number of exact overlapping variants with this variant.", true);
        bcf_hdr_append_info_with_backup_naming(odr->hdr, "FUZZY_OVERLAPS", "1", "Integer", "Number of fuzzy overlapping variants with this variant.", true);
        bcf_hdr_sync(odr->hdr);

        a_odw = NULL;
        b_odw = NULL;
        if (write_partition)
        {
            a_odw = new BCFOrderedWriter(a_vcf_file);
            a_odw->link_hdr(odr->hdr);
            a_odw->write_hdr();
            b_odw = new BCFOrderedWriter(b_vcf_file);
            b_odw->link_hdr(obom->odr->hdr);
            b_odw->write_hdr();
        }

        ///////
        //stats
        ///////
        no_variants = 0;

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexps[0].c_str());
        filter_exists = fexps[0]=="" ? false : true;

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    /**
     * Increments the EXACT_OVERLAPS count of a variant record.
     */
    void increment_exact_overlap(bcf_hdr_t* h, bcf1_t* v, int32_t k)
    {
        int32_t n = 0;
        int32_t *count = NULL;
        bcf_unpack(v, BCF_UN_INFO);
        if (bcf_get_info_int32(h, v, "EXACT_OVERLAPS", &count, &n)>0)
        {
            count[0] = std::min(127, count[0]+k);
            bcf_update_info_int32(h, v, "EXACT_OVERLAPS", count, n);
            free(count);
        }
        else
        {
            int32_t c = 1;
            bcf_update_info_int32(h, v, "EXACT_OVERLAPS", &c, 1);
        }
    }

    /**
     * Increments the FUZZY_OVERLAPS count of a variant record.
     */
    void increment_fuzzy_overlap(bcf_hdr_t* h, bcf1_t* v, int32_t k)
    {
        int32_t n = 0;
        int32_t *count = NULL;
        bcf_unpack(v, BCF_UN_INFO);
        if (bcf_get_info_int32(h, v, "FUZZY_OVERLAPS", &count, &n)>0)
        {
            count[0] = std::min(127, count[0]+k);
            bcf_update_info_int32(h, v, "FUZZY_OVERLAPS", count, n);
            free(count);
        }
        else
        {
            int32_t c = 1;
            bcf_update_info_int32(h, v, "FUZZY_OVERLAPS", &c, 1);
        }
    }


    /**
     * Updates joint allele distribution.
     */
    void update_joint_allele_dist(int32_t no_indel_alleles, int32_t no_tandem_repeat_alleles)
    {
        if (joint_allele_dist.size()<=no_indel_alleles)
        {
            std::vector<int32_t> vec;
            for (uint32_t i = joint_allele_dist.size(); i<=no_indel_alleles; ++i)
            {
                joint_allele_dist.push_back(vec);
            }
        }

        if (joint_allele_dist[no_indel_alleles].size()<=no_tandem_repeat_alleles)
        {
            for (uint32_t i = joint_allele_dist[no_indel_alleles].size(); i<=no_tandem_repeat_alleles; ++i)
            {
                joint_allele_dist[no_indel_alleles].push_back(0);
            }
        }

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
            int32_t rid = bcf_get_rid(v);
            int32_t no_tr_alleles = bcf_get_n_allele(v);
            int32_t beg1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end1(v);

            if (obom->overlaps_with(rid, beg1-left_window, end1+right_window, overlap_vars, b_odw))
            {
                int32_t no_overlaps = overlap_vars.size();
                update_joint_allele_dist(no_tr_alleles, no_overlaps);
                
                bool exact = false;
                int32_t partial_overlap = 0;
                int32_t exact_overlap = 0;
                for (uint32_t j=0; j<overlap_vars.size(); ++j)
                {
                    //check for exactness
                    if (obom->is_exact_match(rid, beg1, end1, overlap_vars[j]))
                    {
                        exact = true;
                        ++exact_overlap;
                    }
                    else
                    {
                        ++partial_overlap;
                    }
                }
                
                increment_exact_overlap(h, v, exact_overlap);
                increment_fuzzy_overlap(h, v, partial_overlap);
                
                if (exact)
                {
                    ++stats.ab1;
                }
                else
                {
                    ++stats.ap;
                }
            }
            else
            {
                ++stats.a;
            }
            
            if (write_partition)
            {
                a_odw->write(v);
            }

            
            ++no_variants;
        }

        obom->flush(b_odw);

        stats.ab2 += obom->get_no_exact_overlaps();
        stats.bp += obom->get_no_fuzzy_overlaps();
        stats.b += obom->get_no_nonoverlaps();

        odr->close();
        obom->close();
        
        if (write_partition)
        {    
            a_odw->close();
            b_odw->close();
        }
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
        if (write_partition)
        {
            std::clog << "         [w] write_partition    true (partitions will be written to a.bcf and b.bcf\n";
            std::clog << "         [a] output VCF file a  " << a_vcf_file << "\n";
            std::clog << "         [b] output VCF file b  " << b_vcf_file << "\n";    
        }    
        else
        {
            std::clog << "         [w] write_partition    false\n";

        }
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::clog << "stats:\n";
        fprintf(stderr, "\n"); 
        fprintf(stderr, "    A :  %d\n",no_variants);  
        fprintf(stderr, "    B :  %d\n",obom->no_variants);    
        fprintf(stderr, "\n");    
        fprintf(stderr, "    A-B  %10d \n", stats.a);
        fprintf(stderr, "    A-B~ %10d \n", stats.ap);
        fprintf(stderr, "    A&B1 %10d \n", stats.ab1);
        fprintf(stderr, "    A&B2 %10d \n", stats.ab2);
        fprintf(stderr, "    B-A~ %10d \n", stats.bp);
        fprintf(stderr, "    B-A  %10d \n", stats.b);
        fprintf(stderr, "\n");

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

void fuzzy_partition(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.fuzzy_partition();
    igor.print_stats();
};
