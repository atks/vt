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

#include "partition.h"

namespace
{

class OverlapStats
{
    public:

    uint32_t a,ab,b,a_ins,ab_ins,b_ins,a_del,ab_del,b_del,a_ts,ab_ts,b_ts,a_tv,ab_tv,b_tv;

    OverlapStats()
    {
        a = 0;
        ab = 0;
        b = 0;

        a_ts = 0;
        a_tv = 0;
        ab_ts = 0;
        ab_tv = 0;
        b_ts = 0;
        b_tv = 0;

        a_ins = 0;
        a_del = 0;
        ab_ins = 0;
        ab_del = 0;
        b_ins = 0;
        b_del = 0;
    };
};

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string filters;
    std::vector<std::string> input_vcf_files;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;

    /////////
    //stats//
    /////////
    std::vector<OverlapStats> stats;

    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "partition variants. check the overlap of variants between 2 data sets.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_filters("f", "filters", "Filter (e.g. AF>0.3) ", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf><in2.vcf>", "2 input VCF files for comparison", true, "files", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
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
        sr = new BCFSyncedReader(input_vcf_files, intervals, false);

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");

        ////////////////////////
        //stats initialization//
        ////////////////////////
        OverlapStats s;
        stats.resize(32, s);

    }

    void partition()
    {
        //for combining the alleles
        std::vector<bcfptr*> current_recs;

        Variant variant;
        std::vector<int32_t> presence(2);

        while(sr->read_next_position(current_recs))
        {
            bcf1_t *v = current_recs[0]->v;
            bcf_hdr_t *h = current_recs[0]->h;
            std::string chrom = bcf_get_chrom(h,v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);
            int32_t vtype = vm->classify_variant(h, v, variant);

            //check existence
            for (uint32_t i=0; i<current_recs.size(); ++i)
            {
                ++presence[current_recs[i]->file_index];
            }

            int32_t ins = 0;
            int32_t del = 0;
            int32_t ts = 0;
            int32_t tv = 0;

            update_overlap_stats(stats, presence, 0, ts, tv, ins, del);

            if (bcf_get_n_allele(v)==2)
            {
                if (vtype == VT_SNP || vtype == VT_MNP  )
                {
                    ts = variant.alleles[0].ts;
                    tv = variant.alleles[0].tv;
                }
                else if (vtype == VT_INDEL)
                {
                    ins = variant.alleles[0].ins;
                    del = 1-ins;
                }

                update_overlap_stats(stats, presence, vtype, ts, tv, ins, del);

                if (vtype & VT_CLUMPED)
                {
                    
                }
                else if (vtype==VT_SNP)
                {
                    
                }
                else if (vtype==VT_MNP)
                {
                   
                }
                else if (vtype==VT_INDEL) //strictly simple indels
                {
                    
                }
                else if (vtype==(VT_SNP|VT_INDEL))
                {
                    
                }
                else if (vtype==(VT_MNP|VT_INDEL))
                {
                   
                }
                else if (vtype==(VT_SNP|VT_MNP|VT_INDEL))
                {
                  
                }
                else if (vtype==(VT_SNP|VT_MNP))
                {
                  
                }
                else if (vtype==VT_REF) //MNPs that are not real MNPs
                {
                  
                }
                else
                {
                    
                }

            }

            presence[0] = 0;
            presence[1] = 0;
        }
    };

    void update_overlap_stats(std::vector<OverlapStats>& stats, std::vector<int32_t>& presence, int32_t vtype, int32_t ts,  int32_t tv, int32_t ins, int32_t del)
    {
        //update overlap stats
        if (presence[0] && !presence[1])
        {
            ++stats[vtype].a;

            stats[vtype].a_ts += ts;
            stats[vtype].a_tv += tv;
            stats[vtype].a_ins += ins;
            stats[vtype].a_del += del;
        }
        else if (presence[0] && presence[1])
        {
            ++stats[vtype].ab;

            stats[vtype].ab_ts += ts;
            stats[vtype].ab_tv += tv;
            stats[vtype].ab_ins += ins;
            stats[vtype].ab_del += del;
        }
        else if (!presence[0] && presence[1])
        {
            ++stats[vtype].b;

            stats[vtype].b_ts += ts;
            stats[vtype].b_tv += tv;
            stats[vtype].b_ins += ins;
            stats[vtype].b_del += del;
        }        
    }

    void print_options()
    {
        std::clog << "partition v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF file a   " << input_vcf_files[0] << "\n";
        std::clog << "             input VCF file b   " << input_vcf_files[1] << "\n";
        print_int_op("         [i] intervals          ", intervals);
        std::clog << "\n";
   }

    void print_stats()
    {
        fprintf(stderr, "    A:  %10d variants\n", stats[0].a+stats[0].ab);
        fprintf(stderr, "    B:  %10d variants\n", stats[0].ab+stats[0].b);
        fprintf(stderr, "\n");

        int32_t types[5] = {0, VT_SNP, VT_MNP, VT_INDEL, VT_CLUMPED};
        kstring_t s;

        for (int32_t j=0; j<5; ++j)
        {
            int32_t i = types[j];
            vm->vtype2string(i, &s);
            if (types[j] == 0)
            {
                fprintf(stderr, "    ALL\n");
                fprintf(stderr, "    A-B %10d\n", stats[i].a);
                fprintf(stderr, "    A&B %10d\n", stats[i].ab);
                fprintf(stderr, "    B-A %10d\n", stats[i].b);
                fprintf(stderr, "    of A     %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
                fprintf(stderr, "    of B     %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
                fprintf(stderr, "\n");
            }
            else if (types[j] == VT_SNP || types[j] == VT_MNP)
            {
                fprintf(stderr, "    %s\n", s.s);
                fprintf(stderr, "    A-B %10d [%.2f]\n", stats[i].a,  (float)stats[i].a_ts/(stats[i].a_tv));
                fprintf(stderr, "    A&B %10d [%.2f]\n", stats[i].ab, (float)stats[i].ab_ts/stats[i].ab_tv);
                fprintf(stderr, "    B-A %10d [%.2f]\n", stats[i].b,  (float)stats[i].b_ts/(stats[i].b_tv));
                fprintf(stderr, "    of A     %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
                fprintf(stderr, "    of B     %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
                fprintf(stderr, "\n");
            }
            else
            {
                fprintf(stderr, "    %s\n", s.s);
                fprintf(stderr, "    A-B %10d [%.2f]\n", stats[i].a,  (float)stats[i].a_ins/(stats[i].a_del));
                fprintf(stderr, "    A&B %10d [%.2f]\n", stats[i].ab, (float)stats[i].ab_ins/stats[i].ab_del);
                fprintf(stderr, "    B-A %10d [%.2f]\n", stats[i].b,  (float)stats[i].b_ins/(stats[i].b_del));
                fprintf(stderr, "    of A     %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
                fprintf(stderr, "    of B     %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
                fprintf(stderr, "\n");
            }
        }
    };

    ~Igor()
    {
    };

    private:

};

}

void partition(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.partition();
    igor.print_stats();
}
