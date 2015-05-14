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

    uint32_t a,ab,b,a_ins,ab_ins,b_ins,a_del,ab_del,b_del;
    uint32_t a_ts,ab_ts,b_ts,a_tv,ab_tv,b_tv;
        
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
    bool write_partition;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    OverlapStats stats;

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
            std::string desc = "partition variants from two data sets.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter", false, "", "str", cmd);
            TCLAP::SwitchArg arg_write_partition("w", "w", "write partitioned variants to file", cmd, false);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf><in2.vcf>", "2 input VCF files for comparison", true, "files", cmd);

            cmd.parse(argc, argv);

            fexp = arg_fexp.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
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
        sr = new BCFSyncedReader(input_vcf_files, intervals, SYNC_BY_VAR);

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
    }

    void partition()
    {
        //for combining the alleles
        std::vector<bcfptr*> crecs;
        Variant variant;
        std::vector<int32_t> presence(2, 0);

        BCFOrderedWriter *a;
        BCFOrderedWriter *ab[2];
        BCFOrderedWriter *b;
        if (write_partition)
        {
            a = new BCFOrderedWriter("a-b.bcf");
            a->link_hdr(sr->hdrs[0]);
            a->write_hdr();
            ab[0] = new BCFOrderedWriter("a&b1.bcf");
            ab[0]->link_hdr(sr->hdrs[0]);
            ab[0]->write_hdr();
            ab[1] = new BCFOrderedWriter("a&b2.bcf");
            ab[1]->link_hdr(sr->hdrs[1]);
            ab[1]->write_hdr();
            b = new BCFOrderedWriter("b-a.bcf");
            b->link_hdr(sr->hdrs[1]);
            b->write_hdr();
        }

        while(sr->read_next_position(crecs))
        {
            int32_t vtype = vm->classify_variant(crecs[0]->h, crecs[0]->v, variant);

            //check existence
            for (int32_t i=0; i<crecs.size(); ++i)
            {
                if (filter_exists && !filter.apply(crecs[i]->h,crecs[i]->v,&variant))
                {
                    continue;
                }

                ++presence[crecs[i]->file_index];
            }

            int32_t ts = 0;
            int32_t tv = 0;
            int32_t ins = 0;
            int32_t del = 0;
            std::vector<Allele>& alleles = variant.alleles;
            for (int32_t i=0; i<alleles.size(); ++i)
            {
                ts += alleles[i].ts;
                tv += alleles[i].tv;
                ins += alleles[i].dlen ? variant.alleles[i].ins : 0;
                del += alleles[i].dlen ? 1-variant.alleles[i].ins : 0;
            }

            update_overlap_stats(presence, ts, tv, ins, del);

            if (write_partition)
            {
                if (presence[0])
                {
                    if (presence[1])
                    {
                        ab[crecs[0]->file_index]->write(crecs[0]->v);
                        ab[crecs[1]->file_index]->write(crecs[1]->v);
                    }
                    else
                    {
                        a->write(crecs[0]->v);
                    }
                }
                else if (presence[1])
                {
                    b->write(crecs[0]->v);
                }
            }

            presence[0] = 0;
            presence[1] = 0;
        }

        if (write_partition)
        {
            a->close();
            ab[0]->close();
            ab[1]->close();
            b->close();
        }
    };

    void update_overlap_stats(std::vector<int32_t>& presence, int32_t ts,  int32_t tv, int32_t ins, int32_t del)
    {
        //update overlap stats
        if (presence[0] && !presence[1])
        {
            ++stats.a;

            stats.a_ts += ts;
            stats.a_tv += tv;
            stats.a_ins += ins;
            stats.a_del += del;
        }
        else if (presence[0] && presence[1])
        {
            ++stats.ab;

            stats.ab_ts += ts;
            stats.ab_tv += tv;
            stats.ab_ins += ins;
            stats.ab_del += del;
        }
        else if (!presence[0] && presence[1])
        {
            ++stats.b;

            stats.b_ts += ts;
            stats.b_tv += tv;
            stats.b_ins += ins;
            stats.b_del += del;
        }
    }

    void print_options()
    {
        std::clog << "partition v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF file a   " << input_vcf_files[0] << "\n";
        std::clog << "             input VCF file b   " << input_vcf_files[1] << "\n";
        print_str_op("         [f] filter             ", fexp);
        print_boo_op("         [w] write_partition    ", write_partition);
        print_int_op("         [i] intervals          ", intervals);
        std::clog << "\n";
   }

    void print_stats()
    {
        fprintf(stderr, "    A:  %10d variants\n", stats.a+stats.ab);
        fprintf(stderr, "    B:  %10d variants\n", stats.ab+stats.b);
        fprintf(stderr, "\n");
        fprintf(stderr, "                   ts/tv  ins/del\n");
        fprintf(stderr, "    A-B %10d [%.2f] [%.2f]\n", stats.a,  (float)stats.a_ts/(stats.a_tv), (float)stats.a_ins/(stats.a_del));
        fprintf(stderr, "    A&B %10d [%.2f] [%.2f]\n", stats.ab, (float)stats.ab_ts/stats.ab_tv, (float)stats.ab_ins/(stats.ab_del));
        fprintf(stderr, "    B-A %10d [%.2f] [%.2f]\n", stats.b,  (float)stats.b_ts/(stats.b_tv), (float)stats.b_ins/(stats.b_del));
        fprintf(stderr, "    of A     %4.1f%%\n", 100*(float)stats.ab/(stats.a+stats.ab));
        fprintf(stderr, "    of B     %4.1f%%\n", 100*(float)stats.ab/(stats.b+stats.ab));
        fprintf(stderr, "\n");
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
