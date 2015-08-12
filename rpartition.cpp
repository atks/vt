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

#include "rpartition.h"

namespace
{

class OverlapStats
{
    public:

    uint32_t a,ab,b;
        
    OverlapStats()
    {
        a = 0;
        ab = 0;
        b = 0;
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
    float rfrac;

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
            std::string desc = "reciprocal partition of VNTRs between two data sets.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter", false, "", "str", cmd);
            TCLAP::ValueArg<float> arg_rfrac("r", "r", "reciprocal overlap ", false, 0.9, "float", cmd);
            TCLAP::SwitchArg arg_write_partition("w", "w", "write partitioned variants to file", cmd, false);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf><in2.vcf>", "2 input VCF files for comparison", true, "files", cmd);

            cmd.parse(argc, argv);

            fexp = arg_fexp.getValue();
            rfrac = arg_rfrac.getValue();
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

    void rpartition()
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
                if (vtype != VT_VNTR)
                {
                    continue;
                }
                
                if (filter_exists && !filter.apply(crecs[i]->h,crecs[i]->v,&variant))
                {
                    continue;
                }

                ++presence[crecs[i]->file_index];
            }

            update_overlap_stats(presence);

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

    void update_overlap_stats(std::vector<int32_t>& presence)
    {
        //update overlap stats
        if (presence[0] && !presence[1])
        {
            ++stats.a;
        }
        else if (presence[0] && presence[1])
        {
            ++stats.ab;
        }
        else if (!presence[0] && presence[1])
        {
            ++stats.b;
        }
    }

    void print_options()
    {
        std::clog << "rpartition v" << version << "\n";
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
        fprintf(stderr, "    A:  %10d VNTRs\n", stats.a+stats.ab);
        fprintf(stderr, "    B:  %10d VNTRs\n", stats.ab+stats.b);
        fprintf(stderr, "\n");
        fprintf(stderr, "                   ts/tv  ins/del\n");
        fprintf(stderr, "    A-B %10d\n", stats.a);
        fprintf(stderr, "    A&B %10d\n", stats.ab);
        fprintf(stderr, "    B-A %10d\n", stats.b);
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

void rpartition(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.rpartition();
    igor.print_stats();
}
