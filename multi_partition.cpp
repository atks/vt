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

#include "multi_partition.h"

namespace
{

class OverlapStats
{
    public:

    uint32_t no_datasets;
    std::vector<uint32_t> total;
    std::vector<uint32_t> no;
    std::vector<uint32_t> ts;
    std::vector<uint32_t> tv;
    std::vector<uint32_t> ins;
    std::vector<uint32_t> del;

    /**
     * Constructor.
     */
    OverlapStats(uint32_t no_datasets=2)
    {
        resize(no_datasets);
    };

    /**
     * Resize for no_datasets.
     */
    void resize(uint32_t no_datasets)
    {
        this->no_datasets = no_datasets;
        resize();
    };

    /**
     * Resize for no_datasets.
     */
    void resize()
    {
        no.resize(no_datasets, 0);
        no.resize((1 << no_datasets), 0);
        ts.resize((1 << no_datasets), 0);
        tv.resize((1 << no_datasets), 0);
        ins.resize((1 << no_datasets), 0);
        del.resize((1 << no_datasets), 0);
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
            std::string desc = "partition variants from any number of data sets.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf><in2.vcf>...", "multiple input VCF files for comparison", true, "files", cmd);

            cmd.parse(argc, argv);

            fexp = arg_fexp.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            input_vcf_files = arg_input_vcf_files.getValue();

            if (input_vcf_files.size()<2)
            {
                std::cerr << "error: require 2 or more VCF files\n";
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
        stats.resize(input_vcf_files.size());
    }

    void partition()
    {
        //for combining the alleles
        std::vector<bcfptr*> crecs;
        Variant variant;

        std::vector<BCFOrderedWriter *> out;

        while(sr->read_next_position(crecs))
        {
            int32_t vtype = vm->classify_variant(crecs[0]->h, crecs[0]->v, variant);

            uint32_t presence = 0;
            for (uint32_t i=0; i<crecs.size(); ++i)
            {
                if (filter_exists && !filter.apply(crecs[i]->h,crecs[i]->v,&variant))
                {
                    continue;
                }

                presence |= 1 << crecs[i]->file_index;

            }

            uint32_t ts = 0;
            uint32_t tv = 0;
            uint32_t ins = 0;
            uint32_t del = 0;
            std::vector<Allele>& alleles = variant.alleles;
            for (uint32_t i=0; i<alleles.size(); ++i)
            {
                ts += alleles[i].ts;
                tv += alleles[i].tv;
                ins += alleles[i].dlen ? variant.alleles[i].ins : 0;
                del += alleles[i].dlen ? 1-variant.alleles[i].ins : 0;
            }

            ++stats.no[presence];
            stats.ts[presence] += ts;
            stats.tv[presence] += tv;
            stats.ins[presence] += ins;
            stats.del[presence] += del;
        }
    };

    void print_options()
    {
        std::clog << "multi_partition v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF file a   " << input_vcf_files[0] << "\n";
        char c = 'b';
        for (uint32_t i=1; i<input_vcf_files.size(); ++i)
        {
            std::clog << "             input VCF file " << c << "   " << input_vcf_files[i] << "\n";
            ++c;
        }
        print_str_op("         [f] filter             ", fexp);
        print_int_op("         [i] intervals          ", intervals);
        std::clog << "\n";
   }

    void print_stats()
    {
        //print size of variants
        char c = 'A';
        for (uint32_t i=1; i<=stats.no_datasets; ++i)
        {
            uint32_t k = 1<<(i-1);
            uint32_t no = 0;
            for (uint32_t j=1; j<=stats.no.size(); ++j)
            {
                if (j&k) no+=stats.no[j];
            }

            fprintf(stderr, "    %c:  %10d variants\n", c, no);
            ++c;
        }
        fprintf(stderr, "\n");

        //print partitions
        fprintf(stderr, "    %s         no  [ts/tv] [ins/del]\n", std::string(stats.no_datasets,' ').c_str());
        uint32_t unique_total = 0;
        for (uint32_t i=1; i<stats.no.size(); ++i)
        {
            std::string partition;
            for (uint32_t j=1; j<=stats.no_datasets; ++j)
            {
                if (i&(1<<(j-1)))
                {
                    partition.append(1, 'A'+(j-1));
                }
                else
                {
                    partition.append(1, '-');
                }
            }

            unique_total += stats.no[i];
            fprintf(stderr, "    %s %10d  [%.2f]  [%.2f]\n", partition.c_str(), stats.no[i], (float)stats.ts[i]/stats.tv[i], (float)stats.ins[i]/stats.del[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "    Unique variants     : %10d\n", unique_total);
        fprintf(stderr, "    Overall concordance : %10.2f%% (#intersection/#union)\n", (float) stats.no[stats.no.size()-1]/unique_total*100);
        fprintf(stderr, "\n");
    };

    ~Igor()
    {
    };

    private:

};
}

void multi_partition(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.partition();
    igor.print_stats();
}
