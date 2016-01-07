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

#include "profile_vntrs.h"

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
    OrderedBCFOverlapMatcher *obom;

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
        //reference data set//
        //////////////////////
//# This file contains information on how to process reference data sets.
//# dataset - name of data set, this label will be printed.
//# type    - True Positives (TP) and False Positives (FP).
//#           overlap percentages labeled as (Precision, Sensitivity) and (False Discovery Rate, Type I Error) respectively.
//#         - annotation.
//#           file is used for GENCODE annotation of frame shift and non frame shift Indels.
//# filter  - filter applied to variants for this particular data set.
//# path    - path of indexed BCF file.
//#dataset     type            filter                       path
//trf_lobstr   TP              VTYPE==VNTR                  /net/fantasia/home/atks/ref/vt/grch37/trf.lobstr.sites.bcf

        dataset_labels.push_back("data");
        dataset_types.push_back(REF);
        dataset_fexps.push_back(fexp);
        dataset_info_site_tags.push_back("");
        dataset_info_gt_tags.push_back("");

        htsFile *hts = hts_open(ref_data_sets_list.c_str(), "r");
        if (!hts)
        {
            fprintf(stderr, "[E:%s:%d %s] Reference file cannot be opened %s\n", __FILE__, __LINE__, __FUNCTION__, ref_data_sets_list.c_str());
            exit(1);
        }

        kstring_t s = {0,0,0};
        std::vector<std::string> vec;
        while (hts_getline(hts, '\n', &s)>=0)
        {
            if (s.s[0] == '#')
                continue;

            std::string line(s.s);
            split(vec, " ", line);

            if (vec[1] == "Truth" || vec[1] == "False" || vec[1] == "BroadKB")
            {
                dataset_labels.push_back(vec[0]);
                if (vec[1]=="Truth")
                {
                    dataset_types.push_back(TRUE);
                }
                else if (vec[1]=="False")
                {
                    dataset_types.push_back(FALSE);
                }
                else if (vec[1]=="BroadKB")
                {
                    dataset_types.push_back(BROADKB);
                }
                dataset_fexps.push_back(vec[2]);
                input_vcf_files.push_back(vec[3]);
            }
            else if (vec[1] == "cds_annotation")
            {
                cds_bed_file = vec[3];
            }
            else
            {
                fprintf(stderr, "[E:%s:%d %s] Reference data set type %s not recognized\n", __FILE__, __LINE__, __FUNCTION__, vec[1].c_str());
                exit(1);
            }
        }
        hts_close(hts);
        if (s.m) free(s.s);
            
        //////////////////////
        //i/o initialization//
        //////////////////////
        sr = new BCFSyncedReader(input_vcf_files, intervals, SYNC_BY_VAR);
        obom = new OrderedBCFOverlapMatcher(vcf_file);

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
        std::clog << "profile_vntrs v" << version << "\n";
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

void profile_vntrs(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.rpartition();
    igor.print_stats();
}
