/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#include "xcmp.h"

namespace
{

class OverlapStats
{
    public:

    int32_t a,ab,b,e,a_ts,ab_ts,b_ts,e_ts,a_tv,ab_tv,b_tv,e_tv,a_ins,ab_ins,b_ins,e_ins,a_del,ab_del,b_del,e_del;

    OverlapStats()
    {
        a = 0;
        ab = 0;
        b = 0;
        e = 0;

        a_ts = 0;
        a_tv = 0;
        ab_ts = 0;
        ab_tv = 0;
        b_ts = 0;
        b_tv = 0;
        e_ts = 0;
        e_tv = 0;

        a_ins = 0;
        a_del = 0;
        ab_ins = 0;
        ab_del = 0;
        b_ins = 0;
        b_del = 0;
        e_ins = 0;
        e_del = 0;
    };
};

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_vcf_files;
    std::string output_tabulate_dir;
    std::string output_pdf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    //////////////////////////////////////////////
    //reference file info : to store in an object?
    //////////////////////////////////////////////
    int32_t no_files;
    std::vector<std::string> dataset_labels;
    std::vector<std::string> dataset_fexps;
    std::string cds_bed_file;
    std::string cplx_bed_file;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    bcf1_t *v;

    //////////
    //filter//
    //////////
    std::string fexp;
    std::vector<Filter> filters;
    std::vector<bool> filter_exists;
    int32_t no_filters;

    /////////
    //stats//
    /////////
    std::vector<std::vector<OverlapStats> > stats;

    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;
    OrderedRegionOverlapMatcher *orom_lcplx;
    OrderedRegionOverlapMatcher *orom_gencode_cds;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Cross compares VCF files";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_tabulate_dir("x", "x", "output latex directory []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output pdf file []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_labels("l", "l", "Comma delimited labels for the files", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_filters(dataset_fexps, arg_fexp.getValue());
            output_tabulate_dir = arg_output_tabulate_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
            split(dataset_labels, ",", arg_labels.getValue());
            split(dataset_fexps, ",", arg_fexp.getValue());
            input_vcf_files = arg_input_vcf_files.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        no_files = input_vcf_files.size();
        
        if (dataset_labels.size()==0)
        {
            kstring_t s = {0,0,0};

            for (size_t i=0; i<no_files; ++i)
            {
                s.l = 0;
                kputs("data", &s);
                kputw(i, &s);
                dataset_labels.push_back(std::string(s.s));
            }

            if (s.m) free(s.s);
        }

        /////////////////////////
        //filter initialization//
        /////////////////////////
        if (dataset_fexps.size()==0)
        {
            for (size_t i=0; i<no_files; ++i)
            {
                filters.push_back(Filter(""));
                filter_exists.push_back(false);
            }
        }
        else if (dataset_fexps.size()==1)
        {
            for (size_t i=0; i<no_files; ++i)
            {
                filters.push_back(Filter(dataset_fexps[0]));
                filter_exists.push_back(dataset_fexps[0]!="");
            }
        }
        else if (dataset_fexps.size()==no_files)
        {
            for (size_t i=0; i<no_files; ++i)
            {
                filters.push_back(Filter(dataset_fexps[i]));
                filter_exists.push_back(dataset_fexps[i]!="");
            }
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Number of filter expressions should either be 1 or the same as the number of input vcf files\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }

        std::vector<OverlapStats> rstats(no_files);
        for (size_t i=0; i<no_files; ++i)
        {
            stats.push_back(rstats);
        }

        //////////////////////
        //i/o initialization//
        //////////////////////
        sr = new BCFSyncedReader(input_vcf_files, intervals, SYNC_BY_VAR);

        ///////////////////////
        //tool initialization//
        ///////////////////////

        ////////////////////////
        //stats initialization//
        ////////////////////////

    }

    void xcmp()
    {
        //for combining the alleles
        std::vector<bcfptr*> current_recs;
        std::vector<Interval*> overlaps;
        Variant variant;
        std::vector<int32_t> presence(no_files, 0);
      
        while(sr->read_next_position(current_recs))
        {
            //check first variant
            bcf1_t *v = current_recs[0]->v;
            bcf_hdr_t *h = current_recs[0]->h;
            int32_t vtype = vm->classify_variant(h, v, variant);

            std::string chrom = bcf_get_chrom(h,v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);

            for (size_t i=0; i<no_files; ++i)
                presence[i]=0;

            //check existence
            int32_t no_filtered = 0;
            for (size_t i=0; i<current_recs.size(); ++i)
            {
                int32_t index = current_recs[i]->file_index;

                if (filter_exists[index])
                {
                    if (!filters[index].apply(current_recs[i]->h,current_recs[i]->v,&variant))
                    {
                        ++no_filtered;
                        continue;
                    }
                }

                ++presence[index];
            }
            
            if (no_filtered==current_recs.size())
            {
                continue;
            }

            //annotate
            for (size_t i=0; i<no_files; ++i)
            {
                int32_t counts = 0;

                for (size_t j=0; j<no_files; ++j)
                {
                    if (presence[i] && !presence[j])
                    {
                        ++stats[i][j].a;
                        stats[i][j].a_ts += variant.ts;
                        stats[i][j].a_tv += variant.tv;
                        stats[i][j].a_ins += variant.ins;
                        stats[i][j].a_del += variant.del;
                    }
                    else if (presence[i] && presence[j])
                    {
                        ++stats[i][j].ab;
                        stats[i][j].ab_ts += variant.ts;
                        stats[i][j].ab_tv += variant.tv;
                        stats[i][j].ab_ins += variant.ins;
                        stats[i][j].ab_del += variant.del;
                    }
                    else if (!presence[i] && presence[j])
                    {
                        ++stats[i][j].b;
                        stats[i][j].b_ts += variant.ts;
                        stats[i][j].b_tv += variant.tv;
                        stats[i][j].b_ins += variant.ins;
                        stats[i][j].b_del += variant.del;
                    }
                    else
                    {
                        //not in either, do nothing
                    }

                    counts += presence[j];
                }

                if (counts && presence[i]==counts)
                {
                    ++stats[i][i].e;
                    stats[i][i].e_ts += variant.ts;
                    stats[i][i].e_tv += variant.tv;
                    stats[i][i].e_ins += variant.ins;
                    stats[i][i].e_del += variant.del;
                }
            }
        }
    };

    void print_options()
    {
        std::clog << "xcmp v" << version << "\n\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File                 " << "" << "\n";
        print_str_op("         [f] filter                         ", fexp);
        print_str_op("         [x] output tabulate directory      ", output_tabulate_dir);
        print_str_op("         [y] output pdf file                ", output_pdf_file);
        print_int_op("         [i] intervals                      ", intervals);
        std::clog << "\n";
   }

    void print_pdf()
    {
        if (output_tabulate_dir=="") return;

        append_cwd(output_tabulate_dir);

        //generate file
        int32_t ret = mkdir(output_tabulate_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


        std::string filepath = output_tabulate_dir + "/tabulate.tex";
        FILE *out = fopen(filepath.c_str(), "w");

        std::string g2s[3] = {"R/R","R/A","A/A"};

//        fprintf(out, "\\PassOptionsToPackage{table}{xcolor}\n");
//        fprintf(out, "\\documentclass{beamer}\n");
//        fprintf(out, "\\begin{document}\n");
//        fprintf(out, "\n");
//        fprintf(out, "\\begin{frame}{Data set summary}\n");
//        fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
//        fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
//        fprintf(out, "\\begin{tabular}{rrrr}\n");
//        fprintf(out, "\\rowcolor{blue!50}\n");
//        fprintf(out, "No. Indels & Ins/Del & Insertions & Deletions \\\\ \n");
//        fprintf(out, "%d & %.1f & %d & %d \\\\ \n", stats[0].a, (float)stats[0].a_ins/(stats[0].a_del), stats[0].a_ins, stats[0].a_del);
//        fprintf(out, "\\end{tabular}}\n");
//        fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
//        fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
//        fprintf(out, "\\begin{tabular}{rrr}\n");
//        fprintf(out, "\\rowcolor{blue!50}\n");
//        fprintf(out, "Frameshift Indel Proportion (\\%%) & FS & NFS \\\\ \n");
//        //fprintf(out, "%.2f & %d & %d \\\\ \n", (float)fs/(fs+nfs), fs, nfs);
//        fprintf(out, "\\end{tabular}}\n");
//        fprintf(out, "\\end{frame}\n");
//
//        for (size_t i=1; i<dataset_labels.size(); ++i)
//        {
//            fprintf(out, "\n");
//            fprintf(out, "\\begin{frame}{Data set summary}\n");
//            fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
//            fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
//            fprintf(out, "\\begin{tabular}{rrrrr}\n");
//            fprintf(out, "\\rowcolor{blue!50}\n");
//            fprintf(out, "%s & no. indels & ins/del & ins & del\\\\ \n", dataset_labels[i].c_str());
//            fprintf(out, "A-B & %d & %.1f & %d & %d\\\\ \n",  stats[i].a, (float)stats[i].a_ins/(stats[i].a_del), stats[i].a_ins, stats[i].a_del);
//            fprintf(out, "A\\&B & %d & %.1f & %d & %d\\\\ \n",  stats[i].ab, (float)stats[i].ab_ins/(stats[i].ab_del), stats[i].ab_ins, stats[i].ab_del);
//            fprintf(out, "B-A & %d & %.1f & %d & %d\\\\ \n",  stats[i].b, (float)stats[i].b_ins/(stats[i].b_del), stats[i].b_ins, stats[i].b_del);
//            fprintf(out, " &  &  & &  \\\\ \n");
//            fprintf(out, " Precision & %.2f\\%% &  &  & \\\\ \n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
//            fprintf(out, " Sensitivity & %.2f\\%% &  &  &  \\\\ \n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
//            fprintf(out, "\\end{tabular}}\n");
//            fprintf(out, "\\end{frame}\n");
//        }

        fprintf(out, "\n");
        fprintf(out, "\\end{document}\n");

        fclose(out);

        std::string cmd = "cd "  + output_tabulate_dir + "; pdflatex tabulate.tex > run.log; mv tabulate.pdf " + output_pdf_file;
        std::cerr << cmd << "\n";
        int32_t sys_ret = system(cmd.c_str());
    };

    void print_stats()
    {
        //print overlap
        fprintf(stderr, "\noverlap");
        for (size_t j=0; j<dataset_labels.size(); ++j)
        {
            fprintf(stderr, "    %6s", dataset_labels[j].c_str());
        }
        fprintf(stderr, "    %6s", "ex");
        fprintf(stderr, "    n");
        fprintf(stderr, "\n");

        for (size_t i=0; i<dataset_labels.size(); ++i)
        {
            fprintf(stderr, "%7s", dataset_labels[i].c_str());
            for (size_t j=0; j<dataset_labels.size(); ++j)
            {
                fprintf(stderr, "    %5.1f%%", (float)stats[i][j].ab/(stats[i][j].ab+stats[i][j].b)*100);
            }
            fprintf(stderr, "    %5.1f%%", (float)stats[i][i].e/(stats[i][i].ab)*100);
            fprintf(stderr, "    %d", stats[i][i].ab);
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
        
        //print ts/tv ratios
        fprintf(stderr, "ts/tv  ");
        for (size_t j=0; j<dataset_labels.size(); ++j)
        {
            fprintf(stderr, "    %6s", dataset_labels[j].c_str());
        }
        fprintf(stderr, "    %6s", "ex");
        fprintf(stderr, "    n");
        fprintf(stderr, "\n");

        for (size_t i=0; i<dataset_labels.size(); ++i)
        {
            fprintf(stderr, "%7s", dataset_labels[i].c_str());
            for (size_t j=0; j<dataset_labels.size(); ++j)
            {
                fprintf(stderr, "    %6.2f", (float)stats[i][j].ab_ts/stats[i][j].ab_tv);
            }
            fprintf(stderr, "    %6.2f", (float)stats[i][i].e_ts/(stats[i][i].e_tv));
            fprintf(stderr, "    %d", stats[i][i].ab);
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
        
        //print ins/del ratios
        fprintf(stderr, "ins/del");
        for (size_t j=0; j<dataset_labels.size(); ++j)
        {
            fprintf(stderr, "    %6s", dataset_labels[j].c_str());
        }
        fprintf(stderr, "    %6s", "ex");
        fprintf(stderr, "    n");
        fprintf(stderr, "\n");

        for (size_t i=0; i<dataset_labels.size(); ++i)
        {
            fprintf(stderr, "%7s", dataset_labels[i].c_str());
            for (size_t j=0; j<dataset_labels.size(); ++j)
            {
                fprintf(stderr, "    %6.2f", (float)stats[i][j].ab_ins/stats[i][j].ab_del);
            }
            fprintf(stderr, "    %6.2f", (float)stats[i][i].e_ins/(stats[i][i].e_del));
            fprintf(stderr, "    %d", stats[i][i].ab);
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
    };

    ~Igor()
    {
    };

    private:
};

}

void xcmp(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.xcmp();
    igor.print_stats();
    igor.print_pdf();
}
