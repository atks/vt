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

    int32_t a,ab,b,a_ins,ab_ins,b_ins,a_del,ab_del,b_del;

    OverlapStats()
    {
        a = 0;
        ab = 0;
        b = 0;

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
    std::vector<std::string> input_vcf_files;
    std::string input_vcf_file_list;
    std::string output_tabulate_dir;
    std::string output_pdf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    //////////////////////////////////////////////
    //reference file info : to store in an object?
    //////////////////////////////////////////////
    std::vector<std::string> dataset_labels;
    std::vector<std::string> dataset_types;
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
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "VTYPE==INDEL", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_tabulate_dir("x", "x", "output latex directory []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output pdf file []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_filters(dataset_fexps, arg_fexp.getValue());
            output_tabulate_dir = arg_output_tabulate_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        dataset_labels.push_back("data");
        dataset_types.push_back("ref");
        dataset_fexps.push_back(fexp);

        /////////////////////////
        //filter initialization//
        /////////////////////////
        for (size_t i=0; i<dataset_fexps.size(); ++i)
        {
            filters.push_back(Filter(dataset_fexps[i]));
            filter_exists.push_back(dataset_fexps[i]!="");
        }
        no_filters = filters.size();

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
        int32_t no_overlap_files = input_vcf_files.size();
        std::vector<int32_t> presence(no_overlap_files, 0);
        stats.resize(no_overlap_files);

        while(sr->read_next_position(current_recs))
        {
            //check first variant
            bcf1_t *v = current_recs[0]->v;
            bcf_hdr_t *h = current_recs[0]->h;
            int32_t vtype = vm->classify_variant(h, v, variant);

            std::string chrom = bcf_get_chrom(h,v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);

            for (size_t i=0; i<no_overlap_files; ++i)
                presence[i]=0;

            //check existence
            for (size_t i=0; i<current_recs.size(); ++i)
            {
                int32_t index = current_recs[i]->file_index;

                if (filter_exists[index])
                {
                    if (!filters[index].apply(current_recs[i]->h,current_recs[i]->v,&variant))
                    {
                        continue;
                    }
                }

                ++presence[index];
            }

            //annotate
            if (presence[0])
            {



            }

            int32_t ins = variant.alleles[0].ins;
            int32_t del = 1-ins;

//            if (presence[0])
//            {
//                ++stats[0].a;
//                stats[0].a_ins += ins;
//                stats[0].a_del += del;
//            }

            //update overlap stats
//            for (size_t i=1; i<no_overlap_files; ++i)
//            {
//                if (presence[0] && !presence[i])
//                {
//                    ++stats[i].a;
//                    stats[i].a_ins += ins;
//                    stats[i].a_del += del;
//                }
//                else if (presence[0] && presence[i])
//                {
//                    ++stats[i].ab;
//                    stats[i].ab_ins += ins;
//                    stats[i].ab_del += del;
//                }
//                else if (!presence[0] && presence[i])
//                {
//                    ++stats[i].b;
//                    stats[i].b_ins += ins;
//                    stats[i].b_del += del;
//                }
//                else
//                {
//                    //not in either, do nothing
//                }
//
//                presence[i]=0;
//            }

            presence[0] = 0;
        }
    };

    void print_options()
    {
        std::clog << "cross_compare v" << version << "\n\n";
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
        fprintf(stderr, "\n");
        fprintf(stderr, "  %s\n", "data set");
        fprintf(stderr, "\n");

        for (size_t i=1; i<dataset_labels.size(); ++i)
        {
//            fprintf(stderr, "  %s\n", dataset_labels[i].c_str());
//            fprintf(stderr, "    A-B %10d [%.2f]\n", stats[i].a,  (float)stats[i].a_ins/(stats[i].a_del));
//            fprintf(stderr, "    A&B %10d [%.2f]\n", stats[i].ab, (float)stats[i].ab_ins/stats[i].ab_del);
//            fprintf(stderr, "    B-A %10d [%.2f]\n", stats[i].b,  (float)stats[i].b_ins/(stats[i].b_del));
//
//            if (dataset_types[i]=="TP")
//            {
//                fprintf(stderr, "    Precision    %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
//                fprintf(stderr, "    Sensitivity  %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
//            }
//            else
//            {
//                fprintf(stderr, "    FDR          %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
//                fprintf(stderr, "    Type I Error %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
//            }
//            fprintf(stderr, "\n");
        }
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
