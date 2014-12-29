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

#include "profile_na12878.h"

namespace
{

/**
 * Specific to Broad's KB.
 */
class BroadKBStats
{
    public:

    int32_t tp, fp, tn, fn;

    BroadKBStats()
    {
        tp = 0;
        fp = 0;
        tn = 0;
        fn = 0;
    };
};

class OverlapStats
{
    public:

    int32_t a,ab,b,a_ts,ab_ts,b_ts,a_tv,ab_tv,b_tv,a_ins,ab_ins,b_ins,a_del,ab_del,b_del;

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

class ConcordanceStats
{
    public:

    int32_t geno[4][4];

    ConcordanceStats()
    {
        for (int32_t i=0; i<4; ++i)
        {
            for (int32_t j=0; j<4; ++j)
            {
                geno[i][j] = 0;
            }
        }
    };
};

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::vector<std::string> input_vcf_files;
    std::string ref_fasta_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////////////////////
    //reference data sets//
    ///////////////////////
    std::string ref_data_sets_list;
    std::vector<std::string> dataset_labels;
    std::vector<std::string> dataset_types;
    std::vector<std::string> dataset_fexps;
    std::string cds_bed_file;

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
    BroadKBStats kbstats;
    std::vector<OverlapStats> stats;
    std::vector<ConcordanceStats> concordance;
    int32_t no_variants;
    int32_t no_positive_variants;
    int32_t no_negative_variants;
    int32_t ts;
    int32_t tv;
    int32_t syn;
    int32_t nsyn;
    int32_t ins;
    int32_t del;
    int32_t fs;
    int32_t nfs;

    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;
    OrderedRegionOverlapMatcher *orom_gencode_cds;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "profile NA12878";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_data_sets_list("g", "g", "file containing list of reference datasets []", false, "", "file", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            ref_fasta_file = arg_ref_fasta_file.getValue();
            fexp = arg_fexp.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_data_sets_list = arg_ref_data_sets_list.getValue();
            input_vcf_file = arg_input_vcf_file.getValue();

            ///////////////////////
            //parse input VCF files
            ///////////////////////
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
//#
//# dataset - name of data set, this label will be printed.
//# type    - True Positives (TP) and False Positives (FP)
//#           overlap percentages labeled as (Precision, Sensitivity) and (False Discovery Rate, Type I Error) respectively
//#         - annotation
//#           file is used for GENCODE annotation of frame shift and non frame shift Indels
//# filter  - filter applied to variants for this particular data set
//# path    - path of indexed BCF file
//#dataset              type         filter    path
//broad.kb              TP           PASS      /net/fantasia/home/atks/dev/vt/bundle/public/grch37/broad.kb.241365variants.genotypes.bcf
//illumina.platinum     TP           PASS      /net/fantasia/home/atks/dev/vt/bundle/public/grch37/NA12878.illumina.platinum.5284448variants.genotypes.bcf
//gencode.v19           annotation   .         /net/fantasia/home/atks/dev/vt/bundle/public/grch37/gencode.v19.annotation.gtf.gz

        input_vcf_files.push_back(input_vcf_file);
        dataset_labels.push_back("data");
        dataset_types.push_back("ref");
        dataset_fexps.push_back(fexp);

        htsFile *hts = hts_open(ref_data_sets_list.c_str(), "r");
        kstring_t s = {0,0,0};
        std::vector<std::string> vec;
        while (hts_getline(hts, '\n', &s)>=0)
        {
            if (s.s[0] == '#')
                continue;

            std::string line(s.s);
            split(vec, " ", line);

            if (vec[1] == "TP" || vec[1] == "FP")
            {
                dataset_labels.push_back(vec[0]);
                dataset_types.push_back(vec[1]);
                dataset_fexps.push_back(vec[2]);
                input_vcf_files.push_back(vec[3]);
            }
            else if (vec[1] == "cds_annotation")
            {
                cds_bed_file = vec[3];
            }
            else
            {
                std::cerr << "Reference data set type: \"" << vec[1] << "\" not recognised\n";
                exit(1);
            }
        }
        hts_close(hts);
        if (s.m) free(s.s);

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
        sr = new BCFSyncedReader(input_vcf_files, intervals, false);

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip(ref_fasta_file);
        orom_gencode_cds = new OrderedRegionOverlapMatcher(cds_bed_file);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_positive_variants = 0;
        no_negative_variants = 0;
        fs = 0;
        nfs = 0;
    }

    void profile_chm1()
    {
        //for combining the alleles
        std::vector<bcfptr*> current_recs;
        std::vector<Interval*> overlaps;
        Variant variant;
        int32_t no_overlap_files = input_vcf_files.size();
        std::vector<int32_t> presence(no_overlap_files, 0);
        std::vector<bcfptr*> presence_bcfptr(no_overlap_files, NULL);
        stats.resize(no_overlap_files);
        concordance.resize(no_overlap_files);

        int32_t *gts = NULL;
        int32_t n = 0;
        int32_t x1, x2, x, k0;

        int32_t na12878_index[no_overlap_files];

        //get NA12878 gt position from reference file
        for (int32_t i=0; i<no_overlap_files; ++i)
        {
            na12878_index[i] = bcf_hdr_id2int(sr->hdrs[i], BCF_DT_SAMPLE, "NA12878");
            if (na12878_index[i]==-1)
            {
                fprintf(stderr, "[E:%s:%d %s] NA12878 not found in %s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files[i].c_str());
                exit(1);
            }
        }

        int32_t discordance_filter = 0;

        while(sr->read_next_position(current_recs))
        {
            //check first variant
            bcf1_t *v = current_recs[0]->v;
            bcf_hdr_t *h = current_recs[0]->h;
            int32_t vtype = vm->classify_variant(h, v, variant);

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
                presence_bcfptr[index] = current_recs[i];
            }

            std::string chrom = bcf_get_chrom(h,v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);

            //annotate
            if (presence[0])
            {
                if (orom_gencode_cds->overlaps_with(chrom, start1, end1))
                {
                    if (abs(variant.alleles[0].dlen)%3!=0)
                    {
                        ++fs;
                    }
                    else
                    {
                        ++nfs;
                    }
                }
            }

            int32_t ts = variant.alleles[0].ts;
            int32_t tv = variant.alleles[0].tv;
            int32_t ins = variant.alleles[0].dlen ? variant.alleles[0].ins : 0;
            int32_t del = variant.alleles[0].dlen ? 1-ins : 0;

            if (presence[0])
            {
                ++stats[0].a;
                stats[0].a_ts += ts;
                stats[0].a_tv += tv;
                stats[0].a_ins += ins;
                stats[0].a_del += del;

                bcf_unpack(v, BCF_UN_STR);
                int k = bcf_get_genotypes(presence_bcfptr[0]->h, presence_bcfptr[0]->v, &gts, &n);
                x1 = bcf_gt_allele(gts[na12878_index[0]*2]);
                x2 = bcf_gt_allele(gts[na12878_index[0]*2+1]);

                x = x1+x2;

                if (x<0)
                {
                    x = 3;
                }
            }

            //update overlap stats
            for (size_t i=1; i<no_overlap_files; ++i)
            {
                //******************
                //for BROAD KB stats
                //******************
                if (i==1)
                {
                    int32_t y1 = -1;
                    int32_t y2 = -1;
                    int32_t y = -2;
                    int32_t xt = -2;

                    char* dst = NULL;
                    int32_t ndst = 0;

                    if (presence[1])
                    {
                        bcf_unpack(presence_bcfptr[i]->v, BCF_UN_IND);
                        k0 = bcf_get_genotypes(presence_bcfptr[1]->h, presence_bcfptr[1]->v, &gts, &n);

                        bcf_get_info_string(presence_bcfptr[1]->h, presence_bcfptr[1]->v, "TruthStatus", &dst, &ndst);

                        y1 = bcf_gt_allele(gts[na12878_index[i]*2]);
                        y2 = bcf_gt_allele(gts[na12878_index[i]*2+1]);
                        if (y1>0 || y2>0)
                        {
                            y = 1;
                        }
                        else if (y1==0 || y2==0)
                        {
                            y = 0;
                        }
                        else if (y1<0 || y2<0)
                        {
                            y = -1;
                        }
                    }

                    if (presence[0])
                    {
                        if (x1>0 || x2>0)
                        {
                            xt = 1;
                        }
                        else if (x1==0 || x2==0)
                        {
                            xt = 0;
                        }
                        else if (x1<0 || x2<0)
                        {
                            xt = -1;
                        }
                        ++no_variants;

                        if (xt>0) ++no_positive_variants;
                        if (xt==0) ++no_negative_variants;
                    }

                    if (presence[0] && presence[1])
                    {
                        if (strcmp(dst, "TRUE_POSITIVE")==0)
                        {
                            if (xt>0 && y>0)
                            {
                                ++kbstats.tp;
                            }
                            else if (xt>0 && y==0)
                            {
                                ++kbstats.fp;
                            }
                            else if (xt==0 && y>0)
                            {
                                ++kbstats.fn;
                            }
                        }
                        else if (strcmp(dst, "FALSE_POSITIVE")==0)
                        {
                            if (xt>0)
                            {
                                ++kbstats.fp;
                            }
                            else if (xt==0)
                            {
                                ++kbstats.tn;
                            }
                        }

                    }
                    else if (presence[0] && !presence[1])
                    {
                        //++kbstats.fp; potentially false positive BUT knowledge base is curated manually.
                    }
                    else if (!presence[0] && presence[1])
                    {
                        if (strcmp(dst, "TRUE_POSITIVE")==0)
                        {
                            if (y>0)
                            {
                                ++kbstats.fn;
                            }
                            else if (y==0)
                            {
                                ++kbstats.tn;
                            }
                        }
                        else if (strcmp(dst, "FALSE_POSITIVE")==0)
                        {
                            ++kbstats.tn;
                        }
                    }
                    else
                    {
                        //++kbstats.tn; potentially true negative BUT knowledge base is curated manually.
                    }

                    if (ndst) free(dst);
                }

                if (presence[0] && !presence[i])
                {
                    ++stats[i].a;
                    stats[i].a_ts += ts;
                    stats[i].a_tv += tv;
                    stats[i].a_ins += ins;
                    stats[i].a_del += del;
                }
                else if (presence[0] && presence[i])
                {
                    ++stats[i].ab;
                    stats[i].ab_ts += ts;
                    stats[i].ab_tv += tv;
                    stats[i].ab_ins += ins;
                    stats[i].ab_del += del;

                    bcf_unpack(presence_bcfptr[i]->v, BCF_UN_IND);
                    int k = bcf_get_genotypes(presence_bcfptr[i]->h, presence_bcfptr[i]->v, &gts, &n);

                    int32_t y1 = bcf_gt_allele(gts[na12878_index[i]*2]);
                    int32_t y2 = bcf_gt_allele(gts[na12878_index[i]*2+1]);
                    int32_t y = y1+y2;
                    if (y<0)
                    {
                        y = 3;
                    }

                    ++concordance[i].geno[x][y];
                }
                else if (!presence[0] && presence[i])
                {
                    ++stats[i].b;
                    stats[i].b_ts += ts;
                    stats[i].b_tv += tv;
                    stats[i].b_ins += ins;
                    stats[i].b_del += del;
                }
                else
                {
                    //not in either, do nothing
                }

                presence[i]=0;
                presence_bcfptr[i]=NULL;
            }

            presence[0] = 0;
        }
    };

    void print_options()
    {
        std::clog << "profile_na12878 v" << version << "\n\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File                 " << input_vcf_file << "\n";
        print_str_op("         [f] filter                         ", fexp);
        std::clog << "         [g] reference data sets list file  " << ref_data_sets_list << "\n";
        std::clog << "         [r] reference FASTA file           " << ref_fasta_file << "\n";
        print_int_op("         [i] intervals                      ", intervals);
        std::clog << "\n\n";
   }

    void print_stats()
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "  %s\n", "data set");
        fprintf(stderr, "  No Variants : %10d [%.2f] [%.2f]\n", stats[0].a, (float)stats[0].a_ts/stats[0].a_tv, (float)stats[0].a_ins/(stats[0].a_del));
        fprintf(stderr, "       FS/NFS : %10.2f (%d/%d)\n", (float)fs/(fs+nfs), fs, nfs);
        fprintf(stderr, "\n");

        fprintf(stderr, "  Broad KB\n");
        fprintf(stderr, "    TP  %10d \n", kbstats.tp);
        fprintf(stderr, "    FP  %10d \n", kbstats.fp);
        fprintf(stderr, "    TN  %10d \n", kbstats.tn);
        fprintf(stderr, "    FN  %10d \n", kbstats.fn);
        fprintf(stderr, "    P   %10d \n", kbstats.tp+kbstats.fp);
        fprintf(stderr, "    N   %10d \n", kbstats.tn+kbstats.fn);
        fprintf(stderr, "    +   %10d \n", no_positive_variants);
        fprintf(stderr, "    -   %10d \n", no_negative_variants);
        fprintf(stderr, "  TDR %5.2f (TP/(TP+FP))\n", (float)kbstats.tp/(kbstats.tp+kbstats.fp)*100);
        fprintf(stderr, "  FNR %5.2f (FN/(TN+FN))\n", (float)kbstats.fn/(kbstats.fn+kbstats.tn)*100);
        fprintf(stderr, "\n");


        for (int32_t i=1; i<dataset_labels.size(); ++i)
        {
            fprintf(stderr, "  %s\n", dataset_labels[i].c_str());
            fprintf(stderr, "                   ts/tv  ins/del\n");
            fprintf(stderr, "    A-B %10d [%.2f] [%.2f]\n", stats[i].a,  (float)stats[i].a_ts/stats[i].a_tv, (float)stats[i].a_ins/stats[i].a_del);
            fprintf(stderr, "    A&B %10d [%.2f] [%.2f]\n", stats[i].ab, (float)stats[i].ab_ts/stats[i].ab_tv, (float)stats[i].ab_ins/stats[i].ab_del);
            fprintf(stderr, "    B-A %10d [%.2f] [%.2f]\n", stats[i].b,  (float)stats[i].b_ts/stats[i].b_tv, (float)stats[i].b_ins/stats[i].b_del);

            if (dataset_types[i]=="TP")
            {
                fprintf(stderr, "    Precision    %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
                fprintf(stderr, "    Sensitivity  %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
            }
            else
            {
                fprintf(stderr, "    FDR          %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
                fprintf(stderr, "    Type I Error %4.1f%%\n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
            }
            fprintf(stderr, "\n");
        }

        for (int32_t i=1; i<dataset_labels.size(); ++i)
        {
            int32_t (&geno)[4][4] = concordance[i].geno;
            int32_t total = 0;
            int32_t concordance = 0;
            for (int32_t i=0; i<3; ++i)
            {
                for (int32_t j=0; j<3; ++j)
                {
                    total += geno[i][j];
                    if (i==j)
                    {
                        concordance += geno[i][j];
                    }
                }
            }

            int32_t discordance = total-concordance;

            fprintf(stderr, "  %s\n", dataset_labels[i].c_str());
            fprintf(stderr, "                R/R       R/A       A/A       ./.\n");
            fprintf(stderr, "    R/R    %8d  %8d  %8d  %8d\n", geno[0][0], geno[0][1], geno[0][2], geno[0][3]);
            fprintf(stderr, "    R/A    %8d  %8d  %8d  %8d\n", geno[1][0], geno[1][1], geno[1][2], geno[1][3]);
            fprintf(stderr, "    A/A    %8d  %8d  %8d  %8d\n", geno[2][0], geno[2][1], geno[2][2], geno[2][3]);
            fprintf(stderr, "    ./.    %8d  %8d  %8d  %8d\n", geno[3][0], geno[3][1], geno[3][2], geno[3][3]);
            fprintf(stderr, "\n");

            fprintf(stderr, "    Total genotype pairs :  %8d\n", total);
            fprintf(stderr, "    Concordance          :  %8.2f%% (%d)\n", (float)concordance/total*100, concordance);
            fprintf(stderr, "    Discordance          :  %8.2f%% (%d)\n", (float)discordance/total*100, discordance);

            fprintf(stderr, "\n");
        }
    };

    ~Igor()
    {
    };

    private:
};

}

void profile_chm1(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_chm1();
    igor.print_stats();
}
