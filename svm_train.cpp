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

#include "profile_snps.h"

namespace
{
class OverlapStats
{
    public:

    uint32_t a,ab,b,a_ts,ab_ts,b_ts,a_tv,ab_tv,b_tv;

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
    std::string ref_data_sets_list;
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
    kstring_t line;

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
    std::vector<OverlapStats> stats;
    int32_t no_snps;
    int32_t nonsyn;
    int32_t syn;
    int32_t lcplx;

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
            std::string desc = "profile SNPs";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_ref_data_sets_list("g", "g", "file containing list of reference datasets []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "VTYPE==SNP", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            ref_fasta_file = arg_ref_fasta_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
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
//#dataset               type             filter                                 path
//1000g                  TP               N_ALLELE==2&&VTYPE==SNP                /net/fantasia/home/atks/ref/vt/grch37/1000G.v5.snps.indels.complex.svs.sites.bcf
//dbsnp                  TP               N_ALLELE==2&&VTYPE==SNP                /net/fantasia/home/atks/ref/vt/grch37/dbSNP138.snps.indels.complex.sites.bcf
//GENCODE_V19            cds_annotation   .                                      /net/fantasia/home/atks/ref/vt/grch37/gencode.v19.cds.bed.gz
//DUST                   cplx_annotation  .                                      /net/fantasia/home/atks/ref/vt/grch37/mdust.bed.gz

        input_vcf_files.push_back(input_vcf_file);
        dataset_labels.push_back("data");
        dataset_types.push_back("ref");
        dataset_fexps.push_back(fexp);

        htsFile *hts = hts_open(ref_data_sets_list.c_str(), "r");
        if (!hts)
        {
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

            if (vec[1] == "TP" || vec[1] == "FP")
            {
                dataset_labels.push_back(vec[0]);
                dataset_types.push_back(vec[1]);
                dataset_fexps.push_back(vec[2]);
                input_vcf_files.push_back(vec[3]);;
            }
            else if (vec[1] == "cds_annotation")
            {
                cds_bed_file = vec[3];
            }
            else if (vec[1] == "cplx_annotation")
            {
                cplx_bed_file = vec[3];
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
        orom_lcplx = new OrderedRegionOverlapMatcher(cplx_bed_file);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_snps = 0;
        nonsyn = 0;
        syn = 0;
    }

    void profile_snps()
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

//            for (size_t i=0; i<no_overlap_files; ++i)
//                presence[i]=0;
//                
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
                std::string chrom = bcf_get_chrom(h,v);
                int32_t start1 = bcf_get_pos1(v);
                int32_t end1 = bcf_get_end_pos1(v);
                
                if (orom_lcplx->overlaps_with(chrom, start1, end1))
                {
                    ++lcplx;
                }

                if (orom_gencode_cds->overlaps_with(chrom, start1, end1))
                {
                    ++nonsyn;
                }

                ++no_snps;
            }

            int32_t ts = variant.alleles[0].ts;
            int32_t tv = variant.alleles[0].tv;

            if (presence[0])
            {
                ++stats[0].a;
                stats[0].a_ts += ts;
                stats[0].a_tv += tv;
            }

            //update overlap stats
            for (size_t i=1; i<no_overlap_files; ++i)
            {
                if (presence[0] && !presence[i])
                {
                    ++stats[i].a;
                    stats[i].a_ts += ts;
                    stats[i].a_tv += tv;
                }
                else if (presence[0] && presence[i])
                {
                    ++stats[i].ab;
                    stats[i].ab_ts += ts;
                    stats[i].ab_tv += tv;
                }
                else if (!presence[0] && presence[i])
                {
                    ++stats[i].b;
                    stats[i].b_ts += ts;
                    stats[i].b_tv += tv;
                }
                else
                {
                    //not in either, do nothing
                }

                presence[i]=0;
            }

            presence[0] = 0;
        }
    };

    void print_options()
    {
        std::clog << "profile_snps v" << version << "\n\n";
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
        fprintf(stderr, "     No. SNPs          : %10d [%.2f]\n", stats[0].a,  (float)stats[0].a_ts/(stats[0].a_tv));
        //fprintf(stderr, "        SYN/NONSYN     : %10.2f (%d/%d)\n", (float)nonsyn/(nonsyn+syn), nonsyn, syn);
        fprintf(stderr, "        Low complexity : %10.2f (%d/%d)\n", (float)lcplx/stats[0].a, lcplx, stats[0].a);
        fprintf(stderr, "\n");

        for (int32_t i=1; i<dataset_labels.size(); ++i)
        {
            fprintf(stderr, "  %s\n", dataset_labels[i].c_str());
            fprintf(stderr, "    A-B %10d [%.2f]\n", stats[i].a,  (float)stats[i].a_ts/(stats[i].a_tv));
            fprintf(stderr, "    A&B %10d [%.2f]\n", stats[i].ab, (float)stats[i].ab_ts/stats[i].ab_tv);
            fprintf(stderr, "    B-A %10d [%.2f]\n", stats[i].b,  (float)stats[i].b_ts/(stats[i].b_tv));

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
    };

    ~Igor()
    {
    };

    private:
};

}

void profile_snps(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_snps();
    igor.print_stats();
}
