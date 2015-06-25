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

#include "merge.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_vcf_files;
    std::string input_vcf_file_list;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    bool print;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    ///////////////
    //general use//
    ///////////////
    kstring_t variant;

    /////////
    //stats//
    /////////
    uint32_t no_samples;
    uint32_t no_candidate_snps;
    uint32_t no_candidate_indels;
    uint32_t no_candidate_snpindels;
    uint32_t no_other_variant_types;

    /////////
    //tools//
    /////////
    LogTool *lt;
    VariantManip * vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Merges VCF files. If duplicate variants are observed for a sample, they are dropped.  \n\
                Assumes that genotype fields are the same for all individuals for each record\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);

            cmd.parse(argc, argv);

            input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            const std::vector<std::string>& v = arg_input_vcf_files.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print = arg_print.getValue();

            if (input_vcf_files.size()==0)
            {
                fprintf(stderr, "[E:%s:%d %s] no input vcf files.\n", __FILE__, __LINE__, __FUNCTION__);
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
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bcf_hdr_append(odw->hdr, "##fileformat=VCFv4.1");
        bcf_hdr_transfer_contigs(sr->hdrs[0], odw->hdr);

        ///////////////
        //general use//
        ///////////////
        variant = {0,0,0};
        no_samples = sr->nfiles;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
        no_candidate_indels = 0;
        no_candidate_snpindels = 0;
        no_other_variant_types = 0;

        /////////
        //tools//
        /////////
        lt = new LogTool();
        vm = new VariantManip();
    }

    void merge()
    {
        int32_t nfiles = sr->get_nfiles();

        //add all sample names to output vcf header and yell if there are more than one occurence of a sample
        int32_t no_samples = 0;
        for (int32_t i=0; i<nfiles; ++i)
        {
            for (int32_t j=0; j<bcf_hdr_nsamples(sr->hdrs[i]); ++j)
            {
                if (bcf_hdr_id2int(odw->hdr, BCF_DT_SAMPLE, bcf_hdr_get_sample_name(sr->hdrs[i], j))==-1)
                {
                   // std::cerr << "adding " << bcf_hdr_get_sample_name(sr->hdrs[i], j) << "\n";

                    bcf_hdr_add_sample(odw->hdr, bcf_hdr_get_sample_name(sr->hdrs[i], j));
                    ++no_samples;
                }
                else
                {
                    fprintf(stderr, "[E:%s:%d %s] %s is present more than once.\n", __FILE__, __LINE__, __FUNCTION__, bcf_hdr_get_sample_name(sr->hdrs[i], j));
                    exit(1);
                }
            }
        }
        bcf_hdr_add_sample(odw->hdr, NULL);

        if (no_samples)
        {
            bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            bcf_hdr_append(odw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes\">");
            bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");
            bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allele Depth\">");
            bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
        }
        odw->write_hdr();

        std::vector<bcfptr*> current_recs;
        int32_t *cgt = (int32_t*) malloc(no_samples*2*sizeof(int32_t));
        int32_t *pls = (int32_t*) malloc(no_samples*3*sizeof(int32_t));
        int32_t *dps = (int32_t*) malloc(no_samples*sizeof(int32_t));
        int32_t *ads = (int32_t*) malloc(no_samples*3*sizeof(int32_t));
        int32_t *gqs = (int32_t*) malloc(no_samples*sizeof(int32_t));
        int ncount =0;

        std::vector<bcfptr*> sample2record(no_samples, NULL);
        std::vector<int32_t> sample2index(no_samples, -1);

        while(sr->read_next_position(current_recs))
        {
            if (no_samples)
            {
                bcf1_t *v = current_recs[0]->v;
                bcf_hdr_t *h = current_recs[0]->h;
                bcf1_t *nv = odw->get_bcf1_from_pool();
                bcf_set_chrom(odw->hdr, nv, bcf_get_chrom(h, v));
                bcf_set_pos1(nv, bcf_get_pos1(v));
                bcf_update_alleles(odw->hdr, nv, const_cast<const char**>(bcf_get_allele(v)), bcf_get_n_allele(v));
                bcf_set_n_sample(nv, no_samples);

                int32_t ngt = 0;
                for (size_t i=0; i<current_recs.size(); ++i)
                {
                    int32_t file_index = current_recs[i]->file_index;
                    bcf1_t *v = current_recs[i]->v;
                    bcf_hdr_t *h = current_recs[i]->h;

                    int32_t *gt = NULL;
                    int32_t n = 0;
                    int32_t ploidy = bcf_get_genotypes(h, v, &gt, &n);

                    ploidy /= bcf_hdr_nsamples(h);

                    int32_t *pl = NULL;
                    int32_t n_pl=0;
                    bcf_get_format_int32(h, v, "PL", &pl, &n_pl);

                    int32_t *dp = NULL;
                    int32_t n_dp=0;
                    bcf_get_format_int32(h, v, "DP", &dp, &n_dp);

                    int32_t *ad = NULL;
                    int32_t n_ad=0;
                    bcf_get_format_int32(h, v, "AD", &ad, &n_ad);

                    int32_t *gq = NULL;
                    int32_t n_gq=0;
                    bcf_get_format_int32(h, v, "GQ", &gq, &n_gq);

                    for (int32_t j=0; j<bcf_hdr_nsamples(h); ++j)
                    {
                        int32_t k = bcf_hdr_id2int(odw->hdr, BCF_DT_SAMPLE, bcf_hdr_get_sample_name(sr->hdrs[file_index], j));
                        cgt[k*2] = gt[j*2];
                        cgt[k*2+1] = gt[j*2+1];

                        if (n_pl)
                        {
                            pls[k*3] = pl[j*3];
                            pls[k*3+1] = pl[j*3+1];
                            pls[k*3+2] = pl[j*3+2];
                        }
                        else
                        {
                            pls[k*3] = bcf_int32_missing;
                            pls[k*3+1] = bcf_int32_vector_end;
                        }

                        if (n_dp)
                        {
                            dps[k] = dp[j];
                        }
                        else
                        {
                            dps[k] = bcf_int32_missing;
                        }

                        if (n_ad)
                        {
                            ads[k*3] = ad[j*3];
                            ads[k*3+1] = ad[j*3+1];
                            ads[k*3+2] = ad[j*3+2];
                        }
                        else
                        {
                            ads[k*3] = bcf_int32_missing;
                            ads[k*3+1] = bcf_int32_vector_end;
                        }

                        if (n_gq)
                        {
                            gqs[k] = gq[j];
                        }
                        else
                        {
                            gqs[k] = bcf_int32_missing;
                        }

                        ++ngt;
                    }

                    free(gt); //is this necessary - check later
                    free(pl);
                    free(dp);
                }

                bcf_update_genotypes(odw->hdr,nv,cgt,ngt*2);
                bcf_update_format_int32(odw->hdr,nv,"PL",pls,ngt*3);
                bcf_update_format_int32(odw->hdr,nv,"DP",dps,ngt);
                bcf_update_format_int32(odw->hdr,nv,"AD",ads,ngt*3);
                bcf_update_format_int32(odw->hdr,nv,"GQ",gqs,ngt);

                odw->write(nv);
            }
            else
            {
                odw->write(current_recs[0]->v);
            }
        }

        free(cgt);

        sr->close();
        odw->close();
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "merge v" << version << "\n\n";
        print_ifiles("options:     input VCF file        ", input_vcf_files);
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::cerr << "stats: Total Number of Candidate SNPs                 " << no_candidate_snps << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool merge(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.merge();
    igor.print_stats();
    return igor.print;
}

