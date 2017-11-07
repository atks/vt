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

#include "merge_genotypes.h"

#define SINGLE     0
#define AGGREGATED 1

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
    std::string candidate_sites_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    ///////////////
    //general use//
    ///////////////
    std::vector<int32_t> file_types;
    kstring_t variant;

    /////////
    //stats//
    /////////
    uint32_t no_samples;
    uint32_t no_candidate_snps;
    uint32_t no_candidate_indels;

    /////////
    //tools//
    /////////
    VariantManip * vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc =
"Merge GT, PL, DP, ADF and ADR across samples.\n\
Extracts only the naive genotypes based on best guess genotypes.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_candidate_sites_vcf_file("c", "c", "candidate sites VCF file with annotation [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", true, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);

            cmd.parse(argc, argv);

            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            candidate_sites_vcf_file = arg_candidate_sites_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
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
        input_vcf_files.insert(input_vcf_files.begin(), candidate_sites_vcf_file);
        sr = new BCFSyncedReader(input_vcf_files, intervals, false);

        fprintf(stderr, "[I:%s:%d %s] Initializaing %zd VCF files ...\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files.size());

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bcf_hdr_append(odw->hdr, "##fileformat=VCFv4.2");
        bcf_hdr_transfer_contigs(sr->hdrs[0], odw->hdr);
        bcf_hdr_append(odw->hdr, "##QUAL=Maximum variant score of the alternative allele likelihood ratio: -10 * log10 [P(Non variant)/P(Variant)] amongst all individuals.");

        //add samples to output merged file
        for (int32_t i=1; i<sr->hdrs.size(); ++i)
        {
            bcf_hdr_add_sample(odw->hdr, bcf_hdr_get_sample_name(sr->hdrs[i], 0) );
        }

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

        /////////
        //tools//
        /////////
        vm = new VariantManip();
    }

    void merge_genotypes()
    {
        int32_t *NSAMPLES = NULL;
        int32_t no_NSAMPLES = 0;
        int32_t *GT = NULL;
        int32_t *PL = NULL;
        int32_t *AD = NULL;
        int32_t *ADF = NULL;
        int32_t *ADR = NULL;
        int32_t *BQSUM = NULL;
        int32_t *DP = NULL;
        int32_t no_GT = 0;
        int32_t no_PL = 0;
        int32_t no_AD = 0;
        int32_t no_ADF = 0;
        int32_t no_ADR = 0;
        int32_t no_BQSUM = 0;
        int32_t no_DP = 0;

        bcf1_t* nv = bcf_init();
        Variant var;
        std::vector<bcfptr*> current_recs;

        while(sr->read_next_position(current_recs))
        {
            int32_t vtype;

//            std::cerr << "no recs " << current_recs.size() << "\n";

            //for each file
            for (uint32_t i=0; i<current_recs.size(); ++i)
            {
                int32_t file_index = current_recs[i]->file_index;
                bcf1_t *v = current_recs[i]->v;
                bcf_hdr_t *h = current_recs[i]->h;

//                bcf_print(h, v);

                if (i==0)
                {
                    //update variant information
                    bcf_clear(nv);
                    bcf_set_chrom(odw->hdr, nv, bcf_get_chrom(h, v));
                    bcf_set_pos1(nv, bcf_get_pos1(v));
                    bcf_update_alleles(odw->hdr, nv, const_cast<const char**>(bcf_get_allele(v)), bcf_get_n_allele(v));
                    vtype = vm->classify_variant(odw->hdr, nv, var);
                }

                float variant_score = bcf_get_qual(v);

                if (bcf_float_is_missing(variant_score))
                {
                    variant_score = 0;
                }

                if (vtype==VT_SNP)
                {
                    //GT:PL:DP:AD:ADF:ADR:BQ:MQ:CY:ST:AL:NM
                    if (bcf_get_format_int32(h, v, "GT", &GT, &no_GT) > 0 &&
                        bcf_get_format_int32(h, v, "PL", &PL, &no_PL) > 0 &&
                        bcf_get_format_int32(h, v, "DP", &DP, &no_DP) > 0 &&
                        bcf_get_format_int32(h, v, "ADF", &ADF, &no_ADF) > 0 &&
                        bcf_get_format_int32(h, v, "ADR", &ADR, &no_ADR) > 0)
                    {

                    }
                    //GT:BQSUM:DP
                    else if (bcf_get_format_int32(h, v, "GT", &GT, &no_GT) < 0 ||
                             bcf_get_format_int32(h, v, "BQSUM", &BQSUM, &no_BQSUM) < 0 ||
                             bcf_get_format_int32(h, v, "DP", &DP, &no_DP) < 0)
                    {
                        fprintf(stderr, "[E:%s:%d %s] cannot get format values GT, BQSUM, DP from %s\n", __FILE__, __LINE__, __FUNCTION__, sr->file_names[file_index].c_str());
                        exit(1);
                    }
                    else
                    {
                        fprintf(stderr, "[E:%s:%d %s] cannot get format values GT:PL:DP:ADF:ADR or GT:BQSUM:DP from %s\n", __FILE__, __LINE__, __FUNCTION__, sr->file_names[file_index].c_str());
                        exit(1);

                    }
                }
                else if (vtype == VT_INDEL)
                {

                }
                else if (vtype == VT_VNTR)
                {

                }
                else
                {
                    continue;
                }
            }

//            if (max_variant_score_gt_cutoff)
//            {
//                bcf_update_info_int32(odw->hdr, nv, "NSAMPLES", &no_samples, 1);
//                bcf_update_info_string(odw->hdr, nv, "SAMPLES", samples.c_str());
//                bcf_update_info_int32(odw->hdr, nv, "E", &e[0], no_samples);
//                bcf_update_info_int32(odw->hdr, nv, "N", &n[0], no_samples);
//                bcf_update_info_int32(odw->hdr, nv, "ESUM", &esum, 1);
//                bcf_update_info_int32(odw->hdr, nv, "NSUM", &nsum, 1);
//                bcf_set_qual(nv, max_variant_score);
//
//                odw->write(nv);
//
////                bcf_print(odw->hdr, nv);
//
//                if (vtype == VT_SNP)
//                {
//                    ++no_candidate_snps;
//                }
//                else if (vtype == VT_INDEL)
//                {
//                    ++no_candidate_indels;
//                }
//            }

//            exit(1);
        }

        sr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "merge_candidate_variants v" << version << "\n\n";
        std::clog << "options: [L] input VCF file list         " << input_vcf_file_list << " (" << input_vcf_files.size() << " files)\n";
        std::clog << "         [c] candidate sites VCF file  " << candidate_sites_vcf_file << "\n";
        std::clog << "         [o] output VCF file             " << output_vcf_file << "\n";
        print_int_op("         [i] intervals                   ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: Total Number of Candidate SNPs                 " << no_candidate_snps << "\n";
        std::clog << "       Total Number of Candidate Indels               " << no_candidate_indels << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {
    };

    private:
};

}

void merge_genotypes(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.merge_genotypes();
    igor.print_stats();
}

