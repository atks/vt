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

#include "compute_features2.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    std::string arg_sample_list;
    char** samples;
    int32_t *imap;
    int32_t nsamples;
    bool print_sites_only;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    int32_t no_samples;
    int32_t no_variants;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    Estimator *est;

    Igor(int argc, char ** argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Compute features for variants in vt pipeline from separate genotype files";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "intervals", "Intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "interval-list", "File containing list of intervals", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF/VCF.GZ/BCF file [-]", false, "-", "str", cmd);
            TCLAP::SwitchArg arg_print_sites_only("s", "s", "print site information only without genotypes [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "", "file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print_sites_only = arg_print_sites_only.getValue();
            fexp = arg_fexp.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        //convert to paste like set up
        
        
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        if (print_sites_only)
        {
            odw->link_hdr(bcf_hdr_subset(odr->hdr, 0, 0, 0));
        }
        else
        {
            odw->link_hdr(odr->hdr);
        }

        bcf_hdr_append(odw->hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=GN,Number=1,Type=Integer,Description=\"Total Number of Genotypes Counts\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=GF,Number=G,Type=Float,Description=\"Genotype Frequency\">\n");
        
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWEAF,Number=A,Type=Float,Description=\"Genotype likelihood based MLE Allele Frequency assuming HWE\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWEGF,Number=G,Type=Float,Description=\"Genotype likelihood based MLE Genotype Frequency assuming HWE\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Genotype likelihood based MLE Allele Frequency\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=MLEGF,Number=G,Type=Float,Description=\"Genotype likelihood based MLE Genotype Frequency\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWE_LLR,Number=1,Type=Float,Description=\"Genotype likelihood based Hardy Weinberg ln(Likelihood Ratio)\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWE_LPVAL,Number=1,Type=Float,Description=\"Genotype likelihood based Hardy Weinberg Likelihood Ratio Test Statistic ln(p-value)\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWE_DF,Number=1,Type=Integer,Description=\"Degrees of freedom for Genotype likelihood based Hardy Weinberg Likelihood Ratio Test Statistic\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=FIC,Number=1,Type=Float,Description=\"Genotype likelihood based Inbreeding Coefficient\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Genotype likelihood based Allele Balance\">\n");
        
        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_samples = bcf_hdr_nsamples(odr->hdr);
        no_variants = 0;

        if (!no_samples)
        {
            fprintf(stderr, "[%s:%d %s] No samples in VCF file: %s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
            exit(1);
        }
        
        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");
        est = new Estimator();
    }

    void compute_features()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odr->hdr;
        Variant variant;

        int32_t *gts = NULL;
        int32_t *pls = NULL;
        int32_t *dps = NULL;
        int32_t n_gts = 0;
        int32_t n_pls = 0;
        int32_t n_dps = 0;

        odw->write_hdr();

        while(odr->read(v))
        {
            variant.clear();
            bool printed = false;

            vm->classify_variant(h, v, variant);
            if (filter_exists && !filter.apply(h,v,&variant))
            {
                continue;
            }

            bcf_unpack(v, BCF_UN_ALL);
            int32_t ploidy = bcf_get_genotypes(odr->hdr, v, &gts, &n_gts);
            ploidy /= no_samples;
            
            if (!n_gts)
            {
                continue;
            }
        
            bcf_get_format_int32(odr->hdr, v, "PL", &pls, &n_pls);
            int32_t no_alleles = bcf_get_n_allele(v);

            if (!n_pls)
            {
                continue;
            }

            float qual = 0;
            int32_t n = 0;
            est->compute_qual(pls, no_samples, ploidy, no_alleles, qual, n);
            if (n)
            {
                bcf_set_qual(v, qual);
            }
            
            int32_t g[ploidy];
            for (int32_t i=0; i<ploidy; ++i) g[i]=0;
            int32_t AC[no_alleles];
            float AF[no_alleles];
            for (int32_t i=0; i<no_alleles; ++i) AC[i]=0;
            int32_t AN=0;
            int32_t NS=0;
            int32_t no_genotypes = bcf_an2gn(no_alleles);
            int32_t GC[no_genotypes];
            int32_t GN=0;
            float GF[no_genotypes];
            
            est->compute_af(gts, no_samples, ploidy, no_alleles, AC, AN, AF, GC, GN, GF, NS);

            if (NS)
            {
                int32_t* AC_PTR = &AC[1];
                bcf_update_info_int32(odw->hdr, v, "AC", AC_PTR, no_alleles-1);
                bcf_update_info_int32(odw->hdr, v, "AN", &AN, 1);
                float* AF_PTR = &AF[1];
                bcf_update_info_float(odw->hdr, v, "AF", AF_PTR, no_alleles-1);
                if (GN)            
                {
                    bcf_update_info_int32(odw->hdr, v, "GC", GC, no_genotypes);
                    bcf_update_info_int32(odw->hdr, v, "GN", &GN, 1);
                    bcf_update_info_float(odw->hdr, v, "GF", GF, no_genotypes);
                }
                bcf_update_info_int32(odw->hdr, v, "NS", &NS, 1);
            }
 
            float MLE_HWE_AF[no_alleles];
            float MLE_HWE_GF[no_genotypes];
            n = 0;
            est->compute_gl_af_hwe(pls, no_samples, ploidy,no_alleles, MLE_HWE_AF, MLE_HWE_GF,  n, 1e-20);
            if (n)
            {
                float* MLE_HWE_AF_PTR = &MLE_HWE_AF[1];
                bcf_update_info_float(odw->hdr, v, "HWEAF", MLE_HWE_AF_PTR, no_alleles-1);
                bcf_update_info_float(odw->hdr, v, "HWEGF", &MLE_HWE_GF, no_genotypes);
            }

            float MLE_AF[no_alleles];
            float MLE_GF[no_genotypes];
            n = 0;
            est->compute_gl_af(pls, no_samples, ploidy,no_alleles, MLE_AF, MLE_GF,  n, 1e-20);
            if (n)
            {
                float* MLE_AF_PTR = &MLE_AF[1];
                bcf_update_info_float(odw->hdr, v, "MLEAF", MLE_AF_PTR, no_alleles-1);
                bcf_update_info_float(odw->hdr, v, "MLEGF", &MLE_GF, no_genotypes);
            }

            float lrts;
            float logp;
            int32_t df;
            n = 0;
            est->compute_hwe_lrt(pls, no_samples, ploidy,
                                 no_alleles, MLE_HWE_GF, MLE_GF, n,
                                 lrts, logp, df);
            if (n)
            {
                bcf_update_info_float(odw->hdr, v, "HWE_LLR", &lrts, 1);
                bcf_update_info_float(odw->hdr, v, "HWE_LPVAL", &logp, 1);
                bcf_update_info_int32(odw->hdr, v, "HWE_DF", &df, 1);
            }

            float f;
            n = 0;
            est->compute_gl_fic(pls, no_samples, ploidy,
                               MLE_HWE_AF, no_alleles, MLE_GF,
                               f, n);
            if (n)
            {
                bcf_update_info_float(odw->hdr, v, "FIC", &f, 1);
            }

            bcf_get_format_int32(odr->hdr, v, "DP", &dps, &n_dps);
            if (n_dps)
            {
                float ab;
                n = 0;
                est->compute_gl_ab(pls, no_samples, ploidy,
                                   dps, MLE_GF, no_alleles, 
                                   ab, n);
                                   
                if (n)
                {
                    bcf_update_info_float(odw->hdr, v, "AB", &ab, 1);
                }
            }
            
            if (print_sites_only)
            {
                bcf_subset(odw->hdr, v, 0, 0);
            }

            odw->write(v);
            ++no_variants;
        }

        if(n_gts) free(gts);
        if(n_pls) free(pls);
        if(n_dps) free(dps);
            
        odw->close();
    };

    void print_options()
    {
        std::clog << "compute_features v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File    " << input_vcf_file << "\n";
        print_str_op("         [f] filter            ", fexp);
        print_int_op("         [i] Intervals         ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: variants   : " << no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {
    };

    private:
};

}

void compute_features2(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.compute_features();
    igor.print_stats();
}

