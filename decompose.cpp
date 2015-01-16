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

#include "decompose.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string ref_fasta_file;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    kstring_t s;
    kstring_t new_alleles;
    kstring_t old_alleles;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    uint32_t new_no_variants;
    uint32_t no_biallelic;
    uint32_t no_multiallelic;
    uint32_t no_additional_biallelic;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "decomposes multialleic variants into biallelic in a VCF file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
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
        odr = new BCFOrderedReader(input_vcf_file, intervals);

        odw = new BCFOrderedWriter(output_vcf_file, 1000);
        odw->set_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description=\"Original chr:pos:ref:alt encoding\">\n");
        odw->write_hdr();

        s = {0,0,0};
        old_alleles = {0,0,0};
        new_alleles = {0,0,0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        new_no_variants = 0;
        no_biallelic = 0;
        no_multiallelic = 0;
        no_additional_biallelic = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void decompose()
    {
        v = odw->get_bcf1_from_pool();
        Variant variant;

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            int32_t n_allele = bcf_get_n_allele(v);

            if (n_allele > 2)
            {
                ++no_multiallelic;
                no_additional_biallelic += n_allele-1;

                old_alleles.l = 0;
                bcf_variant2string(odw->hdr, v, &old_alleles);

                int32_t rid = bcf_get_rid(v);
                int32_t pos1 = bcf_get_pos1(v);
                char** allele = bcf_get_allele(v);

                char** alleles = (char**) malloc(n_allele*sizeof(char*));
                for (size_t i=0; i<n_allele; ++i)
                {
                    alleles[i] = strdup(allele[i]);
                }


                int32_t *gt = NULL;
                int32_t *pl = NULL;
                float *gl = NULL;
                int32_t n = 0;
                size_t no_samples = bcf_hdr_nsamples(odr->hdr);
                bool has_GT = false;
                bool has_PL = false;
                bool has_GL = false;
                int32_t ploidy = 0;
                int32_t n_genotype;
                
                int32_t *gts;
                int32_t *pls;
                float *gls;

                if (no_samples)
                {
                    int32_t ret = bcf_get_genotypes(odr->hdr, v, &gt, &n);
                    if (ret>0) has_GT = true;
                    ploidy = n/bcf_hdr_nsamples(odr->hdr);
                    n_genotype = bcf_ap2g(n_allele, ploidy);
                    if (ret>0) gts = (int32_t*) malloc(no_samples*ploidy*sizeof(int32_t));

                    ret = bcf_get_format_int32(odr->hdr, v, "PL", &pl, &n);
                    if (ret>0) has_PL = true;
                    if (ret>0) pls = (int32_t*) malloc(no_samples*n_genotype*sizeof(int32_t));

                    ret = bcf_get_format_float(odr->hdr, v, "GL", &gl, &n);
                    if (ret>0) has_GL = true;   
                    if (ret>0) gls = (float*) malloc(no_samples*n_genotype*sizeof(float));

//                    std::cerr << "return value    : " << ret << "\n";
//                    std::cerr << "n bytes written : " << n << "\n";
//                    std::cerr << "no samples      : " << bcf_hdr_nsamples(odr->hdr) << "\n";
//                    std::cerr << "ploidy          : " << ploidy << "\n";
//                    std::cerr << "GT              : " << has_GT << "\n";
//                    std::cerr << "PL              : " << has_PL << "\n";
//                    std::cerr << "GL              : " << has_GL << "\n";
                }

                for (size_t i=1; i<n_allele; ++i)
                {
                    bcf_set_rid(v, rid);
                    bcf_set_pos1(v, pos1);
                    new_alleles.l=0;
                    kputs(alleles[0], &new_alleles);
                    kputc(',', &new_alleles);
                    kputs(alleles[i], &new_alleles);

                    bcf_update_alleles_str(odw->hdr, v, new_alleles.s);
                    bcf_update_info_string(odw->hdr, v, "OLD_MULTIALLELIC", old_alleles.s);

                    if (no_samples)
                    {
                        //get array genotypes
                        for (size_t j=0; j<no_samples; ++j)
                        {
                            if (has_GT)
                            {
                                for (size_t k=0; k<ploidy; ++k)
                                {
                                    int32_t g = bcf_gt_allele(gt[j*ploidy+k]);
                                    
                                    if (g<=0)
                                    {
                                    }
                                    else if (g!=i)
                                    {
                                        g = -1;
                                    }
                                    else
                                    {
                                        g = 1;
                                    }
                                    
                                    gts[j*ploidy+k] = bcf_gt_unphased(g);                                    
                                }
                            }

                            if (has_PL)
                            {
                                if (pl[j*n_genotype]!=bcf_int8_missing &&
                                    pl[j*n_genotype]!=bcf_int16_missing &&
                                    pl[j*n_genotype]!=bcf_int32_missing)
                                {
                                    int32_t plref = pl[j*n_genotype];
                                    int32_t plhet = pl[j*n_genotype+bcf_alleles2gt(0,i)];
                                    int32_t plalt = pl[j*n_genotype+bcf_alleles2gt(i,i)];

                                    for (size_t k=0; k<n_genotype; ++k)
                                    {
                                        pls[j*n_genotype+k] = plref;
                                    } 
                                }
                                else
                                {
                                    for (size_t k=0; k<n_genotype; ++k)
                                    {
                                        pls[j*n_genotype+k] = bcf_int32_missing;
                                    }
                                }
                            }

                            if (has_GL)
                            {
                                if (!bcf_float_is_missing(gl[j*n_genotype]))
                                {
                                    for (size_t k=0; k<n_genotype; ++k)
                                    {
                                        gls[j*n_genotype+k] = gl[j*n_genotype+k];
                                    }    

                                }
                                else
                                {
                                    for (size_t k=0; k<n_genotype; ++k)
                                    {
                                        bcf_float_set_missing(gls[j*n_genotype+k]);
                                    }  
                                }
                            }
                        }

                        //remove other format values except for GT, PL and GL
                        if (i==1)
                        {
                            bcf_fmt_t *fmt = v->d.fmt;
                            for (size_t j = 0; j < v->n_fmt; ++j) 
                            {
                                const char* tag = odw->hdr->id[BCF_DT_ID][fmt[j].id].key;
                               
                                if (strcmp(tag,"GT")&&strcmp(tag,"PL")&&strcmp(tag,"GL"))
                                {
                                    bcf_update_format_int32(odw->hdr, v, tag, 0, 0);
                                }
                            }
                        }    

                        if (has_GT) bcf_update_genotypes(odw->hdr, v, gts, no_samples*ploidy);
                        if (has_PL) bcf_update_format_int32(odw->hdr, v, "PL", pls, no_samples*n_genotype);
                        if (has_GL) bcf_update_format_float(odw->hdr, v, "GL", gls, no_samples*n_genotype);
                    }

                    odw->write(v);
                    v = odw->get_bcf1_from_pool();
                    ++new_no_variants;
                }

                for (size_t i=0; i<n_allele; ++i)
                {
                    free(alleles[i]);
                }
                free(alleles);
            }
            else
            {
                ++no_biallelic;

                odw->write(v);
                v = odw->get_bcf1_from_pool();
                ++new_no_variants;
            }

            ++no_variants;
        }

        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "decompose v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. variants                 : " << no_variants << "\n";
        std::clog << "       no. biallelic variants       : " << no_biallelic << "\n";
        std::clog << "       no. multiallelic variants    : " << no_multiallelic << "\n";
        std::clog << "\n";
        std::clog << "       no. additional biallelics    : " << no_additional_biallelic << "\n";
        std::clog << "       new no. variants             : " << new_no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void decompose(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.decompose();
    igor.print_stats();
};
