/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#include "peek.h"

namespace
{

KHASH_MAP_INIT_INT(32, char)

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
    bcf1_t *v;

    /////////
    //stats//
    /////////
    uint32_t no_samples;
    uint32_t no_chromosomes;
    uint32_t no_observed_variants;
    uint32_t no_classified_variants;
    uint32_t no_ref;
    uint32_t no_snp2;
    uint32_t no_snp3;
    uint32_t no_snp4;
    uint32_t no_mnp2;
    uint32_t no_mnp_multi;
    uint32_t no_indel2;
    uint32_t no_ins2;
    uint32_t no_del2;
    uint32_t no_indel_multi;
    uint32_t no_snpindel2;
    uint32_t no_snpins2;
    uint32_t no_snpdel2;
    uint32_t no_snpindel_multi;
    uint32_t no_mnpindel2;
    uint32_t no_mnpins2;
    uint32_t no_mnpdel2;
    uint32_t no_mnpindel_multi;
    uint32_t no_clumped2;
    uint32_t no_clumped_multi;
    
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
            std::string desc = "Summarizes the variants in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_fasta_file = arg_ref_fasta_file.getValue();
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
		v = bcf_init1();

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_samples = 0;
        no_chromosomes = 0;
        no_observed_variants = 0;
        no_classified_variants = 0;
        no_ref = 0;
        no_snp2 = 0;
        no_snp3 = 0;
        no_snp4 = 0;
        no_mnp2 = 0;
        no_mnp_multi = 0;
        no_indel2 = 0;
        no_ins2 = 0;
        no_del2 = 0;
        no_indel_multi = 0;
        no_snpindel2 = 0;
        no_snpins2 = 0;
        no_snpdel2 = 0;
        no_snpindel_multi = 0;
        no_mnpindel2 = 0;
        no_mnpins2 = 0;
        no_mnpdel2 = 0;
        no_mnpindel_multi = 0;
        no_clumped2 = 0;
        no_clumped_multi = 0;
        no_clumped2 = 0;
        no_clumped_multi = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip(ref_fasta_file);
    }

    void peek()
    {
        no_samples = bcf_hdr_get_n_sample(odr->hdr);

        int ret, is_missing;
    	khiter_t k;
    	khash_t(32) *h = kh_init(32);

        Variant variant;

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            
            if ((k = kh_get(32, h, bcf_get_rid(v))) == kh_end(h))
            {
                kh_put(32, h, bcf_get_rid(v), &ret);
                kh_value(h, k) = 1; //not really necessary.
                ++no_chromosomes;   
            }

            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            
            if (vtype & VT_CLUMPED)
            {
                if (bcf_get_n_allele(v)==2)
                {
                    ++no_clumped2;
                }
                else if (bcf_get_n_allele(v)>=3)
                {
                    ++no_clumped_multi;
                }
                
                ++no_classified_variants;
            }
            else if (vtype==VT_SNP)
            {
                if (bcf_get_n_allele(v)==2)
                {
                    ++no_snp2;
                }
                else if (bcf_get_n_allele(v)==3)
                {
                    ++no_snp3;
                }
                else if (bcf_get_n_allele(v)==4)
                {
                    ++no_snp4;
                }
                
                ++no_classified_variants;
            }
            else if (vtype==VT_MNP)
            {
                if (bcf_get_n_allele(v)==2)
                {
                    ++no_mnp2;
                }
                else if (bcf_get_n_allele(v)>=3)
                {
                    ++no_mnp_multi;
                }
                
                ++no_classified_variants;
            }
            else if (vtype==VT_INDEL) //strictly simple indels
            {
                if (bcf_get_n_allele(v)==2)
                {
                    ++no_indel2;
                    if (variant.alleles[0].dlen>0)
                    {
                        ++no_ins2;
                    }
                    else
                    {
                        ++no_del2;
                    }
                }
                else if (bcf_get_n_allele(v)>=3)
                {      
                    ++no_indel_multi;
                }
                
                ++no_classified_variants;
            }
            else if (vtype==(VT_SNP|VT_INDEL))
            {
                if (bcf_get_n_allele(v)==2)
                {
                    ++no_snpindel2;
                    if (variant.alleles[0].dlen>0)
                    {
                        ++no_snpins2;
                    }
                    else
                    {
                        ++no_snpdel2;
                    }
                }
                else if (bcf_get_n_allele(v)>=3)
                {
                    ++no_snpindel_multi;
                }
                                
                ++no_classified_variants;
            }
            else if (vtype&(VT_MNP|VT_INDEL))
            {
                if (bcf_get_n_allele(v)==2)
                {
                    ++no_mnpindel2;
                    if (variant.alleles[0].dlen>0)
                    {
                        ++no_mnpins2;
                    }
                    else
                    {
                        ++no_mnpdel2;
                    }
                }
                else if (bcf_get_n_allele(v)>=3)
                {
                    ++no_mnpindel_multi;
                }
             
                ++no_classified_variants;
            }
            else if (vtype==VT_REF) //MNPs that are not real MNPs
            {
                ++no_ref;
                ++no_classified_variants;
            }
            else
            {
                std::cerr << "UNCLASSIFIED\n";
                bcf_print_lite(odr->hdr, v);
                std::cerr << vm->vtype2string(vtype) << "\n"; 
            }
            
            ++no_observed_variants;
        }

        kh_destroy(32, h);
        odr->close();
    };

    void print_options()
    {
        std::clog << "peek v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_ref_op("         [r] reference FASTA file  ", ref_fasta_file);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "stats: No. of samples                : %10d\n", no_samples);
        fprintf(stderr, "       No. of chromosomes            : %10d\n", no_chromosomes);
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. of SNPs                   : %10d\n", no_snp2+no_snp3+no_snp4);
        fprintf(stderr, "           biallelic                 : %15d\n", no_snp2);
        fprintf(stderr, "           3 alleles                 : %15d\n", no_snp3);
        fprintf(stderr, "           4 alleles                 : %15d\n", no_snp4);
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. of MNPs                   : %10d\n", no_mnp2+no_mnp_multi);
        fprintf(stderr, "           biallelic                 : %15d\n", no_mnp2);
        fprintf(stderr, "           multiallelic              : %15d\n", no_mnp_multi);
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. Indels                    : %10d\n", no_indel2+no_indel_multi);
        fprintf(stderr, "           biallelic                 : %15d\n", no_indel2);
        fprintf(stderr, "               insertions            : %15d\n", no_ins2);
        fprintf(stderr, "               deletions             : %15d\n", no_del2);
        fprintf(stderr, "           multiallelic              : %15d\n", no_indel_multi);
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. SNP/Indels                : %10d\n", no_snpindel2+no_snpindel_multi);
        fprintf(stderr, "           biallelic                 : %15d\n", no_snpindel2);
        fprintf(stderr, "               insertions            : %15d\n", no_snpins2);
        fprintf(stderr, "               deletions             : %15d\n", no_snpdel2);
        fprintf(stderr, "           multiallelic              : %15d\n", no_snpindel_multi);   
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. MNP/Indels                : %10d\n", no_mnpindel2+no_mnpindel_multi);
        fprintf(stderr, "           biallelic                 : %15d\n", no_mnpindel2);
        fprintf(stderr, "               insertions            : %15d\n", no_mnpins2);
        fprintf(stderr, "               deletions             : %15d\n", no_mnpdel2);
        fprintf(stderr, "           multiallelic              : %15d\n", no_mnpindel_multi);   
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. of clumped variants       : %10d\n", no_clumped2+no_clumped_multi);
        fprintf(stderr, "           biallelic                 : %15d\n", no_clumped2);
        fprintf(stderr, "           multiallelic              : %15d\n", no_clumped_multi);
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. of reference              : %10d\n", no_ref);
        fprintf(stderr, "\n");
        fprintf(stderr, "       No. of observed variants      : %10d\n", no_observed_variants);
        fprintf(stderr, "       No. of unclassified variants  : %10d\n", no_observed_variants-no_classified_variants);
        fprintf(stderr, "\n");
    };

    ~Igor() {};

    private:
};

}

void peek(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.peek();
    igor.print_stats();
};
