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
    int32_t no_chromosomes;
    int32_t no_observed_variants;
    int32_t no_classified_variants;

    int32_t **VAR_COUNT;
    int32_t **VAR_TS;
    int32_t **VAR_TV;
    int32_t **VAR_INS;
    int32_t **VAR_DEL;

    int32_t no_snp3;
    int32_t no_snp4;

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
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            fexp = arg_fexp.getValue();
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

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_samples = 0;
        no_chromosomes = 0;
        no_observed_variants = 0;
        no_classified_variants = 0;

        VAR_COUNT = new int*[6];
        VAR_TS = new int*[6];
        VAR_TV = new int*[6];
        VAR_INS = new int*[6];
        VAR_DEL = new int*[6];

        for (int32_t N_ALLELES=0; N_ALLELES<6; ++N_ALLELES)
        {
            VAR_COUNT[N_ALLELES] = new int[32];
            VAR_TS[N_ALLELES] = new int[32];
            VAR_TV[N_ALLELES] = new int[32];
            VAR_INS[N_ALLELES] = new int[32];
            VAR_DEL[N_ALLELES] = new int[32];
        }

        int32_t vtypes[] = {VT_REF, VT_SNP, VT_MNP, VT_INDEL, VT_SNP|VT_MNP, VT_SNP|VT_INDEL, VT_MNP|VT_INDEL, VT_SNP|VT_MNP|VT_INDEL, VT_CLUMPED};

        for (int32_t i=0; i<9; ++i)
        {
            for (int32_t N_ALLELES=0; N_ALLELES<6; ++N_ALLELES)
            {
                VAR_COUNT[N_ALLELES][vtypes[i]] = 0;
                VAR_TS[N_ALLELES][vtypes[i]] = 0;
                VAR_TV[N_ALLELES][vtypes[i]] = 0;
                VAR_INS[N_ALLELES][vtypes[i]] = 0;
                VAR_DEL[N_ALLELES][vtypes[i]] = 0;
            }
        }

        no_snp3 = 0;
        no_snp4 = 0;

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

        int32_t vtypes[] = {VT_REF, VT_SNP, VT_MNP, VT_INDEL, VT_SNP|VT_MNP, VT_SNP|VT_INDEL, VT_MNP|VT_INDEL, VT_SNP|VT_MNP|VT_INDEL, VT_CLUMPED};

        Variant variant;

        while (odr->read(v))
        {
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);

            if (filter_exists)
            {
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }

            if ((k = kh_get(32, h, bcf_get_rid(v))) == kh_end(h))
            {
                kh_put(32, h, bcf_get_rid(v), &ret);
                kh_value(h, k) = 1; //not really necessary.
                ++no_chromosomes;
            }

            vtype = vtype&VT_CLUMPED ? VT_CLUMPED : vtype;
            int32_t NO_ALLELES = bcf_get_n_allele(v)>5 ? MULTIALLELIC : bcf_get_n_allele(v)-1;

            ++VAR_COUNT[NO_ALLELES][vtype];
            VAR_TS[NO_ALLELES][vtype] += variant.alleles[0].ts;
            VAR_TV[NO_ALLELES][vtype] += variant.alleles[0].tv;
            VAR_INS[NO_ALLELES][vtype] += variant.alleles[0].ins;
            VAR_DEL[NO_ALLELES][vtype] += variant.alleles[0].del;
            
            bool classified = false;
            for (int32_t i=0; i<9; ++i)
            {
                if (vtype == vtypes[i])
                {
                    classified = true;
                    break;
                }
            }

            if (classified) ++no_classified_variants;
            
            
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
        print_str_op("         [f] filter                ", fexp);
        print_ref_op("         [r] reference FASTA file  ", ref_fasta_file);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        int32_t vtypes[] = {VT_REF, VT_SNP, VT_MNP, VT_INDEL, VT_SNP|VT_MNP, VT_SNP|VT_INDEL, VT_MNP|VT_INDEL, VT_SNP|VT_MNP|VT_INDEL, VT_CLUMPED};
        int32_t no_ref = 0;
        
        for (int32_t i=0; i<9; ++i)
        {
            VAR_COUNT[MULTIALLELIC][vtypes[i]] += VAR_COUNT[TRIALLELIC][vtypes[i]] + VAR_COUNT[TETRAALLELIC][vtypes[i]];
            VAR_COUNT[POLYMORPHIC][vtypes[i]] = VAR_COUNT[BIALLELIC][vtypes[i]]  + VAR_COUNT[MULTIALLELIC][vtypes[i]] ;
            VAR_TS[BIALLELIC][vtypes[i]] += VAR_TS[TRIALLELIC][vtypes[i]] + VAR_TS[TETRAALLELIC][vtypes[i]];
            VAR_TV[TRIALLELIC][vtypes[i]] += VAR_TV[TRIALLELIC][vtypes[i]] + VAR_TV[TETRAALLELIC][vtypes[i]];
            VAR_INS[TETRAALLELIC][vtypes[i]] += VAR_INS[TRIALLELIC][vtypes[i]] + VAR_INS[TETRAALLELIC][vtypes[i]];
            VAR_DEL[MULTIALLELIC][vtypes[i]] += VAR_DEL[TRIALLELIC][vtypes[i]] + VAR_DEL[TETRAALLELIC][vtypes[i]];
        }

        fprintf(stderr, "\n");
        fprintf(stderr, "stats: no. of samples                  : %10d\n", no_samples);
        fprintf(stderr, "       no. of chromosomes              : %10d\n", no_chromosomes);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of SNPs                     : %10d\n", VAR_COUNT[POLYMORPHIC][VT_SNP]);
        fprintf(stderr, "           biallelic (ts/tv)           : %15d (%.2f) [%d/%d]\n", VAR_COUNT[BIALLELIC][VT_SNP], (float)VAR_TS[BIALLELIC][VT_SNP]/VAR_TV[BIALLELIC][VT_SNP], VAR_TS[BIALLELIC][VT_SNP], VAR_TV[BIALLELIC][VT_SNP]);
        fprintf(stderr, "           3 alleles                   : %15d\n", VAR_COUNT[TRIALLELIC][VT_SNP]);
        fprintf(stderr, "           4 alleles                   : %15d\n", VAR_COUNT[TETRAALLELIC][VT_SNP]);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of MNPs                     : %10d\n", VAR_COUNT[POLYMORPHIC][VT_MNP]);
        fprintf(stderr, "           biallelic (ts/tv)           : %15d (%.2f) [%d/%d]\n", VAR_COUNT[BIALLELIC][VT_MNP], (float)VAR_TS[BIALLELIC][VT_MNP]/VAR_TV[BIALLELIC][VT_MNP], VAR_TS[BIALLELIC][VT_MNP], VAR_TV[BIALLELIC][VT_MNP]);
        fprintf(stderr, "           multiallelic                : %15d\n", VAR_COUNT[MULTIALLELIC][VT_MNP]);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. Indels                      : %10d\n", VAR_COUNT[POLYMORPHIC][VT_INDEL]);
        fprintf(stderr, "           biallelic (ins/del)         : %15d (%.2f) [%d/%d]\n", VAR_COUNT[BIALLELIC][VT_INDEL], (float)VAR_INS[BIALLELIC][VT_INDEL]/VAR_DEL[BIALLELIC][VT_INDEL], VAR_INS[BIALLELIC][VT_INDEL], VAR_DEL[BIALLELIC][VT_INDEL]);
        fprintf(stderr, "               insertions              : %15d\n", VAR_INS[BIALLELIC][VT_INDEL]);
        fprintf(stderr, "               deletions               : %15d\n", VAR_DEL[BIALLELIC][VT_INDEL]);
        fprintf(stderr, "           multiallelic                : %15d\n", VAR_COUNT[MULTIALLELIC][VT_INDEL]);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. SNP/MNP                     : %10d\n", VAR_COUNT[POLYMORPHIC][VT_SNP|VT_MNP]);
        fprintf(stderr, "           (multiallelic)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. SNP/Indels                  : %10d\n", VAR_COUNT[POLYMORPHIC][VT_SNP|VT_INDEL]);
        fprintf(stderr, "           biallelic (ts/tv) (ins/del) : %15d (%.2f) [%d/%d] (%.2f) [%d/%d]\n",  VAR_COUNT[BIALLELIC][VT_SNP|VT_INDEL], (float)VAR_TS[BIALLELIC][VT_SNP|VT_INDEL]/VAR_TV[BIALLELIC][VT_SNP|VT_INDEL], VAR_TS[BIALLELIC][VT_SNP|VT_INDEL], VAR_TV[BIALLELIC][VT_SNP|VT_INDEL], (float)VAR_INS[BIALLELIC][VT_SNP|VT_INDEL]/VAR_DEL[BIALLELIC][VT_SNP|VT_INDEL], VAR_INS[BIALLELIC][VT_SNP|VT_INDEL], VAR_DEL[BIALLELIC][VT_SNP|VT_INDEL]);
        fprintf(stderr, "               insertions              : %15d\n", VAR_INS[BIALLELIC][VT_SNP|VT_INDEL]);
        fprintf(stderr, "               deletions               : %15d\n", VAR_DEL[BIALLELIC][VT_SNP|VT_INDEL]);
        fprintf(stderr, "           multiallelic                : %15d\n", VAR_COUNT[MULTIALLELIC][VT_SNP|VT_INDEL]);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. MNP/Indels                  : %10d\n", VAR_COUNT[POLYMORPHIC][VT_MNP|VT_INDEL]);
        fprintf(stderr, "           biallelic (ts/tv) (ins/del) : %15d (%.2f) [%d/%d] (%.2f) [%d/%d]\n",  VAR_COUNT[BIALLELIC][VT_MNP|VT_INDEL], (float)VAR_TS[BIALLELIC][VT_MNP|VT_INDEL]/VAR_TV[BIALLELIC][VT_MNP|VT_INDEL], VAR_TS[BIALLELIC][VT_MNP|VT_INDEL], VAR_TV[BIALLELIC][VT_MNP|VT_INDEL], (float)VAR_INS[BIALLELIC][VT_MNP|VT_INDEL]/VAR_DEL[BIALLELIC][VT_MNP|VT_INDEL], VAR_INS[BIALLELIC][VT_MNP|VT_INDEL], VAR_DEL[BIALLELIC][VT_MNP|VT_INDEL]);
        fprintf(stderr, "               insertions              : %15d\n", VAR_INS[BIALLELIC][VT_MNP|VT_INDEL]);
        fprintf(stderr, "               deletions               : %15d\n", VAR_DEL[BIALLELIC][VT_MNP|VT_INDEL]);
        fprintf(stderr, "           multiallelic                : %15d\n", VAR_COUNT[MULTIALLELIC][VT_MNP|VT_INDEL]);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. SNP/MNP/Indels              : %10d\n", VAR_COUNT[MULTIALLELIC][VT_SNP|VT_MNP|VT_INDEL]);
        fprintf(stderr, "           (multiallelic)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of clumped variants         : %10d\n", VAR_COUNT[POLYMORPHIC][VT_CLUMPED]);
        fprintf(stderr, "           biallelic                   : %15d\n", VAR_COUNT[BIALLELIC][VT_CLUMPED]);
        fprintf(stderr, "           multiallelic                : %15d\n", VAR_COUNT[MULTIALLELIC][VT_CLUMPED]);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of reference                : %10d\n", VAR_COUNT[MONOMORPHIC][VT_REF]);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of observed variants      : %10d\n", no_observed_variants);
        fprintf(stderr, "       no. of unclassified variants  : %10d\n", no_observed_variants-no_classified_variants);
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
