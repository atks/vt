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

#define MONOMORPHIC     0
#define BIALLELIC       1
#define TRIALLELIC      2
#define TETRAALLELIC    3
#define GE_PENTAALLELIC 4
#define GE_TRIALLELIC   5
#define GE_TETRAALLELIC 6
#define POLYMORPHIC     7

#define NO_VARIANT_CATEGORIES       40
#define NO_BASIC_VARIANT_CATEGORIES 16
#define NO_ALLELE_CATEGORIES         8
#define NO_MOTIF_LEN_CATEGORIES     10

#define VT_NAIVE_CLUMPED   17
#define VT_BLKSUB          18
#define VT_CPLXSUB         19

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
    std::string output_tabulate_dir;
    std::string output_pdf_file;
    std::vector<GenomeInterval> intervals;
    std::string ref_fasta_file;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    bcf1_t *v;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    int32_t no_samples;
    int32_t no_chromosomes;
    int32_t no_records;
    int32_t no_reference;
    int32_t no_classified_variants;
    int32_t no_unclassified_variants;

    int32_t VAR_COUNT[NO_ALLELE_CATEGORIES][NO_VARIANT_CATEGORIES];
    int32_t VAR_TS[NO_ALLELE_CATEGORIES][NO_VARIANT_CATEGORIES];
    int32_t VAR_TV[NO_ALLELE_CATEGORIES][NO_VARIANT_CATEGORIES];
    int32_t VAR_INS[NO_ALLELE_CATEGORIES][NO_VARIANT_CATEGORIES];
    int32_t VAR_DEL[NO_ALLELE_CATEGORIES][NO_VARIANT_CATEGORIES];
    int32_t VAR_MOTIF_LEN[NO_MOTIF_LEN_CATEGORIES];

    int32_t no_snp3;
    int32_t no_snp4;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    SVTree* sv;

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
            TCLAP::ValueArg<std::string> arg_output_tabulate_dir("x", "x", "output latex directory []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output pdf file [summary.pdf]", false, "", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            fexp = arg_fexp.getValue();
            output_tabulate_dir = arg_output_tabulate_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
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
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_samples = 0;
        no_chromosomes = 0;
        no_records = 0;
        no_reference = 0;
        no_classified_variants = 0;
        no_unclassified_variants = 0;

        for (int32_t vtype=0; vtype<NO_VARIANT_CATEGORIES; ++vtype)
        {
            for (int32_t no_alleles=0; no_alleles<NO_ALLELE_CATEGORIES; ++no_alleles)
            {
                VAR_COUNT[no_alleles][vtype] = 0;
                VAR_TS[no_alleles][vtype] = 0;
                VAR_TV[no_alleles][vtype] = 0;
                VAR_INS[no_alleles][vtype] = 0;
                VAR_DEL[no_alleles][vtype] = 0;
            }
        }

        for (int32_t i=0; i<NO_MOTIF_LEN_CATEGORIES; ++i)
        {
            VAR_MOTIF_LEN[i] = 0;
        }

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip(ref_fasta_file);
        sv = new SVTree();
    }

    void peek()
    {
        no_samples = bcf_hdr_get_n_sample(odr->hdr);

        int ret, is_missing;
        khiter_t k;
        khash_t(32) *h = kh_init(32);

        int32_t vtypes[] = {VT_REF, VT_SNP, VT_MNP, VT_INDEL, VT_SNP|VT_MNP, VT_SNP|VT_INDEL, VT_MNP|VT_INDEL, VT_SNP|VT_MNP|VT_INDEL, VT_CLUMPED, VT_VNTR, VT_SV};

        Variant variant;

        while (odr->read(v))
        {
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);

//            bcf_print(odr->hdr, v);

            if (filter_exists)
            {
                if (!filter.apply(odr->hdr, v, &variant, false))
                {
                    continue;
                }
            }

            //observing chromosomes
            if ((k = kh_get(32, h, bcf_get_rid(v)))==kh_end(h))
            {
                kh_put(32, h, bcf_get_rid(v), &ret);
                kh_value(h, k) = 1; //not really necessary.
                ++no_chromosomes;
            }

            int32_t no_alleles = bcf_get_n_allele(v)-1>=GE_PENTAALLELIC ? GE_PENTAALLELIC : bcf_get_n_allele(v)-1;
            ++VAR_COUNT[no_alleles][vtype];

            for (size_t i=0; i<no_alleles; ++i)
            {
                VAR_TS[no_alleles][vtype] += variant.alleles[i].ts;
                VAR_TV[no_alleles][vtype] += variant.alleles[i].tv;
                VAR_INS[no_alleles][vtype] += variant.alleles[i].ins;
                VAR_DEL[no_alleles][vtype] += variant.alleles[i].del;
            }

            if (vtype==VT_VNTR)
            {
                ++VAR_COUNT[POLYMORPHIC][VT_VNTR];
                if (variant.vntr.motif.size()<NO_MOTIF_LEN_CATEGORIES)
                {
                    ++VAR_MOTIF_LEN[variant.vntr.motif.size()-1];
                }
                else
                {
                    ++VAR_MOTIF_LEN[NO_MOTIF_LEN_CATEGORIES-1];
                }
            }

            if (vtype==VT_SV)
            {
                ++VAR_COUNT[POLYMORPHIC][VT_SV];
                sv->count(variant);
            }

            if (vtype>0 && vtype<64)
            {
                ++no_classified_variants;
            }
            else if (vtype==0)
            {
                ++no_reference;
            }
            else
            {
                std::cerr << "UNCLASSIFIED : ";
                bcf_print(odr->hdr, v);
                ++no_unclassified_variants;
            }

            ++no_records;
        }

        kh_destroy(32, h);
        odr->close();
    };

    void print_options()
    {
        std::clog << "peek v" << version << "\n\n";

        std::clog << "options:     input VCF file            " << input_vcf_file << "\n";
        print_str_op("         [f] filter                    ", fexp);
        print_ref_op("         [r] reference FASTA file      ", ref_fasta_file);
        print_str_op("         [x] output tabulate directory ", output_tabulate_dir);
        print_str_op("         [y] output pdf file           ", output_pdf_file);
        print_int_op("         [i] intervals                 ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        for (size_t vtype=0; vtype<NO_BASIC_VARIANT_CATEGORIES; ++vtype)
        {
            if (vtype&VT_CLUMPED)
            {
                for (size_t no_alleles=0; no_alleles<NO_ALLELE_CATEGORIES; ++no_alleles)
                {
                    VAR_COUNT[no_alleles][VT_NAIVE_CLUMPED] += VAR_COUNT[no_alleles][vtype];
                    VAR_TS[no_alleles][VT_NAIVE_CLUMPED] += VAR_TS[no_alleles][vtype];
                    VAR_TV[no_alleles][VT_NAIVE_CLUMPED] += VAR_TV[no_alleles][vtype];
                    VAR_INS[no_alleles][VT_NAIVE_CLUMPED] += VAR_INS[no_alleles][vtype];
                    VAR_DEL[no_alleles][VT_NAIVE_CLUMPED] += VAR_DEL[no_alleles][vtype];
                }
            }
        }

        int32_t blksub_vtype[] = {VT_MNP, VT_SNP|VT_MNP, VT_SNP|VT_MNP|VT_CLUMPED, VT_MNP|VT_CLUMPED};
        for (size_t i=0; i<4; ++i)
        {
            int32_t vtype = blksub_vtype[i];

            for (size_t no_alleles=0; no_alleles<NO_ALLELE_CATEGORIES; ++no_alleles)
            {
                VAR_COUNT[no_alleles][VT_BLKSUB] += VAR_COUNT[no_alleles][vtype];
                VAR_TS[no_alleles][VT_BLKSUB] += VAR_TS[no_alleles][vtype];
                VAR_TV[no_alleles][VT_BLKSUB] += VAR_TV[no_alleles][vtype];
                VAR_INS[no_alleles][VT_BLKSUB] += VAR_INS[no_alleles][vtype];
                VAR_DEL[no_alleles][VT_BLKSUB] += VAR_DEL[no_alleles][vtype];
            }
        }

        int32_t cplxsub_vtype[] = {VT_SNP|VT_INDEL, VT_MNP|VT_INDEL, VT_SNP|VT_MNP|VT_INDEL,
                                   VT_CLUMPED, VT_INDEL|VT_CLUMPED, VT_SNP|VT_INDEL|VT_CLUMPED,
                                   VT_MNP|VT_INDEL|VT_CLUMPED, VT_SNP|VT_MNP|VT_INDEL|VT_CLUMPED};
        for (size_t i=0; i<(sizeof(cplxsub_vtype)/sizeof(int32_t)); ++i)
        {
            int32_t vtype = cplxsub_vtype[i];

            for (size_t no_alleles=0; no_alleles<NO_ALLELE_CATEGORIES; ++no_alleles)
            {
                VAR_COUNT[no_alleles][VT_CPLXSUB] += VAR_COUNT[no_alleles][vtype];
                VAR_TS[no_alleles][VT_CPLXSUB] += VAR_TS[no_alleles][vtype];
                VAR_TV[no_alleles][VT_CPLXSUB] += VAR_TV[no_alleles][vtype];
                VAR_INS[no_alleles][VT_CPLXSUB] += VAR_INS[no_alleles][vtype];
                VAR_DEL[no_alleles][VT_CPLXSUB] += VAR_DEL[no_alleles][vtype];
            }
        }

        for (int32_t vtype=0; vtype<NO_VARIANT_CATEGORIES; ++vtype)
        {
            VAR_COUNT[GE_TETRAALLELIC][vtype] = VAR_COUNT[TETRAALLELIC][vtype] + VAR_COUNT[GE_PENTAALLELIC][vtype];
            VAR_COUNT[GE_TRIALLELIC][vtype] = VAR_COUNT[TRIALLELIC][vtype] + VAR_COUNT[GE_TETRAALLELIC][vtype];
            VAR_COUNT[POLYMORPHIC][vtype] = VAR_COUNT[BIALLELIC][vtype] + VAR_COUNT[GE_TRIALLELIC][vtype];

            VAR_TS[GE_TETRAALLELIC][vtype] = VAR_TS[TETRAALLELIC][vtype] + VAR_TS[GE_PENTAALLELIC][vtype];
            VAR_TS[GE_TRIALLELIC][vtype] = VAR_TS[TRIALLELIC][vtype] + VAR_TS[GE_TETRAALLELIC][vtype];
            VAR_TS[POLYMORPHIC][vtype] = VAR_TS[BIALLELIC][vtype] + VAR_TS[GE_TRIALLELIC][vtype];

            VAR_TV[GE_TETRAALLELIC][vtype] = VAR_TV[TETRAALLELIC][vtype] + VAR_TV[GE_PENTAALLELIC][vtype];
            VAR_TV[GE_TRIALLELIC][vtype] = VAR_TV[TRIALLELIC][vtype] + VAR_TV[GE_TETRAALLELIC][vtype];
            VAR_TV[POLYMORPHIC][vtype] = VAR_TV[BIALLELIC][vtype] + VAR_TV[GE_TRIALLELIC][vtype];

            VAR_INS[GE_TETRAALLELIC][vtype] = VAR_INS[TETRAALLELIC][vtype] + VAR_INS[GE_PENTAALLELIC][vtype];
            VAR_INS[GE_TRIALLELIC][vtype] = VAR_INS[TRIALLELIC][vtype] + VAR_INS[GE_TETRAALLELIC][vtype];
            VAR_INS[POLYMORPHIC][vtype] = VAR_INS[BIALLELIC][vtype] + VAR_INS[GE_TRIALLELIC][vtype];

            VAR_DEL[GE_TETRAALLELIC][vtype] = VAR_DEL[TETRAALLELIC][vtype] + VAR_DEL[GE_PENTAALLELIC][vtype];
            VAR_DEL[GE_TRIALLELIC][vtype] = VAR_DEL[TRIALLELIC][vtype] + VAR_DEL[GE_TETRAALLELIC][vtype];
            VAR_DEL[POLYMORPHIC][vtype] = VAR_DEL[BIALLELIC][vtype] + VAR_DEL[GE_TRIALLELIC][vtype];
        }

        fprintf(stderr, "\n");
        fprintf(stderr, "stats: no. of samples                     : %10d\n", no_samples);
        fprintf(stderr, "       no. of chromosomes                 : %10d\n", no_chromosomes);
        fprintf(stderr, "\n");
        fprintf(stderr, "       ========== Micro variants ==========\n");
        fprintf(stderr, "\n");

        int32_t vtypes[] = {VT_SNP, VT_MNP, VT_INDEL, 3, 5, 6, 7, 8, 9, 10, 11 ,12 , 13, 14, 15, 0};

        int32_t total_no_micro_variants = 0;
        for (int32_t i=0; i<16; ++i)
        {
            int32_t vtype = vtypes[i];
            if (!VAR_COUNT[POLYMORPHIC][vtype]) continue;
            fprintf(stderr, "       no. of %-21s       : %10d\n", Variant::vtype2string(vtype).c_str(), VAR_COUNT[POLYMORPHIC][vtype]);
            for (int32_t no_alleles=1; no_alleles<=4; ++no_alleles)
            {
                total_no_micro_variants += VAR_COUNT[no_alleles][vtype];
                if (VAR_COUNT[no_alleles][vtype])
                {
                    if (no_alleles==4)
                    {
                        fprintf(stderr, "           >=%d alleles                    : %15d", no_alleles+1, VAR_COUNT[no_alleles][vtype]);
                    }
                    else
                    {
                        fprintf(stderr, "           %d alleles                      : %15d", no_alleles+1, VAR_COUNT[no_alleles][vtype]);
                    }
                    if (vtype&(VT_SNP|VT_MNP)) fprintf(stderr, " (%.2f) [%d/%d]", (float)VAR_TS[no_alleles][vtype]/VAR_TV[no_alleles][vtype], VAR_TS[no_alleles][vtype],VAR_TV[no_alleles][vtype]);
                    if (vtype&VT_INDEL) fprintf(stderr, " (%.2f) [%d/%d]", (float)VAR_INS[no_alleles][vtype]/VAR_DEL[no_alleles][vtype], VAR_INS[no_alleles][vtype], VAR_DEL[no_alleles][vtype]);
                    fprintf(stderr, "\n");
                }
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "       no. of micro variants              : %10d\n", total_no_micro_variants);
        fprintf(stderr, "\n");
        fprintf(stderr, "       ++++++ Other useful categories +++++\n");
        fprintf(stderr, "\n");
        int32_t other_vtypes[3] = {VT_NAIVE_CLUMPED, VT_BLKSUB, VT_CPLXSUB};

        for (int32_t i=0; i<3; ++i)
        {
            int32_t vtype = other_vtypes[i];

            if (!VAR_COUNT[POLYMORPHIC][vtype])  continue;
            std::string variant_desc = vtype==VT_BLKSUB ? "block substitutions" : (vtype==VT_NAIVE_CLUMPED ? "clumped variants" : "complex substitutions");
            fprintf(stderr, "       no. of %-21s       : %10d\n", variant_desc.c_str(), VAR_COUNT[POLYMORPHIC][vtype]);
            for (int32_t no_alleles=1; no_alleles<=4; ++no_alleles)
            {
                if (VAR_COUNT[no_alleles][vtype])
                {
                    if (no_alleles==4)
                    {
                        fprintf(stderr, "           >=%d alleles                    : %15d", no_alleles+1, VAR_COUNT[no_alleles][vtype]);
                    }
                    else
                    {
                        fprintf(stderr, "           %d alleles                      : %15d", no_alleles+1, VAR_COUNT[no_alleles][vtype]);
                    }
                    fprintf(stderr, " (%.2f) [%d/%d]", (float)VAR_TS[no_alleles][vtype]/VAR_TV[no_alleles][vtype], VAR_TS[no_alleles][vtype],VAR_TV[no_alleles][vtype]);
                    if (vtype!=VT_BLKSUB) fprintf(stderr, " (%.2f) [%d/%d]", (float)VAR_INS[no_alleles][vtype]/VAR_DEL[no_alleles][vtype], VAR_INS[no_alleles][vtype], VAR_DEL[no_alleles][vtype]);
                    fprintf(stderr, "\n");
                }
            }
            fprintf(stderr, "\n");
        }

        fprintf(stderr, "       ============== VNTR ===============\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of VNTRs                       : %10d\n", VAR_COUNT[POLYMORPHIC][VT_VNTR]);
        for (int32_t i=0; i<NO_MOTIF_LEN_CATEGORIES; ++i)
        {
            if (VAR_MOTIF_LEN[i])
            {
                if (i<NO_MOTIF_LEN_CATEGORIES-1)
                {
                    fprintf(stderr, "           no. of %d bp motifs                   : %10d\n", i+1, VAR_MOTIF_LEN[i]);
                }
                else
                {
                    fprintf(stderr, "           no. of >%d bp motifs                   : %10d\n", i+1, VAR_MOTIF_LEN[i]);
                }
            }
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "       ======= Structural variants ========\n");
        fprintf(stderr, "\n");
        std::vector<SVNode*> s = sv->enumerate_dfs();
        fprintf(stderr, "       no. of structural variants         : %10d\n", s[0]->count+s[0]->mcount);
        if (s[0]->count)
        {
            fprintf(stderr, "           2 alleles                      : %15d\n", s[0]->count);
            for (size_t i=1; i<s.size(); ++i)
            {
                if (s[i]->count)
                {
                    fprintf(stderr, "            ");
                    for (int32_t j=0; j<s[i]->depth; ++j) std::cerr << "   ";
                    std::string desc = s[i]->sv_type2string(s[i]->type.s);
                    fprintf(stderr, "%s", desc.c_str());
                    for (int32_t j=0; j<30-desc.size()-s[i]->depth*3; ++j) std::cerr << " ";
                    fprintf(stderr, ":           ");
                    for (int32_t j=0; j<s[i]->depth; ++j) std::cerr << "    ";
                    fprintf(stderr, "%6d\n", s[i]->count);
                }
            }
        }

        if (s[0]->mcount)
        {
            fprintf(stderr, "           >=3 alleles                    : %15d\n", s[0]->mcount);
            for (size_t i=1; i<s.size(); ++i)
            {
                if (s[i]->mcount)
                {
                    fprintf(stderr, "            ");
                    for (int32_t j=0; j<s[i]->depth; ++j) std::cerr << "   ";
                    std::string desc = s[i]->sv_type2string(s[i]->type.s);
                    fprintf(stderr, "%s", desc.c_str());
                    for (int32_t j=0; j<30-desc.size()-s[i]->depth*3; ++j) std::cerr << " ";
                    fprintf(stderr, ":           ");
                    for (int32_t j=0; j<s[i]->depth; ++j) std::cerr << "    ";
                    fprintf(stderr, "%6d\n", s[i]->mcount);
                }
            }
        }

        if (sv->mixed_sv_count)
        {
            fprintf(stderr, "            mixed sv                      :  %15d\n", sv->mixed_sv_count);
        }

        fprintf(stderr, "\n");
        fprintf(stderr, "       ========= General summary ==========\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of reference records                  : %10d\n", no_reference);
        fprintf(stderr, "       no. of classified variants                : %10d\n", no_classified_variants);
        fprintf(stderr, "       no. of unclassified variants              : %10d\n", no_unclassified_variants);
        fprintf(stderr, "\n");
        fprintf(stderr, "       no. of VCF records                        : %10d\n", no_records);
        fprintf(stderr, "\n");
    };

    void print_pdf()
    {
        if (output_tabulate_dir=="")
        {
            return;
        }
        else
        {
            if (output_pdf_file=="")
            {
                output_pdf_file = "summary.pdf";
            }
        }

        append_cwd(output_tabulate_dir);

        //generate file
        int32_t ret = mkdir(output_tabulate_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


        std::string filepath = output_tabulate_dir + "/tabulate.tex";
        FILE *out = fopen(filepath.c_str(), "w");

        fprintf(out, "\\PassOptionsToPackage{table}{xcolor}\n");
        fprintf(out, "\\documentclass{article}\n");
        fprintf(out, "\\usepackage{xcolor}\n");
        fprintf(out, "\\begin{document}\n");

        fprintf(out, "\\begin{table}[H]\n");
        fprintf(out, "\\centering\n");
        //fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
        fprintf(out, "\\begin{tabular}{lcr}\n");
        fprintf(out, "no. of samples & :& %d \\\\ \n", no_samples);
        fprintf(out, "no. of chromosomes & : & %d \\\\ \n", no_chromosomes);
        fprintf(out, "\\end{tabular}\n");
        fprintf(out, "\\end{table}\n");

        fprintf(out, "\\begin{table}[H]\n");
        fprintf(out, "\\centering\n");
        //fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
        fprintf(out, "\\begin{tabular}{lrrrrrrrrrr}\n");
        fprintf(out, "\\multicolumn{3}{r}{} & counts & subcounts & ts/tv & ts & tv & ins/del & ins & del \\\\ \n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{3}{l}{SNPs} & %d & & & & & & & \\\\ \n", VAR_COUNT[POLYMORPHIC][VT_SNP]);
        fprintf(out, "\\multicolumn{3}{r}{2 alleles} & & %d & %.2f & %d & %d & & & \\\\ \n",  VAR_COUNT[BIALLELIC][VT_SNP],
                                                                (float)VAR_TS[BIALLELIC][VT_SNP]/VAR_TV[BIALLELIC][VT_SNP],
                                                                 VAR_TS[BIALLELIC][VT_SNP], VAR_TV[BIALLELIC][VT_SNP]);
        fprintf(out, "\\multicolumn{3}{r}{3 alleles} & & %d & %.2f & %d & %d & & & \\\\ \n", VAR_COUNT[TRIALLELIC][VT_SNP],
                                                                 (float)VAR_TS[TRIALLELIC][VT_SNP]/VAR_TV[TRIALLELIC][VT_SNP],
                                                                 VAR_TS[TRIALLELIC][VT_SNP], VAR_TV[TRIALLELIC][VT_SNP]);
        fprintf(out, "\\multicolumn{3}{r}{4 alleles} & & %d & %.2f & %d & %d & & & \\\\ \n", VAR_COUNT[TETRAALLELIC][VT_SNP],
                                                                 (float)VAR_TS[TETRAALLELIC][VT_SNP]/VAR_TV[TETRAALLELIC][VT_SNP],
                                                                 VAR_TS[TETRAALLELIC][VT_SNP], VAR_TV[TETRAALLELIC][VT_SNP]);
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\multicolumn{3}{l}{MNPs} & %d & & & & & & & \\\\ \n", VAR_COUNT[POLYMORPHIC][VT_MNP]);
        fprintf(out, "&&\\multicolumn{1}{r}{2 alleles} & & %d & %.2f & %d & %d & & & \\\\ \n", VAR_COUNT[BIALLELIC][VT_MNP],
                                                                 (float)VAR_TS[BIALLELIC][VT_MNP]/VAR_TV[BIALLELIC][VT_MNP],
                                                                 VAR_TS[BIALLELIC][VT_MNP], VAR_TV[BIALLELIC][VT_MNP]);
        fprintf(out, "&&\\multicolumn{1}{r}{$\\geq$ 3 alleles} & & %d & %.2f & %d & %d & & & \\\\ \n", VAR_COUNT[GE_TRIALLELIC][VT_MNP],
                                                                 (float)VAR_TS[GE_TRIALLELIC][VT_MNP]/VAR_TV[GE_TRIALLELIC][VT_MNP],
                                                                 VAR_TS[GE_TRIALLELIC][VT_MNP], VAR_TV[GE_TRIALLELIC][VT_MNP]);
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\multicolumn{3}{l}{Indel} & %d & & & & & & & \\\\ \n", VAR_COUNT[POLYMORPHIC][VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{2 alleles} & & %d & & & & %.2f & %d & %d \\\\ \n", VAR_COUNT[BIALLELIC][VT_INDEL],
                                                                 (float)VAR_INS[BIALLELIC][VT_INDEL]/VAR_DEL[BIALLELIC][VT_INDEL],
                                                                 VAR_INS[BIALLELIC][VT_INDEL], VAR_DEL[BIALLELIC][VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{$\\geq$ 3 alleles} & & %d & & & & %.2f & %d & %d \\\\ \n", VAR_COUNT[GE_TRIALLELIC][VT_INDEL],
                                                                 (float)VAR_INS[GE_TRIALLELIC][VT_INDEL]/VAR_DEL[GE_TRIALLELIC][VT_INDEL],
                                                                 VAR_INS[GE_TRIALLELIC][VT_INDEL], VAR_DEL[GE_TRIALLELIC][VT_INDEL]);
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\multicolumn{3}{l}{SNP/MNP} & %d & & & & & & & \\\\ \n", VAR_COUNT[POLYMORPHIC][VT_SNP|VT_MNP]);
        fprintf(out, "&&\\multicolumn{1}{r}{3 alleles} & & %d & %.2f & %d & %d & & & \\\\ \n", VAR_COUNT[TRIALLELIC][VT_SNP|VT_MNP],
                                                                 (float)VAR_TS[TRIALLELIC][VT_SNP|VT_MNP]/VAR_TV[TRIALLELIC][VT_SNP|VT_MNP],
                                                                 VAR_TS[TRIALLELIC][VT_SNP|VT_MNP], VAR_TV[TRIALLELIC][VT_SNP|VT_MNP]);
        fprintf(out, "&&\\multicolumn{1}{r}{$\\geq$ 4 alleles} & & %d & %.2f & %d & %d & & & \\\\ \n", VAR_COUNT[GE_TETRAALLELIC][VT_SNP|VT_MNP],
                                                                 (float)VAR_TS[GE_TETRAALLELIC][VT_SNP|VT_MNP]/VAR_TV[GE_TETRAALLELIC][VT_SNP|VT_MNP],
                                                                 VAR_TS[GE_TETRAALLELIC][VT_SNP|VT_MNP], VAR_TV[GE_TETRAALLELIC][VT_SNP|VT_MNP]);
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\multicolumn{3}{l}{SNP/Indel} & %d & & & & &  & & \\\\ \n", VAR_COUNT[POLYMORPHIC][VT_SNP|VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{2 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[BIALLELIC][VT_SNP|VT_INDEL],
                                                                 (float)VAR_TS[BIALLELIC][VT_SNP|VT_INDEL]/VAR_TV[BIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_TS[BIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_TV[BIALLELIC][VT_SNP|VT_INDEL],
                                                                 (float)VAR_INS[BIALLELIC][VT_SNP|VT_INDEL]/VAR_DEL[BIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_INS[BIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_DEL[BIALLELIC][VT_SNP|VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{$\\geq$ 3 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[GE_TRIALLELIC][VT_SNP|VT_INDEL],
                                                                 (float)VAR_TS[GE_TRIALLELIC][VT_SNP|VT_INDEL]/VAR_TV[GE_TRIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_TS[GE_TRIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_TV[GE_TRIALLELIC][VT_SNP|VT_INDEL],
                                                                 (float)VAR_INS[GE_TRIALLELIC][VT_SNP|VT_INDEL]/VAR_DEL[GE_TRIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_INS[GE_TRIALLELIC][VT_SNP|VT_INDEL],
                                                                  VAR_DEL[GE_TRIALLELIC][VT_SNP|VT_INDEL]);
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\multicolumn{3}{l}{MNP/Indel} & %d &  & & & & & & \\\\ \n", VAR_COUNT[POLYMORPHIC][VT_MNP|VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{biallelic} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[BIALLELIC][VT_MNP|VT_INDEL],
                                                                 (float)VAR_TS[BIALLELIC][VT_MNP|VT_INDEL]/VAR_TV[BIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_TS[BIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_TV[BIALLELIC][VT_MNP|VT_INDEL],
                                                                 (float)VAR_INS[BIALLELIC][VT_MNP|VT_INDEL]/VAR_DEL[BIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_INS[BIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_DEL[BIALLELIC][VT_MNP|VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{$\\geq$ 3 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[GE_TRIALLELIC][VT_MNP|VT_INDEL],
                                                                 (float)VAR_TS[GE_TRIALLELIC][VT_MNP|VT_INDEL]/VAR_TV[GE_TRIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_TS[GE_TRIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_TV[GE_TRIALLELIC][VT_MNP|VT_INDEL],
                                                                 (float)VAR_INS[GE_TRIALLELIC][VT_MNP|VT_INDEL]/VAR_DEL[GE_TRIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_INS[GE_TRIALLELIC][VT_MNP|VT_INDEL],
                                                                 VAR_DEL[GE_TRIALLELIC][VT_MNP|VT_INDEL]);
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\multicolumn{3}{l}{SNP/MNP/Indel} & %d &  & & & & & & \\\\ \n", VAR_COUNT[POLYMORPHIC][VT_SNP|VT_MNP|VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{3 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 (float)VAR_TS[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL]/VAR_TV[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_TS[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_TV[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 (float)VAR_INS[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL]/VAR_DEL[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_INS[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_DEL[TRIALLELIC][VT_SNP|VT_MNP|VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{4 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 (float)VAR_TS[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL]/VAR_TV[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_TS[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_TV[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 (float)VAR_INS[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL]/VAR_DEL[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_INS[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_DEL[TETRAALLELIC][VT_SNP|VT_MNP|VT_INDEL]);
        fprintf(out, "&&\\multicolumn{1}{r}{$\\geq$ 5 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 (float)VAR_TS[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL]/VAR_TV[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_TS[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_TV[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 (float)VAR_INS[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL]/VAR_DEL[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_INS[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL],
                                                                 VAR_DEL[GE_PENTAALLELIC][VT_SNP|VT_MNP|VT_INDEL]);
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\hline\n");
        fprintf(out, "\\multicolumn{11}{l}{} \\\\\n");
        fprintf(out, "\\multicolumn{3}{l}{Clumped} & %d & & & & & & & \\\\ \n",  VAR_COUNT[POLYMORPHIC][VT_CLUMPED]);
        fprintf(out, "&&\\multicolumn{1}{r}{2 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[BIALLELIC][VT_CLUMPED],
                                                                 (float)VAR_TS[BIALLELIC][VT_CLUMPED]/VAR_TV[BIALLELIC][VT_CLUMPED],
                                                                 VAR_TS[BIALLELIC][VT_CLUMPED],
                                                                 VAR_TV[BIALLELIC][VT_CLUMPED],
                                                                 (float)VAR_INS[BIALLELIC][VT_CLUMPED]/VAR_DEL[BIALLELIC][VT_CLUMPED],
                                                                 VAR_INS[BIALLELIC][VT_CLUMPED],
                                                                 VAR_DEL[BIALLELIC][VT_CLUMPED]);
        fprintf(out, "&&\\multicolumn{1}{r}{3 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[TRIALLELIC][VT_CLUMPED],
                                                                 (float)VAR_TS[TRIALLELIC][VT_CLUMPED]/VAR_TV[TRIALLELIC][VT_CLUMPED],
                                                                 VAR_TS[TRIALLELIC][VT_CLUMPED],
                                                                 VAR_TV[TRIALLELIC][VT_CLUMPED],
                                                                 (float)VAR_INS[TRIALLELIC][VT_CLUMPED]/VAR_DEL[TRIALLELIC][VT_CLUMPED],
                                                                 VAR_INS[TRIALLELIC][VT_CLUMPED],
                                                                 VAR_DEL[TRIALLELIC][VT_CLUMPED]);
        fprintf(out, "&&\\multicolumn{1}{r}{4 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[TETRAALLELIC][VT_CLUMPED],
                                                                 (float)VAR_TS[TETRAALLELIC][VT_CLUMPED]/VAR_TV[TETRAALLELIC][VT_CLUMPED],
                                                                 VAR_TS[TETRAALLELIC][VT_CLUMPED],
                                                                 VAR_TV[TETRAALLELIC][VT_CLUMPED],
                                                                 (float)VAR_INS[TETRAALLELIC][VT_CLUMPED]/VAR_DEL[TETRAALLELIC][VT_CLUMPED],
                                                                 VAR_INS[TETRAALLELIC][VT_CLUMPED],
                                                                 VAR_DEL[TETRAALLELIC][VT_CLUMPED]);
        fprintf(out, "&&\\multicolumn{1}{r}{$\\geq$ 5 alleles} & & %d & %.2f & %d & %d & %.2f & %d & %d \\\\ \n", VAR_COUNT[GE_PENTAALLELIC][VT_CLUMPED],
                                                                 (float)VAR_TS[GE_PENTAALLELIC][VT_CLUMPED]/VAR_TV[GE_PENTAALLELIC][VT_CLUMPED],
                                                                 VAR_TS[GE_PENTAALLELIC][VT_CLUMPED],
                                                                 VAR_TV[GE_PENTAALLELIC][VT_CLUMPED],
                                                                 (float)VAR_INS[GE_PENTAALLELIC][VT_CLUMPED]/VAR_DEL[GE_PENTAALLELIC][VT_CLUMPED],
                                                                 VAR_INS[GE_PENTAALLELIC][VT_CLUMPED],
                                                                 VAR_DEL[GE_PENTAALLELIC][VT_CLUMPED]);
        fprintf(out, "\\end{tabular}\n");
        fprintf(out, "\\end{table}\n");

        fprintf(out, "\\begin{table}[H]\n");
        fprintf(out, "\\centering\n");
        //fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
        fprintf(out, "\\begin{tabular}{lcr}\n");
        fprintf(out, "no. of reference& :& %d \\\\ \n", VAR_COUNT[MONOMORPHIC][VT_REF]);
        fprintf(out, "no. of observed & : & %d \\\\ \n", no_records);
        fprintf(out, "no. of unclassified & : & %d \\\\ \n", no_records-no_classified_variants);
        fprintf(out, "\\end{tabular}\n");
        fprintf(out, "\\end{table}\n");

        fprintf(out, "\\end{document}\n");

        fclose(out);

        std::string cmd = "cd "  + output_tabulate_dir + "; pdflatex tabulate.tex > run.log; mv tabulate.pdf " + output_pdf_file;
        std::cerr << cmd << "\n";

        int32_t sys_ret = system(cmd.c_str());
    };

    ~Igor()
    {
        delete odr;
        bcf_destroy(v);

        delete vm;
        delete sv;
    };

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
    igor.print_pdf();
};
