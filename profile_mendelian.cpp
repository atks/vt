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

#include "profile_mendelian.h"

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
    std::string input_ped_file;
    std::string output_tabulate_dir;
    std::string output_pdf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    int32_t min_depth;
    float_t min_gq;
    bool ignore_non_variants;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;

    ///////////////
    //general use//
    ///////////////
    kstring_t variant;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    uint32_t no_trios;
    uint32_t no_variants;
    uint32_t no_failed_min_depth;
    int32_t trio_genotypes[3][3][3];

    /////////
    //tools//
    /////////
    Pedigree *pedigree;
    VariantManip *vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Profile Mendelian Errors.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_ped_file("p", "p", "pedigree file", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_tabulate_dir("x", "x", "output latex directory [tabulate_mendelian]", false, "tabulate_mendelian", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output pdf file [mendelian.pdf]", false, "mendelian.pdf", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_min_depth("d", "d", "minimum depth", false, 0, "str", cmd);
            TCLAP::ValueArg<float> arg_min_gq("q", "q", "minimum genotype quality", false, 2, "str", cmd);
            TCLAP::SwitchArg arg_ignore_non_variants("n", "n", "ignore non variants", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            input_ped_file = arg_input_ped_file.getValue();
            fexp = arg_fexp.getValue();
            output_tabulate_dir = arg_output_tabulate_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
            min_depth = arg_min_depth.getValue();
            min_gq = arg_min_gq.getValue();
            ignore_non_variants = arg_ignore_non_variants.getValue();
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

        if (bcf_hdr_nsamples(odr->hdr)==0)
        {
            fprintf(stderr, "[%s:%d %s] No samples in VCF file: %s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
            exit(1);
        }

        ///////////////////////////
        //ped file initialization//
        ///////////////////////////
        pedigree = new Pedigree(input_ped_file);

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp!="";

        ///////////////
        //general use//
        ///////////////
        variant = {0,0,0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_trios = 0;
        no_variants = 0;
        no_failed_min_depth = 0;

        for (int32_t i=0; i<3; ++i)
        {
            for (int32_t j=0; j<3; ++j)
            {
                for (int32_t k=0; k<3; ++k)
                {
                    trio_genotypes[i][j][k] = 0;
                }
            }
        }

        /////////
        //tools//
        /////////
        vm = new VariantManip();
    }

    void profile_mendelian()
    {
        bcf_hdr_t *h = odr->hdr;
        bcf1_t *v = bcf_init1();

        Variant variant;

        std::vector<Trio> trios;

        for (size_t i=0; i<pedigree->trios.size(); ++i)
        {
            int32_t f = bcf_hdr_id2int(h, BCF_DT_SAMPLE, pedigree->trios[i].father.c_str());
            int32_t m = bcf_hdr_id2int(h, BCF_DT_SAMPLE, pedigree->trios[i].mother.c_str());
            int32_t c = bcf_hdr_id2int(h, BCF_DT_SAMPLE, pedigree->trios[i].child.c_str());

            if (f>=0 && m>=0 && c>=0)
            {
                pedigree->trios[i].father_index = f;
                pedigree->trios[i].mother_index = m;
                pedigree->trios[i].child_index = c;
                trios.push_back(pedigree->trios[i]);
            }
        }

        no_trios = trios.size();

        int32_t missing = 0;
        int32_t mendel_homalt_err = 0;

        int32_t nsample = bcf_hdr_get_n_sample(h);
        int32_t *gts = NULL;
        int32_t *dps = (int32_t *) malloc(nsample*sizeof(int32_t));
        int32_t n = 0;
        int32_t n_dp = nsample;

        while(odr->read(v))
        {
            bcf_unpack(v, BCF_UN_IND);
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);

            //if (bcf_get_n_allele(v)!=2 || vtype!=VT_INDEL || bcf_has_filter(odr->hdr, v, const_cast<char*>("PASS"))!=1)
            if (bcf_get_n_allele(v)!=2)
            {
                continue;
            }

            if (filter_exists)
            {
                vm->classify_variant(odr->hdr, v, variant);
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }

            int k = bcf_get_genotypes(h, v, &gts, &n);
            int r = bcf_get_format_int32(h, v, "DP", &dps, &n_dp);

            if (r==-1)
            {
                r = bcf_get_format_int32(h, v, "NR", &dps, &n_dp);

                if (r==-1)
                {
                    for (uint32_t i=0; i<nsample; ++i)
                    {
                        dps[i] = min_depth;
                    }
                }
            }

            bool variant_used = false;

            for (int32_t i =0; i< trios.size(); ++i)
            {
                int32_t j = trios[i].father_index;
                int32_t f1 = bcf_gt_allele(gts[(j<<1)]);
                int32_t f2 = bcf_gt_allele(gts[(j<<1)+1]);
                int32_t min_dp = dps[j];

                j = trios[i].mother_index;
                int32_t m1 = bcf_gt_allele(gts[(j<<1)]);
                int32_t m2 = bcf_gt_allele(gts[(j<<1)+1]);
                min_dp = dps[j]<min_dp ? dps[j] : min_dp;

                j = trios[i].child_index;
                int32_t c1 = bcf_gt_allele(gts[(j<<1)]);
                int32_t c2 = bcf_gt_allele(gts[(j<<1)+1]);
                min_dp = dps[j]<min_dp ? dps[j] : min_dp;

                if (min_dp<min_depth)
                {
                    ++no_failed_min_depth;
                    continue;
                }

                if (!(f1<0 || f2<0 || m1<0 || m2<0 || c1<0 || c2<0))
                {
                    if (!ignore_non_variants || (f1+f2+m1+m2+c1+c2!=0))
                    {
                        ++trio_genotypes[f1+f2][m1+m2][c1+c2];
                        variant_used = true;
                    }
                }
            }
            if (variant_used) ++no_variants;
        }

        free(gts);
        if (dps) free(dps);
        odr->close();
    };

    void print_options()
    {
        std::clog << "profile_mendelian v" << version << "\n\n";
        std::clog << "options:     input VCF file            " << input_vcf_file << "\n";
        std::clog << "         [p] input PED file            " << input_ped_file << "\n";
        std::clog << "         [d] minimum depth             " << min_depth << "\n";
        print_str_op("         [f] filter                    ", fexp);
        print_str_op("         [x] output tabulate directory ", output_tabulate_dir);
        print_str_op("         [y] output pdf file           ", output_pdf_file);
        print_int_op("         [i] intervals                 ", intervals);
        std::clog << "\n";
    }

    float get_error_rate(int32_t gt[3][3][3], int32_t f, int32_t m, int32_t collapse)
    {
        if (collapse==-1) //ALL
        {
            float total = 0;
            float error_count = 0;

            for (int32_t i=0; i<3; ++i)
            {
                for (int32_t j=0; j<3; ++j)
                {
                    f = i;
                    m = j;
                    total += gt[f][m][0] + gt[f][m][1] + gt[f][m][2];

                    if (f==m && f!=1)//0/0,2/2
                    {
                        error_count += gt[f][m][1] + gt[f][m][2-f];
                    }
                    else if (abs(f-m)==2) //0/2,2/0
                    {
                        error_count += gt[f][m][0] + gt[f][m][2];
                    }
                    else if (abs(f-m)==1)//1/2,2/1,1/0,0/1
                    {
                        error_count += gt[f][m][f==1?2-m:2-f];
                    }
                    else if (f==1 && m==1)//1/1
                    {
                        error_count += 0;
                    }
                }
            }

            return error_count/total*100;
        }
        else if (collapse==0) //no collapse
        {
            float total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
            float error_count = 0;
            if (f==m && f!=1)//0/0,2/2
            {
                error_count = gt[f][m][1] + gt[f][m][2-f];
            }
            else if (abs(f-m)==2) //0/2,2/0
            {
                error_count = gt[f][m][0] + gt[f][m][2];
            }
            else if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                error_count = gt[f][m][f==1?2-m:2-f];
            }
            else if (f==1 && m==1)//1/1
            {
                error_count = 0;
            }

            return error_count/total*100;
        }
        else if (collapse==1) //ignoring father/mother
        {
            float total = 0;
            float error_count = 0;
            if (f==m && f!=1)//0/0,2/2
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                error_count = gt[f][m][1] + gt[f][m][2-f];
            }
            else if (abs(f-m)==2) //0/2,2/0
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                total += gt[m][f][0] + gt[m][f][1] + gt[m][f][2];
                error_count = gt[f][m][0] + gt[f][m][2];
                error_count += gt[m][f][0] + gt[m][f][2];
            }
            else if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                total += gt[m][f][0] + gt[m][f][1] + gt[m][f][2];
                error_count = gt[f][m][f==1?2-m:2-f];
                error_count += gt[m][f][f==1?2-m:2-f];
            }
            else if (f==1 && m==1)//1/1
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                error_count = 0;
            }

            return error_count/total*100;
        }
        else if (collapse==2) //ignoring allele and father/mother
        {
            float total = 0;
            float error_count = 0;
            if (f==m && f!=1)//0/0,2/2
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                total += gt[2-f][2-m][0] + gt[2-f][2-m][1] + gt[2-f][2-m][2];
                error_count = gt[f][m][1] + gt[f][m][2-f];
                error_count += gt[2-f][2-m][1] + gt[2-f][2-m][f];
            }
            else if (abs(f-m)==2) //0/2,2/0
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                total += gt[m][f][0] + gt[m][f][1] + gt[m][f][2];
                error_count = gt[f][m][0] + gt[f][m][2];
                error_count += gt[m][f][0] + gt[m][f][2];
            }
            else if (abs(f-m)==1)// 1/2,2/1,1/0,0/1
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                total += gt[m][f][0] + gt[m][f][1] + gt[m][f][2];
                total += gt[2-f][2-m][0] + gt[2-f][2-m][1] + gt[2-f][2-m][2];
                total += gt[2-m][2-f][0] + gt[2-m][2-f][1] + gt[2-m][2-f][2];

                error_count = gt[f][m][f==1?2-m:2-f];
                error_count += gt[m][f][f==1?2-m:2-f];
                error_count += gt[2-f][2-m][f==1?m:f];
                error_count += gt[2-m][2-f][f==1?m:f];
            }
            else if (f==1 && m==1)//1/1
            {
                total = gt[f][m][0] + gt[f][m][1] + gt[f][m][2];
                error_count = 0;
            }

            return error_count/total*100;
        }

        return NAN;
    };

    float get_homhet_ratio(int32_t gt[3][3][3], int32_t f, int32_t m, int32_t collapse)
    {
        float hom = 0;
        float het = 0;
        if (collapse==0) //ALL
        {
            if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                hom = gt[f][m][f==1?m:f];
                het = gt[f][m][1];
            }
            else if (f==1 && m==1)//1/1
            {
                hom = gt[f][m][0]+gt[f][m][2];
                het = gt[f][m][1];
            }
        }
        else if (collapse==1) //no collapse
        {
            if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                hom = gt[f][m][f==1?m:f];
                hom += gt[m][f][f==1?m:f];
                het = gt[f][m][1];
                het += gt[m][f][1];
            }
            else if (f==1 && m==1)//1/1
            {
                hom = gt[f][m][0]+gt[f][m][2];
                het = gt[f][m][1];
            }
        }
        else if (collapse==2) //no collapse
        {
            if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                hom = gt[f][m][f==1?m:f];
                hom += gt[m][f][f==1?m:f];
                hom += gt[2-f][2-m][f==1?2-m:2-f];
                hom += gt[2-m][2-f][f==1?2-m:2-f];
                het = gt[f][m][1];
                het += gt[m][f][1];
                het += gt[2-f][2-m][1];
                het += gt[2-m][2-f][1];
            }
            else if (f==1 && m==1)//1/1
            {
                hom = gt[f][m][0]+gt[f][m][2];
                het = gt[f][m][1];
            }
        }

        return (het==0 ? NAN : hom/het);
    };

    float get_homhet_proportion(int32_t gt[3][3][3], int32_t f, int32_t m, int32_t collapse)
    {
        float hom = 0;
        float het = 0;
        if (collapse==0) //ALL
        {
            if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                hom = gt[f][m][f==1?m:f];
                het = gt[f][m][1];
            }
            else if (f==1 && m==1)//1/1
            {
                hom = gt[f][m][0]+gt[f][m][2];
                het = gt[f][m][1];
            }
        }
        else if (collapse==1) //no collapse
        {
            if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                hom = gt[f][m][f==1?m:f];
                hom += gt[m][f][f==1?m:f];
                het = gt[f][m][1];
                het += gt[m][f][1];
            }
            else if (f==1 && m==1)//1/1
            {
                hom = gt[f][m][0]+gt[f][m][2];
                het = gt[f][m][1];
            }
        }
        else if (collapse==2) //no collapse
        {
            if (abs(f-m)==1)//1/2,2/1,1/0,0/1
            {
                hom = gt[f][m][f==1?m:f];
                hom += gt[m][f][f==1?m:f];
                hom += gt[2-f][2-m][f==1?2-m:2-f];
                hom += gt[2-m][2-f][f==1?2-m:2-f];
                het = gt[f][m][1];
                het += gt[m][f][1];
                het += gt[2-f][2-m][1];
                het += gt[2-m][2-f][1];
            }
            else if (f==1 && m==1)//1/1
            {
                hom = gt[f][m][0]+gt[f][m][2];
                het = gt[f][m][1];
            }
        }

        return ((het+hom)==0 ? NAN : het/(hom+het)*100);
    };

    void print_pdf()
    {
        append_cwd(output_tabulate_dir);

        //generate file
        int32_t ret = mkdir(output_tabulate_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        std::string filepath = output_tabulate_dir + "/tabulate.tex";
        FILE *out = fopen(filepath.c_str(), "w");

        std::string g2s[3] = {"R/R","R/A","A/A"};

        fprintf(out, "\\PassOptionsToPackage{table}{xcolor}\n");
        fprintf(out, "\\documentclass{beamer}\n");
        fprintf(out, "\\begin{document}\n");
        fprintf(out, "\\begin{frame}{Mendel Error (All parental genotypes)}\n");
        fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
        fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
        fprintf(out, "\\begin{tabular}{ccrrrrrr}\n");
        fprintf(out, "\\rowcolor{blue!50}\n");
        fprintf(out, "father& mother & R/R & R/A & A/A & Error \\%% & HomHet & Het \\%%\\\\ \n");
        for (int32_t i=0; i<3; ++i)
        {
            for (int32_t j=0; j<3; ++j)
            {
                fprintf(out, "%s&%s&%10d&%10d&%10d&%6.2f&%4.2f&%5.2f \\\\ \n", g2s[i].c_str(), g2s[j].c_str(), trio_genotypes[i][j][0], trio_genotypes[i][j][1], trio_genotypes[i][j][2], get_error_rate(trio_genotypes, i, j, 0), get_homhet_ratio(trio_genotypes, i, j, 0), get_homhet_proportion(trio_genotypes, i, j, 0));
            }
        }
        fprintf(out, "\\end{tabular}}\n");
        fprintf(out, "\\end{frame}\n");

        fprintf(out, "\\begin{frame}{Mendel Error (Collapsed genotypes)}\n");
        fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
        fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
        fprintf(out, "\\begin{tabular}{ccrrrrrr}\n");
        fprintf(out, "\\rowcolor{blue!50}\n");
        fprintf(out, "\\multicolumn{2}{c}{Parental}  & R/R & R/A & A/A & Error \\%% & HomHet & Het \\%%\\\\ \n");
        for (int32_t i=0; i<3; ++i)
        {
            for (int32_t j=0; j<3; ++j)
            {
                if (i!=j)
                {
                    if (i<j)
                    {
                        int32_t rr = trio_genotypes[i][j][0] + trio_genotypes[j][i][0];
                        int32_t ra = trio_genotypes[i][j][1] + trio_genotypes[j][i][1];
                        int32_t aa = trio_genotypes[i][j][2] + trio_genotypes[j][i][2];
                        fprintf(out, "%s&%s&%10d&%10d&%10d&%6.2f&%4.2f&%5.2f \\\\ \n", g2s[i].c_str(), g2s[j].c_str(), rr, ra, aa, get_error_rate(trio_genotypes, i, j, 1), get_homhet_ratio(trio_genotypes, i, j, 1), get_homhet_proportion(trio_genotypes, i, j, 1));
                    }
                }
                else
                {
                    fprintf(out, "%s&%s&%10d&%10d&%10d&%6.2f&%4.2f&%5.2f \\\\ \n", g2s[i].c_str(), g2s[j].c_str(), trio_genotypes[i][j][0], trio_genotypes[i][j][1], trio_genotypes[i][j][2], get_error_rate(trio_genotypes, i, j, 1), get_homhet_ratio(trio_genotypes, i, j, 1), get_homhet_proportion(trio_genotypes, i, j, 1));
                }
            }
        }
        fprintf(out, "\n");
        fprintf(out, "\\end{tabular}}\n");
        fprintf(out, "\\end{frame}\n");

        fprintf(out, "\\begin{frame}{Mendel Error (Collapse parental allelotypes)}\n");
        fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
        fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
        fprintf(out, "\\begin{tabular}{ccrrrrrr}\n");
        fprintf(out, "\\rowcolor{blue!50}\n");
        fprintf(out, "\\multicolumn{2}{c}{Parental}  & R/R & R/A & A/A & Error \\%% & HomHet & Het \\%%\\\\ \n");
        int32_t i,j, rr, ra, aa;
        i=0; j=0;
        rr = trio_genotypes[i][j][0] + trio_genotypes[2][2][0];
        ra = trio_genotypes[i][j][1] + trio_genotypes[2][2][1];
        aa = trio_genotypes[i][j][2] + trio_genotypes[2][2][2];
        fprintf(out, "%s&%s&%10d&%10d&%10d&%6.2f&%4.2f&%5.2f \\\\ \n", "HOM", "HOM", rr, ra, aa, get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        i=0; j=1;
        rr = trio_genotypes[i][j][0] + trio_genotypes[j][i][0] + trio_genotypes[2-i][j][0] + trio_genotypes[j][2-i][0];
        ra = trio_genotypes[i][j][1] + trio_genotypes[j][i][1] + trio_genotypes[2-i][j][1] + trio_genotypes[j][2-i][1];
        aa = trio_genotypes[i][j][2] + trio_genotypes[j][i][2] + trio_genotypes[2-i][j][2] + trio_genotypes[j][2-i][2];
        fprintf(out, "%s&%s&%10d&%10d&%10d&%6.2f&%4.2f&%5.2f \\\\ \n", "HOM", "HET", rr, ra, aa, get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        i=1; j=1;
        fprintf(out, "%s&%s&%10d&%10d&%10d&%6.2f&%4.2f&%5.2f \\\\ \n", "HET", "HET", trio_genotypes[i][j][0], trio_genotypes[i][j][1], trio_genotypes[i][j][2], get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        i=0; j=2;
        rr = trio_genotypes[i][j][0] + trio_genotypes[j][i][0];
        ra = trio_genotypes[i][j][1] + trio_genotypes[j][i][1];
        aa = trio_genotypes[i][j][2] + trio_genotypes[j][i][2];
        fprintf(out, "%s&%s&%10d&%10d&%10d&%6.2f&%4.2f&%5.2f \\\\ \n", "HOMREF", "HOMALT", rr, ra, aa, get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        fprintf(out, "\n");
        fprintf(out, "\\end{tabular}}\n");
        fprintf(out, "\\end{frame}\n");

        fprintf(out, "\\begin{frame}{Mendel Error Overall Summary}\n");
        fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
        fprintf(out, "\\begin{tabular}{lr}\n");
        fprintf(out, "total mendelian error & %7.3f \\%% \\\\ \n", get_error_rate(trio_genotypes, -1, -1, -1));
        fprintf(out, "no. of trios     & %d \\\\ \n", no_trios);
        fprintf(out, "no. of variants  & %d \\\\ \n", no_variants);
        fprintf(out, "\n");
        fprintf(out, "\\end{tabular}}\n");
        fprintf(out, "\\end{frame}\n");

        fprintf(out, "\\end{document}\n");
        fprintf(out, "\n");

        fclose(out);

        std::string cmd = "cd "  + output_tabulate_dir + "; pdflatex tabulate.tex > run.log; mv tabulate.pdf " + output_pdf_file;

        int32_t sys_ret = system(cmd.c_str());
    };

    void print_stats()
    {
        std::string g2s[3] = {"R/R","R/A","A/A"};

        fprintf(stderr, "\n");
        fprintf(stderr, "     Mendelian Errors\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "     Father Mother       R/R          R/A          A/A    Error(%%) HomHet    Het(%%)\n");
        for (int32_t i=0; i<3; ++i)
        {
            for (int32_t j=0; j<3; ++j)
            {
                fprintf(stderr, "     %s    %s   %10d   %10d   %10d   %6.2f      %4.2f  %5.2f\n", g2s[i].c_str(), g2s[j].c_str(), trio_genotypes[i][j][0], trio_genotypes[i][j][1], trio_genotypes[i][j][2], get_error_rate(trio_genotypes, i, j, 0), get_homhet_ratio(trio_genotypes, i, j, 0), get_homhet_proportion(trio_genotypes, i, j, 0));
            }
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "     Parental            R/R          R/A          A/A    Error(%%) HomHet    Het(%%)\n");
        for (int32_t i=0; i<3; ++i)
        {
            for (int32_t j=0; j<3; ++j)
            {
                if (i!=j)
                {
                    if (i<j)
                    {
                        int32_t rr = trio_genotypes[i][j][0] + trio_genotypes[j][i][0];
                        int32_t ra = trio_genotypes[i][j][1] + trio_genotypes[j][i][1];
                        int32_t aa = trio_genotypes[i][j][2] + trio_genotypes[j][i][2];
                        fprintf(stderr, "     %s    %s   %10d   %10d   %10d   %6.2f      %4.2f  %5.2f\n", g2s[i].c_str(), g2s[j].c_str(), rr, ra, aa, get_error_rate(trio_genotypes, i, j, 1), get_homhet_ratio(trio_genotypes, i, j, 1), get_homhet_proportion(trio_genotypes, i, j, 1));
                    }
                }
                else
                {
                    fprintf(stderr, "     %s    %s   %10d   %10d   %10d   %6.2f      %4.2f  %5.2f\n", g2s[i].c_str(), g2s[j].c_str(), trio_genotypes[i][j][0], trio_genotypes[i][j][1], trio_genotypes[i][j][2], get_error_rate(trio_genotypes, i, j, 1), get_homhet_ratio(trio_genotypes, i, j, 1), get_homhet_proportion(trio_genotypes, i, j, 1));
                }
            }
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "     Parental            R/R          R/A          A/A    Error(%%) HomHet    Het(%%)\n");
        int32_t i,j, rr, ra, aa;
        i=0; j=0;
        rr = trio_genotypes[i][j][0] + trio_genotypes[2][2][0];
        ra = trio_genotypes[i][j][1] + trio_genotypes[2][2][1];
        aa = trio_genotypes[i][j][2] + trio_genotypes[2][2][2];
        fprintf(stderr, "     %s    %s   %10d   %10d   %10d   %6.2f      %4.2f  %5.2f\n", "HOM", "HOM", rr, ra, aa, get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        i=0; j=1;
        rr = trio_genotypes[i][j][0] + trio_genotypes[j][i][0] + trio_genotypes[2-i][j][0] + trio_genotypes[j][2-i][0];
        ra = trio_genotypes[i][j][1] + trio_genotypes[j][i][1] + trio_genotypes[2-i][j][1] + trio_genotypes[j][2-i][1];
        aa = trio_genotypes[i][j][2] + trio_genotypes[j][i][2] + trio_genotypes[2-i][j][2] + trio_genotypes[j][2-i][2];
        fprintf(stderr, "     %s    %s   %10d   %10d   %10d   %6.2f      %4.2f  %5.2f\n", "HOM", "HET", rr, ra, aa, get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        i=1; j=1;
        fprintf(stderr, "     %s    %s   %10d   %10d   %10d   %6.2f      %4.2f  %5.2f\n", "HET", "HET", trio_genotypes[i][j][0], trio_genotypes[i][j][1], trio_genotypes[i][j][2], get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        i=0; j=2;
        rr = trio_genotypes[i][j][0] + trio_genotypes[j][i][0];
        ra = trio_genotypes[i][j][1] + trio_genotypes[j][i][1];
        aa = trio_genotypes[i][j][2] + trio_genotypes[j][i][2];
        fprintf(stderr, "     %s %s%10d   %10d   %10d   %6.2f      %4.2f  %5.2f\n", "HOMREF", "HOMALT", rr, ra, aa, get_error_rate(trio_genotypes, i, j, 2), get_homhet_ratio(trio_genotypes, i, j, 2), get_homhet_proportion(trio_genotypes, i, j, 2));
        fprintf(stderr, "\n");
        fprintf(stderr, "     total mendelian error : %7.3f%%\n", get_error_rate(trio_genotypes, -1, -1, -1));
        fprintf(stderr, "\n");
        fprintf(stderr, "     no. of trios     : %d\n", no_trios);
        fprintf(stderr, "     no. of variants  : %d\n", no_variants);
        fprintf(stderr, "\n");
        fprintf(stderr, "     no. of trio-sites that fail min depth  : %d\n", no_failed_min_depth);
        fprintf(stderr, "\n");
    };

    ~Igor()
    {
    };

    private:
};

}

void profile_mendelian(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_mendelian();
    igor.print_stats();
    igor.print_pdf();
}

