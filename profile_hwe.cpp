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

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRAFTY OF AFY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRAFTIES OF MERCHAFTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AFD NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR AFY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AF HWE_LPVALTION OF CONTRHWE_LPVALT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include "profile_hwe.h"

namespace
{

KHASH_MAP_INIT_INT(32, int32_t)

struct hwe_t
{
    float hwe_lpval;
    float maf;
    int32_t no_alleles;
    int32_t pass;
};

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_dir;
    std::string output_pdf_file;
    const char* HWE_LPVAL;
    const char* AF;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;

    ///////////////
    //general use//
    ///////////////
    std::vector<hwe_t> hwe_pts;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    int32_t no_variants;

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
            std::string desc = "Plot Allele Frequency Spectrum.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_dir("x", "x", "output directory []", false, "plot_hwe", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output PDF file []", false, "hwe.pdf", "str", cmd);
            TCLAP::ValueArg<std::string> arg_HWE_LPVAL("h", "hwe_lpval", "HWE_LPVAL tag [HWE_LPVAL]", false, "HWE_LPVAL", "str", cmd);
            TCLAP::ValueArg<std::string> arg_AF("a", "af", "AF tag [AF]", false, "AF", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            fexp = arg_fexp.getValue();
            output_dir = arg_output_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
            HWE_LPVAL = strdup(arg_HWE_LPVAL.getValue().c_str());
            AF = strdup(arg_AF.getValue().c_str());
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

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        /////////
        //tools//
        /////////
        vm = new VariantManip();
    }

    void profile_hwe()
    {
        bcf1_t *v = bcf_init1();

        Variant variant;
        float *hwe_lpval=NULL, *af=NULL;
        int32_t n_hwe_lpval=0, n_af=0;

        while(odr->read(v))
        {
            bcf_unpack(v, BCF_UN_ALL);

            if (filter_exists)
            {
                int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }

            int32_t pass = bcf_has_filter(odr->hdr, v, const_cast<char*>("PASS"));

            int32_t ret1 = bcf_get_info_float(odr->hdr, v, HWE_LPVAL, &hwe_lpval, &n_hwe_lpval);
            int32_t ret2 = bcf_get_info_float(odr->hdr, v, AF, &af, &n_af);
            int32_t no_alleles = bcf_get_n_allele(v);
            float maf = 0;

            if (ret1<0||ret2<0)
            {
                continue;
            }

            if (no_alleles==2)
            {
                maf = af[0]>0.5? 1-af[0]: af[0];
            }
            else
            {
                for (size_t i=0; i<n_af; ++i)
                {
                    maf += af[i];
                }

                maf = maf>0.5? 1-maf: maf;
            }

            hwe_t h = {hwe_lpval[0], maf, no_alleles, pass};
            hwe_pts.push_back(h);

            ++no_variants;
        }

        if (n_hwe_lpval) free(hwe_lpval);
        if (n_af) free(af);

        odr->close();
    };

    void print_options()
    {
        std::clog << "plot_hwe v" << version << "\n\n";
        std::clog << "options:     input VCF file         " << input_vcf_file << "\n";
        std::clog << "         [h] HWE_LPVAL tag          " << HWE_LPVAL << "\n";
        std::clog << "         [a] AF tag                 " << AF << "\n";
        print_str_op("         [x] output directory       ", output_dir);
        print_str_op("         [y] output pdf file        ", output_pdf_file);
        print_int_op("         [i] intervals              ", intervals);
        std::clog << "\n";
    }

    void print_pdf()
    {
        append_cwd(output_dir);

        //create directory
        mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        //create data file
        std::string file_path = output_dir + "/data.txt";
        FILE *out = fopen(file_path.c_str(), "w");

        fprintf(out, "hwe_lpval\tmaf\tno_alleles\tpass\n");

        for (size_t i=0; i<hwe_pts.size(); ++i)
        {
            fprintf(out, "%f\t%f\t%d\t%d\n", hwe_pts[i].hwe_lpval, hwe_pts[i].maf, hwe_pts[i].no_alleles, hwe_pts[i].pass);
        }

        fclose(out);

        //create r script
        file_path = output_dir + "/plot.r";
        out = fopen(file_path.c_str(), "w");

        fprintf(out, "setwd(\"%s\")\n", output_dir.c_str());
        fprintf(out, "\n");
        fprintf(out, "data = read.table(\"data.txt\", header=T)\n");
        fprintf(out, "data$hwe_pval=exp(data$hwe_lpval)\n");

        fprintf(out, "plot_hist <- function(data, title, bar_color) {\n");
        fprintf(out, "h = hist(data$hwe_pval, breaks=50, plot=F)\n");
        fprintf(out, "b = barplot(h$counts, log=\"y\", col=bar_color, xlab=\"P-value\", ylab=\"Count\", main=title)\n");
        fprintf(out, "axis(1,c(b[1]-(b[2]-b[1])/2,(b[50]-b[1])/2,b[50]+(b[2]-b[1])/2),c(0,0.5,1.0))\n");
        fprintf(out, "}\n");

        fprintf(out, "pdf(\"%s\",7,8)\n", output_pdf_file.c_str());
        fprintf(out, "\n");
        fprintf(out, "layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==1 & no_alleles==2)\n");
        fprintf(out, "plot_hist(data.subset,\"Passed Biallelic Indels\", rgb(0,0,1,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==1 & no_alleles==2 & maf>0.05)\n");
        fprintf(out, "plot_hist(data.subset,\"Passed Biallelic Indels (MAF>0.05)\", rgb(0,0,1,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==0 & no_alleles==2)\n");
        fprintf(out, "plot_hist(data.subset,\"Failed Biallelic Indels\", rgb(1,0,0,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==0 & no_alleles==2 & maf>0.05)\n");
        fprintf(out, "plot_hist(data.subset,\"Failed Biallelic Indels (MAF>0.05)\", rgb(1,0,0,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==1 & no_alleles!=2)\n");
        fprintf(out, "plot_hist(data.subset,\"Passed Multiallelic Indels\", rgb(0,0,1,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==1 & no_alleles!=2 & maf>0.05)\n");
        fprintf(out, "plot_hist(data.subset,\"Passed Multiallelic Indels (Collapsed MAF>0.05)\", rgb(0,0,1,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==0 & no_alleles!=2)\n");
        fprintf(out, "plot_hist(data.subset,\"Failed Multiallelic Indels\", rgb(1,0,0,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, pass==0 & no_alleles!=2 & maf>0.05)\n");
        fprintf(out, "plot_hist(data.subset,\"Failed Multiallelic Indels (Collapsed MAF>0.05)\", rgb(1,0,0,0.5))\n");
        fprintf(out, "\n");
        fprintf(out, "dev.off()\n");

        fclose(out);

        std::string cmd = "cd "  + output_dir + "; cat plot.r | R --vanilla > run.log";
        int32_t ret = system(cmd.c_str());
    };

    void print_stats()
    {
        fprintf(stderr, "Stats \n");
        fprintf(stderr, "     no. of variants  : %d\n", no_variants);
        fprintf(stderr, "\n");
    };

    ~Igor()
    {
    };

    private:
};

}

void profile_hwe(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_hwe();
    igor.print_stats();
    igor.print_pdf();
}

