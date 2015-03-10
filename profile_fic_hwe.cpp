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

#include "profile_fic_hwe.h"

namespace
{

KHASH_MAP_INIT_INT(32, int32_t)

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
    const char* AC;
    const char* AN;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;

    ///////////////
    //general use//
    ///////////////
    int ret, is_missing;
    khiter_t k;
    khash_t(32) *afs;
    khash_t(32) *pass_afs;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    uint32_t no_variants;

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
            TCLAP::ValueArg<std::string> arg_output_dir("x", "x", "output directory []", false, "plot_afs", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output PDF file []", false, "afs.pdf", "str", cmd);
            TCLAP::ValueArg<std::string> arg_AC("c", "ac", "AC tag [AC]", false, "AC", "str", cmd);
            TCLAP::ValueArg<std::string> arg_AN("n", "an", "AN tag [AN]", false, "AN", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            fexp = arg_fexp.getValue();
            output_dir = arg_output_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
            AC = strdup(arg_AC.getValue().c_str());
            AN = strdup(arg_AN.getValue().c_str());
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

        afs = kh_init(32);
        pass_afs = kh_init(32);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        /////////
        //tools//
        /////////
        vm = new VariantManip();
    }

    void profile_afs()
    {
        bcf1_t *v = bcf_init1();

        Variant variant;
        int32_t *ac=NULL, *an=NULL, n_ac=0, n_an=0;

        while(odr->read(v))
        {
            bcf_unpack(v, BCF_UN_ALL);

            bcf_print(odr->hdr, v);

            if (bcf_get_n_allele(v)!=2)
            {
                continue;
            }

            if (filter_exists)
            {
                int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }

            bool pass = (bcf_has_filter(odr->hdr, v, const_cast<char*>("PASS"))==1);

            bcf_get_info_int32(odr->hdr, v, AC, &ac, &n_ac);
            bcf_get_info_int32(odr->hdr, v, AN, &an, &n_an);

            if (ac[0]>(an[0]>>1)) {ac[0] = an[0]-ac[0];}

            if (ac[0]==0) continue;

            k = kh_get(32, afs, ac[0]);
            if (k==kh_end(afs))
            {
                k = kh_put(32, afs, ac[0], &ret);
                kh_value(afs, k) = 1;
            }
            else
            {
                kh_value(afs, k) = kh_value(afs, k) + 1;
            }

            if (pass)
            {
                k = kh_get(32, pass_afs, ac[0]);
                if (k==kh_end(pass_afs))
                {
                    k = kh_put(32, pass_afs, ac[0], &ret);
                    kh_value(pass_afs, k) = 1;
                }
                else
                {
                    kh_value(pass_afs, k) = kh_value(pass_afs, k) + 1;
                }
            }

            ++no_variants;
        }

        if (n_ac) free(ac);
        if (n_an) free(an);
    };

    void print_options()
    {
        std::clog << "plot_afs v" << version << "\n\n";
        std::clog << "options:     input VCF file         " << input_vcf_file << "\n";
        std::clog << "         [c] AC tag                 " << AC << "\n";
        std::clog << "         [n] AN tag                 " << AN << "\n";
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

        fprintf(out, "mac\tf\tpass\n");
        for (k=kh_begin(afs); k!=kh_end(afs); ++k)
        {
            if (kh_exist(afs, k))
            {
                fprintf(out, "%d\t%d\t0\n", k , kh_value(afs, k));
            }
        }

        for (k=kh_begin(pass_afs); k!=kh_end(pass_afs); ++k)
        {
            if (kh_exist(pass_afs, k))
            {
                fprintf(out, "%d\t%d\t1\n", k , kh_value(pass_afs, k));
            }
        }

        fclose(out);

        //create r script
        file_path = output_dir + "/plot.r";
        out = fopen(file_path.c_str(), "w");

        fprintf(out, "setwd(\"%s\")\n", output_dir.c_str());
        fprintf(out, "\n");
        fprintf(out, "data = read.table(\"data.txt\", header=T)\n");
        fprintf(out, "data.all=subset(data, pass==0)\n");
        fprintf(out, "pdf(\"%s\",7,5)\n", output_pdf_file.c_str());
        fprintf(out, "plot(data.all$mac, data.all$f,  log=\"xy\", pch=20, cex=0.5, col=rgb(1,0,0,0.5), main=\"Folded Allele Frequency Spectrum\", xlab=\"Minor allele counts\", ylab=\"Frequency\", panel.last=grid(equilogs=FALSE, col=\"grey\"))\n");
        fprintf(out, "data.pass=subset(data, pass==1)\n");
        fprintf(out, "points(data.pass$mac, data.pass$f, pch=20, cex=0.5, col=rgb(0,0,1,0.5))\n");
        fprintf(out, "legend(\"topright\", c(\"pass\", \"all\"), col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), pch = 20)\n");
        fprintf(out, "dev.off()\n");

        fclose(out);

        //run script
        std::string cmd = "cd "  + output_dir + "; cat plot.r | R --vanilla > run.log";
        int32_t ret = system(cmd.c_str());
    };

    void print_stats()
    {
        fprintf(stderr, "Stats \n");
        fprintf(stderr, "     no. of variants  : %d\n", no_variants);
        fprintf(stderr, "\n");

        //do a textual histogram print out of afs


    };

    ~Igor()
    {
        odr->close();
        kh_destroy(32, afs);
        kh_destroy(32, pass_afs);
    };

    private:
};

}

void profile_fic_hwe(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_afs();
    igor.print_stats();
    igor.print_pdf();
}

