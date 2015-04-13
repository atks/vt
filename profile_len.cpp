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

#include "profile_len.h"

namespace
{
struct len_t
{
    int32_t dlen;
    float maf;
    int32_t no_alleles;
    float ab;
    int32_t coding;
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
    const char* AF;
    const char* AB;
    const char* GENCODE_NFS;
    const char* GENCODE_FS;

    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;

    ///////////////
    //general use//
    ///////////////
    std::vector<len_t> len_pts;

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
            std::string desc = "Profile length distribution";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_dir("x", "x", "output directory []", false, "plot_len", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output PDF file []", false, "len.pdf", "str", cmd);
            TCLAP::ValueArg<std::string> arg_AF("a", "af", "AF tag [AF]", false, "AF", "str", cmd);
            TCLAP::ValueArg<std::string> arg_AB("b", "ab", "AB tag [AB]", false, "AB", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            fexp = arg_fexp.getValue();
            output_dir = arg_output_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
            AF = strdup(arg_AF.getValue().c_str());
            AB = strdup(arg_AB.getValue().c_str());
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

    void profile_len()
    {
        bcf1_t *v = bcf_init1();

        Variant variant;
        float *af = (float*) malloc(sizeof(float));
        float *ab = (float*) malloc(sizeof(float));
        int32_t n_af=1, n_ab=1;

        while(odr->read(v))
        {
            bcf_unpack(v, BCF_UN_ALL);

            //bcf_print_liten(odr->hdr, v);

            if (bcf_get_n_allele(v)!=2)
            {
                continue;
            }

            vm->classify_variant(odr->hdr, v, variant);

            if (filter_exists)
            {
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }

            int32_t dlen = variant.alleles[0].dlen;
            float maf = -1;
            if (bcf_get_info_float(odr->hdr, v, AF, &af, &n_af)>0)
            {
                if (n_af==1)
                {
                    maf = af[0]>0.5? 1-af[0]: af[0];
                }
                else
                {
                    maf = 0;
                    for (size_t i=0; i<n_af; ++i)
                    {
                        maf += af[i];
                    }

                    maf = maf>0.5? 1-maf: maf;
                }
            }

            int32_t no_alleles = bcf_get_n_allele(v);

            if (bcf_get_info_float(odr->hdr, v, AB, &ab, &n_ab)<0)
            {
                ab[0] = -1;
            }


            int32_t fs = bcf_get_info_flag(odr->hdr,v, "FS", 0, 0);
            int32_t nfs = bcf_get_info_flag(odr->hdr,v, const_cast<char*>("NFS"), 0, 0);
            int32_t coding = fs==1 ? 1 : (nfs==1 ? 2 : -1);

            len_t h = {dlen, maf, no_alleles, ab[0], coding};
            len_pts.push_back(h);

            ++no_variants;
        }

        if (af) free(af);
        if (ab) free(ab);
    };

    void print_options()
    {
        std::clog << "profile_len v" << version << "\n\n";
        std::clog << "options:     input VCF file         " << input_vcf_file << "\n";
        std::clog << "         [a] AF tag                 " << AF << "\n";
        std::clog << "         [b] AB tag                 " << AB << "\n";
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
        fprintf(out, "dlen\tmaf\tno_alleles\tpass\tab\tcoding\n");
        for (size_t i=0; i<len_pts.size(); ++i)
        {
            fprintf(out, "%d\t%f\t%d\t%d\t%f\t%d\n", len_pts[i].dlen, len_pts[i].maf, len_pts[i].no_alleles,
                                                  1, len_pts[i].ab, len_pts[i].coding);
        }
        fclose(out);

        //create r script
        file_path = output_dir + "/plot.r";
        out = fopen(file_path.c_str(), "w");

        fprintf(out, "setwd(\"%s\")\n", output_dir.c_str());
        fprintf(out, "\n");
        fprintf(out, "data = read.table(\"data.txt\", header=T)\n");
        fprintf(out, "\n");
        fprintf(out, "pdf(\"%s\",7,5)\n", output_pdf_file.c_str());
        fprintf(out, "\n");
        fprintf(out, "layout(matrix(c(1,2), 2, 1, byrow = TRUE))\n");
        fprintf(out, "data.subset = subset(data, abs(dlen)<=10 & pass==1 & no_alleles==2)\n");
        fprintf(out, "hist(data.subset$dlen, breaks=seq(-10,10,1),col=rgb(0,0,1,0.5), xlab=\"length\", main=\"Passed Indel Length Distribution\")\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, abs(dlen)<=10 & pass!=1 & no_alleles==2)\n");
        fprintf(out, "hist(data.subset$dlen, breaks=seq(-10,10,1),col=rgb(1,0,0,0.5), xlab=\"length\", main=\"Failed Indel Length Distribution\")\n");
        fprintf(out, "\n");
        fprintf(out, "layout(matrix(c(1,2), 2, 1, byrow = TRUE))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, abs(dlen)<=10 & pass==1 & no_alleles==2 & coding>0)\n");
        fprintf(out, "hist(data.subset$dlen, breaks=seq(-10,10,1),col=rgb(0,0,1,0.5), xlab=\"length\", main=\"Passed Coding Indel Length Distribution\")\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, abs(dlen)<=10 & pass!=1 & no_alleles==2 & coding>0)\n");
        fprintf(out, "hist(data.subset$dlen, breaks=seq(-10,10,1),col=rgb(1,0,0,0.5), xlab=\"length\", main=\"Failed Coding Indel Length Distribution\")\n");
        fprintf(out, "\n");
        fprintf(out, "layout(matrix(c(1,2), 1, 2, byrow = TRUE))\n");
        fprintf(out, "data.subset = subset(data, dlen>0 & pass!=1 & no_alleles==2)\n");
        fprintf(out, "x = table(abs(data.subset$dlen))/nrow(data.subset)\n");
        fprintf(out, "pie(x, main=\"Passed Indel Length Proportion\")\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, dlen<0 & pass!=1 & no_alleles==2)\n");
        fprintf(out, "x = table(abs(data.subset$dlen))/nrow(data.subset)\n");
        fprintf(out, "pie(x, main=\"Failed Indel Length Proportion\")\n");
        fprintf(out, "\n");
        fprintf(out, "layout(matrix(c(1,1), 2, 1, byrow = TRUE))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, abs(dlen)<=30 & pass==1 & no_alleles==2 & ab<1 & ab>0)\n");
        fprintf(out, "boxplot(data.subset$ab~data.subset$dlen, col=rgb(0,0,1,0.5), pch=\".\", main=\"Passed Indel Allele Balance Profile\")\n");
        fprintf(out, "abline(h=0.5, col=\"grey\")\n");
        fprintf(out, "\n");
        fprintf(out, "layout(matrix(c(1,1), 2, 1, byrow = TRUE))\n");
        fprintf(out, "\n");
        fprintf(out, "data.subset = subset(data, abs(dlen)<=30 & pass!=1 & no_alleles==2 & ab<1 & ab>0)\n");
        fprintf(out, "boxplot(data.subset$ab~data.subset$dlen, col=rgb(1,0,0,0.5), pch=\".\", main=\"Failed Indel Allele Balance Profile\")\n");
        fprintf(out, "abline(h=0.5, col=\"grey\")\n");
        fprintf(out, "\n");
        fprintf(out, "dev.off()\n");

        fclose(out);

        //run script
        std::string cmd = "cd "  + output_dir + "; cat plot.r | R --vanilla > run.log";
        //int32_t ret = system(cmd.c_str());
    };

    void print_stats()
    {
        fprintf(stderr, "Stats \n");
        fprintf(stderr, "     no. of variants  : %d\n", no_variants);
        fprintf(stderr, "\n");

        //do a textual histogram print out of len spectrum
        fprintf(stderr, "dlen\tmaf\tno_alleles\tpass\tab\tcoding\n");
        for (size_t i=0; i<len_pts.size(); ++i)
        {
            fprintf(stderr, "%d\t%f\t%d\t%d\t%f\t%d\n", len_pts[i].dlen, len_pts[i].maf, len_pts[i].no_alleles,
                                                  1, len_pts[i].ab, len_pts[i].coding);
        }
        fclose(stderr);

    };

    ~Igor()
    {
        odr->close();
    };

    private:
};

}

void profile_len(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_len();
    igor.print_stats();
    igor.print_pdf();
}

