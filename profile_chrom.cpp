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

#include "profile_chrom.h"

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
    khash_t(32) *chrom;
    khash_t(32) *pass_chrom;

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
            std::string desc = "Plot Chromosome Distribution.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_dir("x", "x", "output directory []", false, "plot_chrom", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_pdf_file("y", "y", "output PDF file []", false, "chrom.pdf", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            fexp = arg_fexp.getValue();
            output_dir = arg_output_dir.getValue();
            output_pdf_file = arg_output_pdf_file.getValue();
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

        chrom = kh_init(32);
        pass_chrom = kh_init(32);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        /////////
        //tools//
        /////////
        vm = new VariantManip();
    }

    void profile_chrom()
    {
        bcf1_t *v = bcf_init1();

        Variant variant;
        
        while(odr->read(v))
        {
            bcf_unpack(v, BCF_UN_ALL);
            
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
            
            int32_t rid = bcf_get_rid(v);
            
            k = kh_get(32, chrom, rid);
            if (k==kh_end(chrom))
            {
                k = kh_put(32, chrom, rid, &ret);
                kh_value(chrom, k) = 1;
            }
            else
            {
                kh_value(chrom, k) = kh_value(chrom, k) + 1;
            }    
            
            if (pass)
            {
                k = kh_get(32, pass_chrom, rid);
                if (k==kh_end(pass_chrom))
                {
                    k = kh_put(32, pass_chrom, rid, &ret);
                    kh_value(pass_chrom, k) = 1;
                }
                else
                {
                    kh_value(pass_chrom, k) = kh_value(pass_chrom, k) + 1;
                } 
            }    
            
            ++no_variants;
        }
       
        
    };

    void print_options()
    {
        std::clog << "plot_chrom v" << version << "\n\n";
        std::clog << "options:     input VCF file         " << input_vcf_file << "\n";
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
        
        fprintf(out, "chromosome\ttotal\tpass\tlength\n");
        
        int32_t nseqs;
        const char** seq_names = bcf_hdr_seqnames(odr->hdr, &nseqs);
        int32_t* len = bcf_hdr_seqlen(odr->hdr, &nseqs);
        
        for (int32_t i=0; i<nseqs; ++i)
        {
            int32_t no_total = 0, no_pass = 0;
            k = kh_get(32, chrom, i);
            if (k!=kh_end(chrom))
            {
                no_total = kh_value(chrom, k);
            } 
            
            k = kh_get(32, pass_chrom, i);
            if (k!=kh_end(pass_chrom))
            {
                no_pass = kh_value(pass_chrom, k);
            } 
            
            fprintf(out, "%s\t%d\t%d\t%d\n", seq_names[i], no_total, no_pass, len[i]);
        }
        if(nseqs) free(seq_names);   
        if(nseqs) free(len); 
        fprintf(out, "\n");        
        fclose(out);

        //create r script
        file_path = output_dir + "/plot.r"; 
        out = fopen(file_path.c_str(), "w");
        
        fprintf(out, "setwd(\"%s\")\n", output_dir.c_str());
        fprintf(out, "\n");
        fprintf(out, "data = read.table(\"data.txt\", header=T)\n");
        fprintf(out, "data.subset = subset(data, total!=0)\n");
        fprintf(out, "pdf(\"%s\",9,5)\n", output_pdf_file.c_str());
        fprintf(out, "par(mar = c(5, 4, 4, 5))\n");
        fprintf(out, "plot(c(1:nrow(data.subset)), data.subset$total, type=\"h\", lwd=10, col=\"grey\", xaxt=\"n\", ylab=\"No. of variants\", xlab=\"Chromosome\")\n");
        fprintf(out, "axis(side=1, at=c(1:nrow(data.subset)), labels=data.subset$chromosome)\n");
        fprintf(out, "par(new = TRUE)\n");
        fprintf(out, "lines(data.subset$pass, type=\"h\", lwd=10, col=rgb(0,0,1,0.5)) \n");
        fprintf(out, "par(new = TRUE)\n");
        fprintf(out, "plot(data.subset$length, type=\"h\", lwd=2, col=\"red\", axes=FALSE, bty = \"n\", xlab = \"\", ylab = \"\") \n");
        fprintf(out, "mtext(\"Length\", side=4, line=3)\n");
        fprintf(out, "legend(\"topright\", c(\"all\", \"pass\", \"chromosome length\"), col = c(\"grey\", rgb(0,0,1,0.5), \"red\"), pch = 20)\n");
        fprintf(out, "dev.off()\n");
        fprintf(out, "\n");
        
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
        
        //do a textual histogram print out of chrom
        int32_t nseqs;
        const char** seq_names = bcf_hdr_seqnames(odr->hdr, &nseqs);
        
        
        int32_t* len = bcf_hdr_seqlen(odr->hdr, &nseqs);
        
        fprintf(stderr, "     chromosome     total#      pass#      length\n");
        for (int32_t i=0; i<nseqs; ++i)
        {
            int32_t no_total = 0, no_pass = 0;
            k = kh_get(32, chrom, i);
            if (k!=kh_end(chrom))
            {
                no_total = kh_value(chrom, k);
            } 
            
            k = kh_get(32, pass_chrom, i);
            if (k!=kh_end(pass_chrom))
            {
                no_pass = kh_value(pass_chrom, k);
            } 
            
            fprintf(stderr, "%15s   %8d   %8d   %9d\n", seq_names[i], no_total, no_pass, len[i]);
        }
        fprintf(stderr, "\n");
        if(nseqs) free(seq_names);   
        if(nseqs) free(len);    
    };

    ~Igor()
    {
        odr->close();
        kh_destroy(32, chrom);
        kh_destroy(32, pass_chrom);
    };

    private:
};
}

void profile_chrom(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_chrom();
    igor.print_stats();
    igor.print_pdf();
}

