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

#include "align.h"
#include "tbx_ordered_reader.h"
#include "bed.h"
#include "pcre2.h"
#include "pregex.h"

namespace
{
    
void print_time(float t)
{
    if (t<60)
    {
        fprintf(stderr, "Alignment time elapsed: %fs\n\n", t);
    }
    else if (t<60*60) //less than an hour
    {
        fprintf(stderr, "Time elapsed: %dm %ds\n\n", ((int32_t)(t/60)), ((int32_t)fmod(t, 60)));
    }
    else if (t<60*60*24) //less than a day
    {
        double m = fmod(t, 60*60); //remaining minutes
        fprintf(stderr, "Time elapsed: %dh %dm %ds\n\n", ((int32_t)(t/(60*60))), ((int32_t)(m/60)), ((int32_t)fmod(m, 60)));
    }
    else if (t<60*60*24*365) //less than a year
    {
        double h = fmod(t, 60*60*24); //remaining hours
        double m = fmod(h, 60*60); //remaining minutes
        fprintf(stderr, "Alignment time elapsed: %dd %dh %dm %.6fs\n\n", ((int32_t)(t/(60*60*24))), ((int32_t)(h/(60*60))), ((int32_t)(m/60)), (fmod(m, 60)));
    }
};

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string method;
    std::vector<std::string> x;


    bool debug;
    uint32_t no;


  
    Igor(int argc, char **argv)
    {
        
    };

    /**
     * n choose r.
     */
    uint32_t choose(uint32_t n, uint32_t r)
    {
        if (r>n)
        {
            return 0;
        }
        else if (r==n)
        {
            return 1;
        }
        else if (r==0)
        {
            return 1;
        }
        else
        {
            if (r>(n>>1))
            {
                r = n-r;
            }

            uint32_t num = n;
            uint32_t denum = 1;

            for (uint32_t i=1; i<r; ++i)
            {
                num *= n-i;
                denum *= i+1;
            }

            return num/denum;
        }
    }

    /**
     * Gets number of genotypes from number of alleles and ploidy.
     */
    uint32_t bcf_ap2g(uint32_t no_allele, uint32_t no_ploidy)
    {
        return choose(no_ploidy+no_allele-1, no_allele-1);
    }

    /**
     * Gets number of genotypes from number of alleles and ploidy.
     */
    uint32_t bcf_g2i(std::string genotype)
    {
        uint32_t index = 0;
        for (uint32_t i=0; i<genotype.size(); ++i)
        {
            uint32_t allele = genotype.at(i)-65;
            index += bcf_ap2g(allele, i+1);
        }
        return index;
    }

    void print_genotypes(uint32_t A, uint32_t P, std::string genotype)
    {
        if (genotype.size()==P)
        {
            std::cerr << no << ") " << genotype << " " << bcf_g2i(genotype) << "\n";
            ++no;
        }
        else
        {
            for (uint32_t a=0; a<A; ++a)
            {
                std::string s(1,(char)(a+65));
                s.append(genotype);
                print_genotypes(a+1, P, s);
            }
        }
    }

    void test(int argc, char ** argv)
    {
        
        
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "vt test  -m detect_motif -s ACTGACT \n";

            std::string version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::UnlabeledMultiArg<std::string> arg_x("ap", "#ploidy #alleles", true, "", cmd);

            cmd.parse(argc, argv);

            x = arg_x.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }

        no = 1;
        
//        print_genotypes(x[0], x[1], "");
//        uint32_t g = bcf_ap2g(x[0], x[1]);
//
//        std::cerr << "A: " << x[0] << " P: " << x[1] << " G: " << g << "\n";

        vcfFile *vcf = bcf_open(x[0].c_str(), "rb");
        bcf_hdr_t *h = bcf_hdr_read(vcf);
        bcf1_t *v = bcf_init();
        
        
        std::cerr << "writing to " << x[1] << "\n";
        vcfFile *ovcf = bcf_open(x[1].c_str(), "wu");
        bcf_hdr_write(ovcf, h);
        
        while (bcf_read(vcf, h, v)>=0)
        {
            //std::cerr << "test\n";
            bcf_write(ovcf, h, v);    
        }
        
        bcf_close(ovcf);
        bcf_close(vcf);
        
    };

    void print_stats()
    {
        std::clog << "\n";

    };

    /**
     * Get a sequence.  User have to free the char* returned.
     */
    char* get_sequence(std::string& chrom, uint32_t pos1, uint32_t len, faidx_t* fai)
    {
        int ref_len = 0;
        
        
        char* seq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos1-1, pos1+len-2, &ref_len);
        
        
        if (!seq || ref_len!=len)
        {
            fprintf(stderr, "[%s:%d %s] failure to extract sequence from fasta file: %s:%d: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos1-1);
            exit(1);
        }
    
        return seq;
    };

    /**
     *  Analyse an mdust file
     *  Compute how much of genome it is off.
     */    
    void analyse_mdust(int argc, char ** argv)
    {
        std::string bed_file;
        std::string ref_fasta_file;
        faidx_t *fai;
        
        try
        {
            std::string desc = "analyse_mdust  -m detect_motif -s ACTGACT \n";

            std::string version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_bed_file("b", "b", "BED file", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference FASTA file", true, "", "string", cmd);

            cmd.parse(argc, argv);

            bed_file = arg_bed_file.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }

        TBXOrderedReader todr(bed_file);

        if (ref_fasta_file!="")
        {
            fai = fai_load(ref_fasta_file.c_str());
            if (fai==NULL)
            {
                fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
                exit(1);
            }
        }
        
        kstring_t s = {0,0,0};
        
        uint32_t no_non_n = 0;
        uint32_t no = 0;
        
        std::string chrom = "";
        
        while (todr.read(&s))
        {
            BEDRecord br(&s);
            //std::cerr << br.chrom << " " << br.end1 << "-" << br.start1 << " (" << (br.end1-br.start1+1)<< ")\n";
        
            if (br.chrom!=chrom)
            {
                chrom = br.chrom;
                if (chrom=="GL000207.1")
                {
                    break;
                }
                std::cerr << chrom << "\n";
            }    
        
            int ref_len = 0;
            char* seq = faidx_fetch_uc_seq(fai, br.chrom.c_str(), br.start1-1, br.end1-1, &ref_len);
            //std::cerr << seq << " (" << ref_len <<")\n";
        
            for (uint32_t i=0; i<ref_len; ++i)
            {
                if (seq[i]!='N')
                {
                    ++no_non_n;
                }   
                
                ++no; 
            }
        
            free(seq);
        }


        std::cerr << no_non_n << "/" << no << "\n";


    };
        
    /**
     *  Analyse an mdust file
     *  Compute how much of genome it is off.
     */    
    void analyse_str(int argc, char ** argv)
    {
        std::string bed_file;
        std::string ref_fasta_file;
        faidx_t *fai;
        
        try
        {
            std::string desc = "analyse_mdust  -m detect_motif -s ACTGACT \n";

            std::string version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref("v", "v", "reference", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> arg_alt("a", "a", "alternative", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference FASTA file", true, "", "string", cmd);

            cmd.parse(argc, argv);

            ref_fasta_file = arg_ref_fasta_file.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }

        TBXOrderedReader todr(bed_file);

        if (ref_fasta_file!="")
        {
            fai = fai_load(ref_fasta_file.c_str());
            if (fai==NULL)
            {
                fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
                exit(1);
            }
        }
        
        


    };

   /**
     *  Analyse an mdust file
     *  Compute how much of genome it is off.
     */    
    void test_pcre2(int argc, char ** argv)
    {
        std::string regex;
        std::string text;
        faidx_t *fai;
        
        try
        {
            std::string desc = "analyse_mdust  -m detect_motif -s ACTGACT \n";

            std::string version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_regex("r", "r", "regular expression", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> arg_text("t", "t", "text", true, "", "string", cmd);

            cmd.parse(argc, argv);

            regex = arg_regex.getValue();
            text = arg_text.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }

        printf("regex = %s\n", regex.c_str());
        printf("text  = %s\n", text.c_str());
        
        PERLregex pregex;
        
        pregex.set(regex);
        std::cerr << "matched: " << pregex.match(text) << "\n";
        
        
//        pcre2_code *re;
//        PCRE2_SPTR pattern;     /* PCRE2_SPTR is a pointer to unsigned code units of */
//        PCRE2_SPTR subject;     /* the appropriate width (8, 16, or 32 bits). */
//        PCRE2_SPTR name_table;
//        
//        int crlf_is_newline;
//        int errornumber;
//        int find_all;
//        int i;
//        int namecount;
//        int name_entry_size;
//        int rc;
//        int utf8;
//        
//        uint32_t option_bits;
//        uint32_t newline;
//        
//        PCRE2_SIZE erroroffset;
//        PCRE2_SIZE *ovector;
//        
//        size_t subject_length;
//        pcre2_match_data *match_data;
//
//        pattern = (PCRE2_SPTR)regex.c_str();
//        subject = (PCRE2_SPTR)text.c_str();
//        subject_length = strlen((char *)subject);
//
//        re = pcre2_compile(
//        pattern,               /* the pattern */
//        PCRE2_ZERO_TERMINATED, /* indicates pattern is zero-terminated */
//        0,                     /* default options */
//        &errornumber,          /* for error number */
//        &erroroffset,          /* for error offset */
//        NULL);  
//        
//        if (re == NULL)
//        {
//            PCRE2_UCHAR buffer[256];
//            pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
//            printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroroffset,
//            buffer);
//        
//        }
//                
//        match_data = pcre2_match_data_create_from_pattern(re, NULL);
//        
//        rc = pcre2_match(
//                    re,                   /* the compiled pattern */
//                    subject,              /* the subject string */
//                    subject_length,       /* the length of the subject */
//                    0,                    /* start at offset 0 in the subject */
//                    0,                    /* default options */
//                    match_data,           /* block for storing the result */
//                    NULL);                /* use default match context */
//      
//        if (rc < 0)
//        {
//          switch(rc)
//            {
//            case PCRE2_ERROR_NOMATCH: printf("No match\n"); break;
//            /*
//            Handle other special cases if you like
//            */
//            default: printf("Matching error %d\n", rc); break;
//            }
//          pcre2_match_data_free(match_data);   /* Release memory used for the match */
//          pcre2_code_free(re);                 /* data and the compiled pattern. */
//        }
//          
//        if (rc>0)
//        {
//            printf("matched!\n");
//        }  
          
    };
    

    ~Igor() {};

    private:
};

}

void test(int argc, char ** argv)
{
    Igor igor(argc, argv);
    
    //igor.analyse_mdust(argc, argv);
    
    //igor.test_pcre2(argc, argv);
    
    //igor.test();
};
