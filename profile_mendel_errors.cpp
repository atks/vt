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

#include "profile_mendel_errors.h"

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
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    int32_t min_depth;
    float_t min_gq;
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    
    ///////////////
    //general use//
    ///////////////
    kstring_t variant;

    /////////
    //stats//
    /////////
    uint32_t no_trios;
    uint32_t no_variants;
    
    uint32_t no_00;
    uint32_t no_00_00;
    uint32_t no_00_01;
    uint32_t no_00_11;
    
    uint32_t no_11;
    uint32_t no_11_00;
    uint32_t no_11_01;
    uint32_t no_11_11;

    /////////
    //tools//
    /////////
    Pedigree *pedigree;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Profile Mendel Errors.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_ped_file("p", "p", "pedigree file", true, "", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_min_depth("d", "d", "minimum depth", false, 5, "str", cmd);
            TCLAP::ValueArg<float> arg_min_gq("q", "q", "minimum genotype quality", false, 2, "str", cmd);
            
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            input_ped_file = arg_input_ped_file.getValue();
            min_depth = arg_min_depth.getValue();
            min_gq = arg_min_gq.getValue();
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

        ///////////////////////////
        //ped file initialization//
        ///////////////////////////
        pedigree = new Pedigree(input_ped_file);

        ///////////////
        //general use//
        ///////////////
        variant = {0,0,0};
        
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_trios = 0;
        no_variants = 0;
        
        no_00 = 0;
        no_00_00 = 0;
        no_00_01 = 0;
        no_00_11 = 0;
        
        no_11 = 0;
        no_11_00 = 0;
        no_11_01 = 0;
        no_11_11 = 0;
        
        /////////
        //tools//
        /////////
    }

    void profile_mendel_errors()
    {
        bcf_hdr_t *h = odr->hdr;
        bcf1_t *v = bcf_init1();
        
        std::vector<Trio>& trios = pedigree->trios;
        no_trios = trios.size();        

        int32_t missing = 0;
        int32_t mendel_homalt_err = 0;
            
        while(odr->read(v))
        {
            int8_t *gts = NULL;
            int32_t n = 0;
            int k = bcf_get_genotypes(h, v, &gts, &n); 
                
            int32_t *dps = NULL;
            n = 0;
            k = bcf_get_format_int(h, v, "DP", &dps, &n); 
         
            int32_t *pls = NULL;
            n = 0;
            k = bcf_get_format_int(h, v, "PL", &pls, &n); 
              
            bool variant_used = false;
            
            for (int32_t i =0; i< trios.size(); ++i)
            {
                int32_t j = bcf_hdr_id2int(h, BCF_DT_SAMPLE, trios[i].father.c_str());
                int8_t *igt = gts+j*2;
                int32_t f1 = bcf_gt_allele(igt[0]);
                int32_t f2 = bcf_gt_allele(igt[1]);
                int32_t fdp = dps[j];
                
                j = bcf_hdr_id2int(h, BCF_DT_SAMPLE, trios[i].mother.c_str());
                igt = gts+j*2;
                int32_t m1 = bcf_gt_allele(igt[0]);
                int32_t m2 = bcf_gt_allele(igt[1]);
                int32_t mdp = dps[j];

                j = bcf_hdr_id2int(h, BCF_DT_SAMPLE, trios[i].child.c_str());
                igt = gts+j*2;
                int32_t c1 = bcf_gt_allele(igt[0]);
                int32_t c2 = bcf_gt_allele(igt[1]);
                int32_t cdp = dps[j];

                std::cout << f1 << "/" << f2 << ":";
                std::cout << m1 << "/" << m2 << ":";
                std::cout << c1 << "/" << c2 << "\n";

                if (f1==1 && f2==0)
                {
                    std::cerr << f1 << "/" << f2 << "\t";
                    std::cerr << m1 << "/" << m2 << "\t";
                    std::cerr << c1 << "/" << c2 << "\n";
                }

                if (f1<0 && f2<0 & m1<0 && m2<0 && c1<0 && c2<0)
                {
                    ++missing;
                }
                
                if (f1==0 && f2==1 & m1==0 && m2==1 && c1==1 && c2==1)
                {
                    ++mendel_homalt_err;
                }
                
//                if (f1==0 && f2==1 & m1==1 && m2==1 && c1==0 && c2==0)
//                {
//                    ++mendel_homalt_err;
//                }
                
                if (!(f1<0 || f2<0 || m1<0 || m2<0 || c1<0 || c2<0) && 
                    (fdp>=min_depth && mdp>=min_depth && cdp>=min_depth) )
                {
                    if ((f1==0 && f2==1 && m1==0 && m2==0) ||
                        (f1==0 && f2==0 && m1==0 && m2==1))
                    {
                        ++no_00;
                        
                        if (c1==0 && c2==0)
                        {
                            ++no_00_00;        
                        }
                        else if (c1==0 && c2==1)
                        {
                            ++no_00_01;        
                        }
                        else if (c1==1 && c2==1)
                        {
                            ++no_00_11;        
                        }
                        else
                        {
                            std::cerr << f1 << "/" << f2 << "\t";
                            std::cerr << m1 << "/" << m2 << "\t";
                            std::cerr << c1 << "/" << c2 << "\n";
                        }
                    }
                    else if ((f1==0 && f2==1 && m1==1 && m2==1) ||
                             (f1==1 && f2==1 && m1==0 && m2==1))
                    {
                        ++no_11;
                        
                        if (c1==0 && c2==0)
                        {
                            ++no_11_00;        
                        }
                        else if (c1==0 && c2==1)
                        {
                            ++no_11_01;        
                        }
                        else if (c1==1 && c2==1)
                        {
                            ++no_11_11;        
                        }
                        else
                        {
                            std::cerr << f1 << "/" << f2 << "\t";
                            std::cerr << m1 << "/" << m2 << "\t";
                            std::cerr << c1 << "/" << c2 << "\n";
                        }
                    } 
                    
                    variant_used = true;
                }
            }

            

            ++no_variants;    
            free(gts);
            free(dps);
            free(pls);
        }

    std::cerr << "MISSING " << missing << "\n";
    std::cerr << "HOMALT MENDEL ERROR " << mendel_homalt_err << "\n";
        odr->close();
    };

    void print_options()
    {
        std::clog << "merge_candidate_variants v" << version << "\n\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [p] input PED file        " << input_ped_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "  Mendelian Errors\n");
        fprintf(stderr, "              0/0     0/1    1/1       Total\n");
        fprintf(stderr, "  0/0+0/1   %5.2f   %5.2f  %5.2f  %10d\n", (float)no_00_00/no_00*100, (float)no_00_01/no_00*100, (float)no_00_11/no_00*100, no_00);
        fprintf(stderr, "  1/1+0/1   %5.2f   %5.2f  %5.2f  %10d\n", (float)no_11_00/no_11*100, (float)no_11_01/no_11*100, (float)no_11_11/no_11*100, no_11);
        fprintf(stderr, "\n");
        fprintf(stderr, "  no. of trios     : %d\n", no_trios);
        fprintf(stderr, "  no. of variants  : %d\n", no_variants);
        fprintf(stderr, "\n");
        fprintf(stderr, "  minimum depth    : %d\n", min_depth);
        fprintf(stderr, "\n");
    };

    ~Igor()
    {
    };

    private:
};

}

void profile_mendel_errors(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_mendel_errors();
    igor.print_stats();
}

