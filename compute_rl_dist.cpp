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

#include "compute_rl_dist.h"

namespace
{
class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string ref_fasta_file;
    std::vector<std::string> x;
    std::string chrom;
    ReferenceSequence *rs;    
        

    bool debug;
    uint32_t no;
  
    Igor(int argc, char **argv)
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
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_chromosome("c", "c", "chromosome fasta file []", true, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_x("ap", "#ploidy #alleles", true, "", cmd);

            cmd.parse(argc, argv);

            ref_fasta_file = arg_ref_fasta_file.getValue();
            x = arg_x.getValue();
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


        ////////////////////////
        //tools initialization//
        ////////////////////////
        rs = new ReferenceSequence(ref_fasta_file);
    }

    /**
     * Gets number of genotypes from number of alleles and ploidy.
     */
    uint32_t compute_rl_dist()
    {
        
         

              return 1;
        
    }

    void print_options()
    {
//        if (!print) return;
//
//        std::clog << "normalize v" << version << "\n";
//        std::clog << "\n";
//        std::clog << "options:     input VCF file                                  " << input_vcf_file << "\n";
//        std::clog << "         [o] output VCF file                                 " << output_vcf_file << "\n";
//        std::clog << "         [w] sorting window size                             " << window_size << "\n";
//        print_str_op("         [f] filter                                          ", fexp);
//        std::clog << "         [n] no fail on reference inconsistency for non SNPs " << (!strict ? "true" : "false") << "\n";
//        std::clog << "         [q] quiet                                           " << (!print ? "true" : "false") << "\n";
//        std::clog << "         [d] debug                                           " << (debug ? "true" : "false")  << "\n";
//        std::clog << "         [r] reference FASTA file                            " << ref_fasta_file << "\n";
//        print_int_op("         [i] intervals                                       ", intervals);
//        std::clog << "\n";
    }

    void print_stats()
    {
//        if (!print) return;
//
//        int32_t no_biallelic_normalized = no_lt+no_rt+no_lt_rt+no_rt_la+no_la;
//        int32_t no_multiallelic_normalized = no_multi_lt+no_multi_rt+no_multi_lt_rt+no_multi_rt_la+no_multi_la;
//        int32_t no_normalized = no_biallelic_normalized + no_multiallelic_normalized;
//
//        std::clog << "\n";
//        std::clog << "stats: biallelic\n";
//        std::clog << "          no. left trimmed                      : " << no_lt << "\n";
//        std::clog << "          no. right trimmed                     : " << no_rt << "\n";
//        std::clog << "          no. left and right trimmed            : " << no_lt_rt << "\n";
//        std::clog << "          no. right trimmed and left aligned    : " << no_rt_la << "\n";
//        std::clog << "          no. left aligned                      : " << no_la << "\n";
//        std::clog << "\n";
//        std::clog << "       total no. biallelic normalized           : " << no_biallelic_normalized << "\n";
//        std::clog << "\n";
//        std::clog << "       multiallelic\n";
//        std::clog << "          no. left trimmed                      : " << no_multi_lt << "\n";
//        std::clog << "          no. right trimmed                     : " << no_multi_rt << "\n";
//        std::clog << "          no. left and right trimmed            : " << no_multi_lt_rt << "\n";
//        std::clog << "          no. right trimmed and left aligned    : " << no_multi_rt_la << "\n";
//        std::clog << "          no. left aligned                      : " << no_multi_la << "\n";
//        std::clog << "\n";
//        std::clog << "       total no. multiallelic normalized        : " << no_multiallelic_normalized << "\n";
//        std::clog << "\n";
//        std::clog << "       total no. variants normalized            : " << no_normalized << "\n";
//        std::clog << "       total no. variants observed              : " << no_variants << "\n";
//        std::clog << "       total no. reference observed             : " << no_refs << "\n";
//        std::clog << "\n";
    };
    
    ~Igor() {};

    private:
};

}

void compute_rl_dist(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.compute_rl_dist();
    igor.print_stats();
};
