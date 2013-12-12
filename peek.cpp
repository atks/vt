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

    /////////
    //stats//
    /////////
    uint32_t no_snp2;
    uint32_t no_snp3;
    uint32_t no_snp4;

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
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
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

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_snp2 = 0;
        no_snp3 = 0;
        no_snp4 = 0;
//        no_mnp2 = 0;
//        no_mnp3+no_mnp4+no_mnp5
        
//        no_mnps << "\n";
//        no_mnps_bi << "\n";
//        std::clog << "            multiallelic      : " << no_mnps_multi << "\n";
//        std::clog << "             length 2         : " << no_mnps2 << "\n";
//        std::clog << "             length 3         : " << no_mnps3 << "\n";
//        std::clog << "             length 4         : " << no_mnps4 << "\n";
//        std::clog << "             length >=5       : " << no_mnps5 << "\n";
        

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip(ref_fasta_file);
    }

    void peek()
    {
        while (odr->read(v))
        {
            int32_t vtype = vm->classify_variant(odr->hdr, v);
            
            if (vtype==VT_SNP)
            {
                if (bcf_get_n_allele(v)==2)
                {
                    ++no_snp2;
                }
                else if (bcf_get_n_allele(v)==3)
                {
                    ++no_snp3;
                }
                else if (bcf_get_n_allele(v)==4)
                {
                    ++no_snp4;
                }    
            }  
              
        }
        
        odr->close();
    };

    void print_options()
    {
        std::clog << "peek v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [r] reference FASTA file  " << ref_fasta_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\nstats: No. of SNPs          : " << no_snp2+no_snp3+no_snp4 << "\n";
        std::clog << "           biallelic          : " << no_snp2 << "\n";
        std::clog << "           multiallelic       : " << no_snp3+no_snp4 << "\n";
        std::clog << "             3 alleles        : " << no_snp3 << "\n";
        std::clog << "             4 alleles        : " << no_snp4 << "\n";
        std::clog << "\n";
//        std::clog << "         No. of MNPs          : " << no_mnp2+no_mnp3+no_mnp4+no_mnp5 << "\n";
//        std::clog << "            biallelic         : " << no_mnp2 << "\n";
//        std::clog << "            multiallelic      : " << no_mnp3+no_mnp4+no_mnp5 << "\n";
//        std::clog << "             length 2         : " << no_mnp2 << "\n";
//        std::clog << "             length 3         : " << no_mnp3 << "\n";
//        std::clog << "             length 4         : " << no_mnp4 << "\n";
//        std::clog << "             length >=5       : " << no_mnp5 << "\n";
//        std::clog << "\n";
//        std::clog << "         No. of Indels        : " << no_indels << "\n";
//        std::clog << "            biallelic         :   " << no_indels_bi << "\n";
//        std::clog << "             insertions       :   " << no_insertions << "\n";
//        std::clog << "             deletions        :   " << no_deletions << "\n";
//        std::clog << "            multiallelic      : " << no_la << "\n";
//        std::clog << "                              : " << no_la << "\n";
//        std::clog << "          No. SNP Indels      : " << no_la << "\n";
//        std::clog << "                               : " << no_la << "\n";
//        std::clog << "          No. Complex Substitutions                     : " << no_la << "\n";
//        std::clog << "             biallelic                      : " << no_la << "\n";
//        std::clog << "             multiallelic                      : " << no_la << "\n";
//        std::clog << "                               : " << no_la << "\n";
//        std::clog << "          No. of SVs                      : " << no_la << "\n";
//        std::clog << "             precise                      : " << no_la << "\n";
//        std::clog << "               biallelic                      : " << no_la << "\n";
//        std::clog << "              multiallelic                    : " << no_la << "\n";
//        std::clog << "             unprecise                  : " << no_la << "\n";
//        std::clog << "              biallelic                      : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//		std::clog << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//        std::clog << "              multiallelic                     : " << no_la << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
//		std::clog << "\n";
    


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
