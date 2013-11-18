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

    kstring_t new_alleles;

    /////////
    //stats//
    /////////
    uint32_t no_normalized;
    uint32_t no_lt;    //# left trimmed
    uint32_t no_lt_la; //# left trimmed and left aligned

    /////////
    //tools//
    /////////
    VariantManip *var_manip;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "normalizes variants in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
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

        new_alleles = {0, 0, 0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_lt = 0;
        no_lt_la = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        var_manip = new VariantManip(ref_fasta_file);
    }

    int32_t classify_variant(bcf_hdr_t *h, bcf1_t *v, std::string& motif, uint32_t& tlen)
    {
        return var_manip->classify_variant(bcf_get_chrom(h, v), bcf_get_pos1(v), bcf_get_allele(v), bcf_get_n_allele(v), motif, tlen);
    }

    void peek()
    {
        uint32_t left_aligned = 0;
        uint32_t left_trimmed = 0;
        uint32_t right_trimmed = 0;
        std::vector<std::string> alleles;
    
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            const char* chrom = odr->get_seqname(v);
            uint32_t pos1 = bcf_get_pos1(v);
            
            alleles.clear();
            for (uint32_t i=0; i<bcf_get_n_allele(v); ++i)
            {
                alleles.push_back(std::string(bcf_get_alt(v, i)));
            }                
            left_aligned = left_trimmed = right_trimmed = 0;
            var_manip->left_align(alleles, pos1, chrom, left_aligned, right_trimmed);
            var_manip->left_trim(alleles, pos1, left_trimmed);
            if (left_trimmed || left_aligned || right_trimmed)
            {
                bcf_set_pos1(v, pos1);
                new_alleles.l = 0;
                for (uint32_t i=0; i<alleles.size(); ++i)
                {
                    if (i) kputc(',', &new_alleles);
                    kputs(alleles[i].c_str(), &new_alleles);
                }
                bcf_update_alleles_str(odr->hdr, v, new_alleles.s);
                
                ++no_normalized;
			}
            
            std::string motif;
            uint32_t tlen;
            int32_t vtype = classify_variant(odr->hdr, v, motif, tlen);
			
            if (bcf_get_n_allele(v)==2)
            {
                bcf_set_pos1(v, pos1);
                new_alleles.l=0;
                for (uint32_t i=0; i<alleles.size(); ++i)
                {
                    if (i) kputc(',', &new_alleles);
                    kputs(alleles[i].c_str(), &new_alleles);
                }
                bcf_update_alleles_str(odr->hdr, v, new_alleles.s);
            }
            else if (bcf_get_n_allele(v)>2)
            {
            }
            else
            {
                //ignore non variants
            }
                

        }
    };

    void print_options()
    {
        std::clog << "peek v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [r] reference FASTA file  " << ref_fasta_file << "\n";
        if (intervals.size()!=0)
        {
            std::clog << "         [i] intervals             " << intervals.size() <<  " intervals\n";
        }
        std::clog << "\n";
    }

    void print_stats()
    {
//        std::clog << "\nstats: No. of SNPs          : " << no_snps  << "\n";
//        std::clog << "           biallelic          :   " << no_snps_bi << "\n";
//        std::clog << "           multiallelic       :   " << no_snps_multi << "\n";
//        std::clog << "             3 alleles        :   " << no_snps3 << "\n";
//        std::clog << "             4 alleles        :   " << no_snps4 << "\n";
//        std::clog << "\n";
//        std::clog << "         No. of MNPs          : " << no_mnps << "\n";
//        std::clog << "            biallelic         : " << no_mnps_bi << "\n";
//        std::clog << "            multiallelic      : " << no_mnps_multi << "\n";
//        std::clog << "             length 2         : " << no_mnps2 << "\n";
//        std::clog << "             length 3         : " << no_mnps3 << "\n";
//        std::clog << "             length 4         : " << no_mnps4 << "\n";
//        std::clog << "             length >=5       : " << no_mnps5 << "\n";
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
