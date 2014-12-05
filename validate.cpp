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

#include "normalize.h"

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
    int32_t window_size;
    bool print;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    bcf1_t *v;



    /////////
    //stats//
    /////////
    uint32_t no_variants;

    uint32_t no_lt;    //# left trimmed

    /////////
    //tools//
    /////////
    faidx_t *fai;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "validates variants in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<int32_t> arg_window_size("w", "w", "window size for local sorting of variants [10000]", false, 10000, "integer", cmd);
            TCLAP::SwitchArg arg_quiet("q", "q", "do not print options and summary []", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print = !arg_quiet.getValue();
            window_size = arg_window_size.getValue();
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



        odw = new BCFOrderedWriter(output_vcf_file, window_size, false);
        odw->link_hdr(odr->hdr);
        odw->write_hdr();


        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        if (ref_fasta_file!="")
        {
            fai = fai_load(ref_fasta_file.c_str());
            if (fai==NULL)
            {
                fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
                exit(1);
            }
        }
    }

    void validate()
    {
        uint32_t left_extended = 0;
        uint32_t left_trimmed = 0;
        uint32_t right_trimmed = 0;

        
        v = odw->get_bcf1_from_pool();
        Variant variant;

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            const char* chrom = odr->get_seqname(v);
            int32_t pos1 = bcf_get_pos1(v);

            int32_t ref_len;
            char* ref = faidx_fetch_uc_seq(fai, chrom, pos1-5, pos1+10, &ref_len);

            bcf_print(odr->hdr, v);



            std::cerr << chrom << ":" << pos1 << ":" << ref << "\n";            
//            
//            
//            std::vector<std::string> alleles;
//            for (size_t i=0; i<bcf_get_n_allele(v); ++i)
//            {
//                
//                
//                ru = faidx_fetch_uc_seq(fai, chrom, pos1, pos1+motif_len-1, &ref_len);
//                
//                char *s = bcf_get_alt(v, i);
//                while (*s)
//                {
//                    *s = toupper(*s);
//                    ++s;
//                }
//                alleles.push_back(std::string(bcf_get_alt(v, i)));
//            }
//            left_extended = left_trimmed = right_trimmed = 0;
//
//     
//
//               
//
//            ++no_variants;
//
//            odw->write(v);
//            v = odw->get_bcf1_from_pool();
        }

        odw->close();
        odr->close();
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "validate v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [w] sorting window size   " << window_size << "\n";
        std::clog << "         [r] reference FASTA file  " << ref_fasta_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;



        std::clog << "\n";
        std::clog << "stats: biallelic\n";
        std::clog << "          no. left trimmed                      : " << "\n";

        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool validate(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.validate();
    igor.print_stats();

    return igor.print;
};
