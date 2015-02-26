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

#include "validate.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string ref_fasta_file;
    bool print;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    bcf1_t *v;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    uint32_t no_unordered;
    uint32_t no_unordered_chrom;
    uint32_t no_inconsistent_ref;

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
            std::string desc = "Checks the following properties of a VCF file\n"
                 "              1. order\n"
                 "              2. reference sequence consistency";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::SwitchArg arg_quiet("q", "q", "do not print invalid records [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print = !arg_quiet.getValue();
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

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_unordered = 0;
        no_unordered_chrom = 0;
        no_inconsistent_ref = 0;
        v = bcf_init1();

        ////////////////////////
        //tools initialization//
        ////////////////////////
        fai = NULL;
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
        const char* last_chrom = NULL;
        int32_t last_rid = -1;
        int32_t last_pos1 = -1;

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            int32_t rid = bcf_get_rid(v);
            const char* chrom = odr->get_seqname(v);
            int32_t pos1 = bcf_get_pos1(v);

            if (rid==last_rid)
            {
                if (last_pos1>pos1)
                {
                    if (print) fprintf(stderr, "[%s:%d %s] UNORDERED: %s: %d after %d\n", __FILE__, __LINE__, __FUNCTION__, chrom, pos1, last_pos1);
                    ++no_unordered;
                }
            }
            else if (last_rid>rid)
            {
                if (print) fprintf(stderr, "[%s:%d %s] UNORDERED CHROM: %s after %s\n", __FILE__, __LINE__, __FUNCTION__, chrom, last_chrom);
                ++no_unordered_chrom;
            }

            last_chrom = chrom;
            last_rid = rid;
            last_pos1 = pos1;

            if (fai)
            {
                const char* chrom = odr->get_seqname(v);

                int32_t ref_len = 0;
                char** alleles = bcf_get_allele(v);
                int32_t len = strlen(alleles[0]);

                char* ref = faidx_fetch_uc_seq(fai, chrom, pos1-1, pos1+len-2, &ref_len);

                if (strcasecmp(ref, alleles[0]))
                {
                    kstring_t var = {0,0,0};
                    bcf_variant2string(odr->hdr, v, &var);
                    if (print) fprintf(stderr, "[%s:%d %s] INCONSISTENT REF: %s:%d %s!=%s(truth) for variant %s\n", __FILE__, __LINE__, __FUNCTION__, chrom, pos1, alleles[0], ref, var.s);
                    if (var.m) free(var.s);
                    ++no_inconsistent_ref;
                }

                if (ref_len)
                {
                    free(ref);
                }
            }

            ++no_variants;
        }

        odr->close();
    };

    void print_options()
    {
        std::clog << "validate v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        print_ref_op("         [r] reference FASTA file  ", ref_fasta_file);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats:    no. unordered                     : " << no_unordered << "\n";
        std::clog << "          no. unordered chrom               : " << no_unordered_chrom << "\n";
        std::clog << "\n";
        if (fai)
        {
            std::clog << "          no. inconsistent REF              : " << no_inconsistent_ref << "\n";
            std::clog << "\n";
        }
        else
        {
            std::clog << "reference consistency not checked.\n\n"; 
        }
        std::clog << "          no. variants                      : " << no_variants << "\n";
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