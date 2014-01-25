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

#include "decompose.h"

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
    BCFOrderedWriter *odw;
    bcf1_t *v;

    kstring_t s;
    kstring_t new_alleles;
    kstring_t old_alleles;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    uint32_t no_biallelic;
    uint32_t no_multiallelic;
    uint32_t no_additional_biallelic;

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
            std::string desc = "decomposes multialleic variants into biallelic in a VCF file, only sites are output.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
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

        odw = new BCFOrderedWriter(output_vcf_file, 100000);
        odw->link_hdr(bcf_hdr_subset(odr->hdr, 0, 0, 0));
        bcf_hdr_append(odw->hdr, "##INFO=<ID=OLD_MULTIALLELIC,Number=1,Type=String,Description=\"Original chr:pos:ref:alt encoding\">\n");
        odw->write_hdr();

        s = {0,0,0};
        old_alleles = {0,0,0};
        new_alleles = {0,0,0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_biallelic = 0;
        no_multiallelic = 0;
        no_additional_biallelic = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void decompose()
    {
        v = odw->get_bcf1_from_pool();
        Variant variant;

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            int32_t n_allele = bcf_get_n_allele(v);
            if (n_allele > 2)
            {
                ++no_multiallelic;
                no_additional_biallelic += n_allele-1;

                old_alleles.l = 0;
                bcf_variant2string(odw->hdr, v, &old_alleles);

                int32_t rid = bcf_get_rid(v);
                int32_t pos1 = bcf_get_pos1(v);
                char** allele = bcf_get_allele(v);
                for (uint32_t i=1; i<n_allele; ++i)
                {
                    bcf_set_rid(v, rid);
                    bcf_set_pos1(v, pos1);
                    new_alleles.l=0;
                    kputs(allele[0], &new_alleles);
                    kputc(',', &new_alleles);
                    kputs(allele[i], &new_alleles);
                    bcf_update_alleles_str(odw->hdr, v, new_alleles.s);
                    bcf_update_info_string(odw->hdr, v, "OLD_MULTIALLELIC", old_alleles.s);
                    bcf_subset(odw->hdr, v, 0, 0);
                    odw->write(v);
                    v = odw->get_bcf1_from_pool();
                }
            }
            else
            {
                ++no_biallelic;

                bcf_subset(odw->hdr, v, 0, 0);
                odw->write(v);
                v = odw->get_bcf1_from_pool();
            }

            ++no_variants;
        }

        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "decompose v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. variants                 : " << no_variants << "\n";
        std::clog << "       no. biallelic variants       : " << no_biallelic << "\n";
        std::clog << "       no. multiallelic variants    : " << no_multiallelic << "\n";
        std::clog << "\n";
        std::clog << "       no. additional biallelics    : " << no_additional_biallelic << "\n";
        std::clog << "       after decomposition\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void decompose(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.decompose();
    igor.print_stats();
};
