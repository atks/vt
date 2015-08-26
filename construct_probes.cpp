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

#include "construct_probes.h"

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
    uint32_t min_flank_length;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    kstring_t s;

    /////////
    //stats//
    /////////
    uint32_t no_probes_generated;    //# left trimmed
    uint32_t no_variants; //# left trimmed and left aligned

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
            std::string desc = "Construct probes for genotyping a variant.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_min_flank_length("f", "f", "minimum flank length [20]", false, 20, "int", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            min_flank_length = arg_min_flank_length.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
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

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=REFPROBE,Number=1,Type=String,Description=\"Probe for Determining Reference Allele\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=ALTPROBE,Number=A,Type=String,Description=\"Probe for Determining Alternate Allele(s)\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=PLEN,Number=1,Type=Integer,Description=End location of 5' flank in the probe\"\">");
        odw->write_hdr();

        s = {0,0,0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_probes_generated = 0;
        no_variants = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        var_manip = new VariantManip(ref_fasta_file);
    }

    void construct_probes()
    {
        v = odw->get_bcf1_from_pool();
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            //ignore alleles with N
            if (strchr(bcf_get_alt(v, 0), 'N') || strchr(bcf_get_alt(v, 1), 'N'))
            {
                continue;
            }

            //bcf_print_liten(odr->hdr,v);

            std::vector<std::string> probes;
            std::vector<std::string> alleles;

            int32_t preambleLength = 0;
            
            //populate alleles
            for (int32_t i=0; i<bcf_get_n_allele(v); ++i)
            {
                alleles.push_back(std::string(bcf_get_alt(v,i)));
            }
            
            var_manip->generate_probes(bcf_get_chrom(odr->hdr, v), bcf_get_pos1(v), 1, alleles, probes, min_flank_length, preambleLength);

            //remove ill defined probes
            bool skip = false;
            for (uint32_t i=1; i<probes.size()-1; ++i)
            {
                if(strchr(probes[i].c_str(), 'N'))
                {
                    skip = true;
                }
            }

            ++no_variants;

            if (skip)
                continue;

            bcf_update_info_string(odw->hdr, v, "REFPROBE", probes[0].c_str());
            bcf_update_info_string(odw->hdr, v, "ALTPROBE", probes[1].c_str());
            bcf_update_info_int32(odw->hdr, v, "PLEN", &preambleLength, 1);

            odw->write(v);
            v = odw->get_bcf1_from_pool();

            ++no_probes_generated;
        }

        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "construct_probes v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [r] reference FASTA file  " << ref_fasta_file << "\n";
        std::clog << "         [f] minimum flank length  " << min_flank_length << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. of probes generated      : " << no_probes_generated << "\n";
        std::clog << "       no. variants                 : " << no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};
}

void construct_probes(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.construct_probes();
    igor.print_stats();
};
