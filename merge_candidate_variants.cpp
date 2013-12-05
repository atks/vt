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

#include "merge_candidate_variants.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string build;
    std::vector<std::string> input_vcf_files;
    std::string input_vcf_file_list;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    /////////
    //stats//
    /////////
    uint32_t no_candidate_snps;
    uint32_t no_candidate_indels;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc =
"Merge candidate variants across samples.\n\
Each VCF file is required to have the FORMAT flags E and N and should have exactly one sample.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", true, "", "str", cmd);

            cmd.parse(argc, argv);

            input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());

            ///////////////////////
            //parse input VCF files
            ///////////////////////
            htsFile *file = hts_open(input_vcf_file_list.c_str(), "r");
            if (file==NULL)
            {
                std::cerr << "cannot open " << input_vcf_file_list.c_str() << "\n";
                exit(1);
            }
            kstring_t *s = &file->line;
            while (hts_getline(file, KS_SEP_LINE, s) >= 0)
            {
                if (s->s[0]!='#')
                {    
                    input_vcf_files.push_back(std::string(s->s));
                }
            }
            hts_close(file);
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
        sr = new BCFSyncedReader(input_vcf_files, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bcf_hdr_append(odw->hdr, "##fileformat=VCFv4.1");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=RL,Number=1,Type=String,Description=\"Length of each read\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=LR,Number=1,Type=String,Description=\"Likelihood Ratio Statistic\">");

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
        no_candidate_indels = 0;
    }

    void merge_candidate_variants()
    {
        //for combining the alleles
        std::vector<bcfptr> current_recs;
        while(sr->read_next_position(current_recs))
        {
            //for each file that contains the next record
            for (uint32_t i=0; i<current_recs.size(); ++i)
            {
                bcf1_t *v = current_recs[i].v;
                std::cerr << "\t[" << current_recs[i].file_index << "]" << v->rid << ":" << (v->pos+1) << ":" << v->d.allele[0] << ":" << v->d.allele[1];
            }
            std::cerr << "\n";
        }

        odw->close();
    };

    void print_options()
    {
        std::clog << "merge_candidate_variants v" << version << "\n\n";
        std::clog << "options:     input VCF file        " << input_vcf_files.size() << " files\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: Total Number of Candidate SNPs     " << no_candidate_snps << "\n";
        std::cerr << "       Total Number of Candidate Indels   " << no_candidate_indels << "\n\n";
    };

    ~Igor()
    {
    };

    private:
};

}

void merge_candidate_variants(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.merge_candidate_variants();
    igor.print_stats();
}

