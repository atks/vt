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

#include "annotate_str.h"

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
    std::string ref_fasta_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    int32_t no_variants_annotated;

    ////////////////
    //common tools//
    ////////////////
    faidx_t *fai;
    VariantManip *vm;
    RFHMM* rfhmm;
    LFHMM* lfhmm;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "annotates variants in a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_fasta_file = arg_ref_fasta_file.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
            abort();
        }
    };

    ~Igor() {};

    void initialize()
    {
        //******************
        //i/o initialization
        //******************
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RU,Number=1,Type=String,Description=\"Repeat unit in a STR or Homopolymer\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RL,Number=1,Type=Integer,Description=\"Repeat Length\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_LFLANK,Number=1,Type=String,Description=\"Right Flank\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RFLANK,Number=1,Type=String,Description=\"Left Flank\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_LFLANKPOS,Number=2,Type=Integer,Description=\"Location of left flank\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RFLANKPOS,Number=2,Type=Integer,Description=\"Location of right flank\">");

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip(ref_fasta_file);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants_annotated = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        fai = fai_load(ref_fasta_file.c_str());
        rfhmm = new RFHMM();
        lfhmm = new LFHMM();

    }

    void print_options()
    {
        std::clog << "annotate_indels v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)     " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_ref_op("         [r] ref FASTA file        ", ref_fasta_file);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants annotated     " << no_variants_annotated << "\n";
        std::clog << "\n";
    }

    void annotate_str()
    {
        odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        std::vector<Interval*> overlaps;
        Variant variant;
        kstring_t s = {0,0,0};

        int32_t lflank_len, ref_genome_len;
        char* lflank, *ru, *rflank, *ref_genome;
        int32_t ru_len;
        std::string qual;
        for (int32_t i=0; i<2048; ++i) qual += 'K';

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            const char* chrom = bcf_get_chrom(odr->hdr,v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);

            //1. reduce variant
            //   ACACACAC=>AC
            //   check concordance
            //   Reference length
            //   accuracy

            //2. run right flank

            if (bcf_get_n_allele(v)!=2 && !(vtype&VT_INDEL))
            {
                continue;
            }

            
            std::cerr << "\n";
            std::cerr << "=========================================================================================================================\n";
            std::cerr << "\n";


            if (!ru_len)
            {
                char* ref = bcf_get_alt(v, 0);
                char* alt = bcf_get_alt(v, 1);

                if (strlen(ref)>strlen(alt))
                {            
                    lflank = faidx_fetch_uc_seq(fai, chrom, start1-10, start1-1, &lflank_len);
                    //bcf_get_info_string(odr->hdr, v, "RU", &ru, &ru_len);
                    ref_genome = faidx_fetch_uc_seq(fai, chrom, start1-10, start1+100, &ref_genome_len);

                    ru = ref;
                    ++ru;
                }
                else // deletion
                {
                    lflank = faidx_fetch_uc_seq(fai, chrom, start1-10, start1-1, &lflank_len);
                    ru = alt;
                    kstring_t str = {(size_t)lflank_len, (size_t)lflank_len, lflank};
                    ++ru;
                    kputs(ru, &str);
                    lflank_len = str.m;
                    lflank = str.s;
                
                    ref_genome = faidx_fetch_uc_seq(fai, chrom, start1, start1+100, &ref_genome_len);
                    str.l=0; str.s=0; str.m=0;
                    kputs(lflank, &str);
                    kputs(ru, &str);
                    kputs(ref_genome, &str);
                    
                    free(ref_genome);
                    ref_genome = str.s;
                    ref_genome_len = str.l;
                }
            }

            bcf_print(odr->hdr, v);

            std::cerr << "lflank    : " << lflank << "\n";
            std::cerr << "RU        : " << ru << "\n";
            std::cerr << "ref_genome: " << ref_genome << "\n";
            std::cerr << "\n";


            lfhmm->set_model(lflank, ru);
            lfhmm->align(ref_genome, qual.c_str());
            lfhmm->print_alignment();

            bcf_print(odr->hdr, v);

            //check if there are at least 10bp to work with
//            rfhmm->initialize(run, rflank);
//            rfhmm->align(ref_genome, qual.c_str());
//            rfhmm->print_alignment();

            //3. run left flank
            //
            //
            //4. try several modes

            if (lflank_len) free(lflank);
            if (ref_genome_len) free(ref_genome);

            ++no_variants_annotated;
            odw->write(v);
            v = odw->get_bcf1_from_pool();
        }

        odw->close();
        odr->close();
    };

    private:
};
}

void annotate_str(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.annotate_str();
    igor.print_stats();
};
