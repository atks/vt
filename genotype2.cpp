/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#include "genotype2.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string sample_id;
    std::string input_vcf_file;
    std::string input_sam_file;
    std::string output_vcf_file;
    std::string ref_fasta_file;
    std::string mode;
    bool debug;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    bcf1_t *v;

    BCFOrderedWriter* odw;

    samFile *isam;
    bam_hdr_t *isam_hdr;
    hts_idx_t *isam_idx;
    bam1_t *srec;

    std::vector<GenomeInterval> intervals;
    
    /////////
    //stats//
    /////////
    uint32_t no_snps_genotyped;
	uint32_t no_indels_genotyped;
    uint32_t noDelRefToAlt;
    uint32_t noDelAltToRef;
    uint32_t noInsRefToAlt;
    uint32_t noInsAltToRef;
    uint32_t readExtendedNo;

    /////////
    //tools//
    /////////
    LogTool lt;
    LHMM1 lhmm_ref, lhmm_alt;
    VariantManip *vm;
    
    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
    	try
    	{
    		std::string desc = "Genotypes SNPs, Indels, VNTRs for each sample.\n";
    		    
       		version = "0.5";
    		TCLAP::CmdLine cmd(desc, ' ', version);
    		VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_input_sam_file("b", "b", "input SAM/BAM/CRAM file", true, "", "string", cmd);
          	TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file", false, "-", "string", cmd);
    		TCLAP::ValueArg<std::string> arg_sample_id("s", "s", "sample ID", true, "", "string", cmd);
    		TCLAP::ValueArg<std::string> arg_mode("m", "m", "mode [x]\n"
                 "              d : iterate by read for dense genotyping.\n"
                 "                 (e.g. 50m variants close to one annother).\n"
                 "              s : iterate by sites for sparse genotyping.\n"
                 "                 (e.g. 100 variants scattered over the genome).\n",
                 false, "d", "str", cmd);
    		TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference FASTA file", false, "/net/fantasia/home/atks/ref/genome/human.g1k.v37.fa", "string", cmd);
    		TCLAP::SwitchArg arg_debug("d", "d", "debug alignments", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

    		cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            input_sam_file = arg_input_sam_file.getValue();
    		output_vcf_file = arg_output_vcf_file.getValue();
    		sample_id = arg_sample_id.getValue();
    		parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_fasta_file = arg_ref_fasta_file.getValue();
    		debug = arg_debug.getValue();
    	}
    	catch (TCLAP::ArgException &e)
    	{
    		std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
    		abort();
    	}

        //////////////////////
        //i/o initialization//
        //////////////////////
        //input vcf
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        
        //input sam


        //output vcf
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->set_hdr(odr->hdr);        
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "SAMPLES");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "NSAMPLES");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "E");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "N");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "ESUM");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "NSUM");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "AF");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "LR");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "REFPROBE");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "ALTPROBE");
        bcf_hdr_remove(odw->hdr, BCF_HL_INFO, "PLEN");
        bcf_hdr_add_sample(odw->hdr, strdup(sample_id.c_str()));
        bcf_hdr_add_sample(odw->hdr, NULL);
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allele Depth\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
        
        odw->write_hdr();
        
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_snps_genotyped = 0;
        readExtendedNo = 0;

        vm = new VariantManip(ref_fasta_file);

    };

 	void print_options()
    {
        std::clog << "genotype2 v" << version << "\n\n";

	    std::clog << "Options: Input VCF File   " << input_vcf_file << "\n";
	    std::clog << "         Input BAM File   " << input_sam_file << "\n";
	    std::clog << "         Output VCF File  " << output_vcf_file << "\n";
	    std::clog << "         Sample ID        " << sample_id << "\n\n";
    }

    void print_stats()
    {
	    std::clog << "Stats: SNPs genotyped     " << no_snps_genotyped << "\n";
	    std::clog << "       Indels genotyped   " << no_indels_genotyped << "\n\n";
        std::clog << "       VNTRs genotyped   " << no_indels_genotyped << "\n\n";
    }

 	~Igor()
    {

    };

    void genotype2()
    {
        
    }

    private:   
};

}

void genotype2(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.genotype2();
    igor.print_stats();
}

