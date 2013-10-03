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
    std::vector<std::string> intervals;
    std::string build;
    std::string ref_fasta_file;

    ///////
    //i/o//
    ///////
    OrderedReader *odr;
    OrderedWriter *odw;
    bcf1_t *v;

    kstring_t s;
    kstring_t old_allele;

    /////////
    //stats//
    /////////
    uint32_t no_lt;    //# left trimmed
	uint32_t no_lt_la; //# left trimmed and left aligned
	uint32_t no_lt_rt; //# left trimmed and right trimmed
	uint32_t no_la;    //# left aligned
	uint32_t no_rt;    //# right trimmed

	uint32_t no_multi_lt;      //# left trimmed
	uint32_t no_multi_lt_la;   //# left trimmed and left aligned
	uint32_t no_multi_lt_rt;   //# left trimmed and right trimmed
	uint32_t no_multi_la;      //# left aligned
	uint32_t no_multi_rt;      //# right trimmed

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
        odr = new OrderedReader(input_vcf_file, intervals);
        odw = new OrderedWriter(output_vcf_file);
        bcf_hdr_append(odr->hdr, "##INFO=<ID=OLD_VARIANT,Number=1,Type=String,Description=\"Original chr:pos:ref:alt encoding\">\n");
        bcf_add_hs37d5_contig_headers(odr->hdr);
        bcf_hdr_fmt_text(odr->hdr);
        odw->set_hdr(odr->hdr);
        odw->write_hdr();

        s.s = 0;
        s.l = s.m = 0;

        old_allele.s = 0;
        old_allele.l = old_allele.m = 0;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_lt = 0;
    	no_lt_la = 0;
    	no_lt_rt = 0;
    	no_la = 0;
    	no_rt = 0;

    	no_multi_lt = 0;
    	no_multi_lt_la = 0;
    	no_multi_lt_rt = 0;
    	no_multi_la = 0;
    	no_multi_rt = 0;

    	////////////////////////
        //tools initialization//
        ////////////////////////
    	var_manip = new VariantManip(ref_fasta_file);
	}

    void normalize()
    {
        v = odw->get_bcf1_from_pool();

        while (odr->read1(v))
        {
            bcf_unpack(v, BCF_UN_INFO);
            bcf_get_pos1(v);
            bcf_set_variant_types(v);

            if (bcf_get_var_type(v) == VCF_INDEL || bcf_get_var_type(v) == VCF_OTHER ||bcf_get_var_type(v) == VCF_MNP )
            {
                const char* chrom = odr->get_seqname(v);
                uint32_t pos1 = bcf_get_pos1(v);
                
                //transfer alleles into a vector (not really a good choice of a data structure ...)
                std::vector<std::string> alleles;
                for (uint32_t i=0; i<bcf_get_n_allele(v); ++i)
                {
                    alleles.push_back(std::string(bcf_get_alt(v, i)));
                }

                uint32_t left_aligned = 0;
    	    	uint32_t left_trimmed = 0;
                uint32_t right_trimmed = 0;

                var_manip->left_align(alleles, pos1, odr->get_seqname(v), left_aligned, right_trimmed);
                var_manip->left_trim(alleles, pos1, left_trimmed);

                bool changed = left_trimmed || left_aligned || right_trimmed;

        	    if (changed)
        	    {
        	        old_allele.l=0; //this is a hack. why the hell was bcf implemented to be not modifiable friendly??
        	        bcf_format_variant(odr->hdr, v, &old_allele);

        	        bcf_set_pos1(v, pos1);
        	        bcf_set_allele(v, alleles);
        	        bcf_add_info(odr->hdr, v, BCF_BT_CHAR, const_cast<char *>("OLD_VARIANT"), (uint8_t*)old_allele.s, old_allele.l);

        	        vcf_format1(odw->hdr, v, &s);
        	        vcf_parse1(&s, odw->hdr, v);

        	        if (bcf_get_n_allele(v)==2)
        	        {
        	            if (left_trimmed)
        	            {
            	            if (left_aligned)
            	            {
            	                ++no_lt_la;
            	            }
            	            else if (right_trimmed)
            	            {
            	                ++no_lt_rt;
            	            }
            	            else
        	                {
        	                    ++no_lt;
        	                }
        	            }
        	            else
        	            {
        	                if (left_aligned)
            	            {
            	                ++no_la;
            	            }
            	            else if (right_trimmed)
            	            {
            	                ++no_rt;
            	            }
        	            }
        	        }
        	        else
    	            {
    	                if (left_trimmed)
        	            {
            	            if (left_aligned)
            	            {
            	                ++no_multi_lt_la;
            	            }
            	            else if (right_trimmed)
            	            {
            	                ++no_multi_lt_rt;
            	            }
            	            else
        	                {
        	                    ++no_multi_lt;
        	                }
        	            }
        	            else
        	            {
        	                if (left_aligned)
            	            {
            	                ++no_multi_la;
            	            }
            	            else if (right_trimmed)
            	            {
            	                ++no_multi_rt;
            	            }
        	            }
    	            }
    	        }
            }

            odw->write1(v);
            v = odw->get_bcf1_from_pool();
        }

        odw->flush();

    };

    void print_options()
    {
        std::clog << "normalize v" << version << "\n\n";

		std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
		std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [r] reference FASTA file  " << ref_fasta_file << "\n";
        std::clog << "         [i] intervals             " << (intervals.size()==0? "NONE" : boost::algorithm::join(intervals, ",")) << "\n";
    	std::clog << "\n"; 
	}

    void print_stats()
    {
        std::clog << "stats: biallelic\n";
        std::clog << "          no. left trimmed                      : " << no_lt << "\n";
        std::clog << "          no. left trimmed and left aligned     : " << no_lt_la << "\n";
        std::clog << "          no. left trimmed and right trimmed    : " << no_lt_rt << "\n";
        std::clog << "          no. left aligned                      : " << no_la << "\n";
        std::clog << "          no. right trimmed                     : " << no_rt << "\n";
        std::clog << "\n";
        std::clog << "       multiallelic\n";
        std::clog << "          no. left trimmed                      : " << no_multi_lt << "\n";
        std::clog << "          no. left trimmed and left aligned     : " << no_multi_lt_la << "\n";
        std::clog << "          no. left trimmed and right trimmed    : " << no_multi_lt_rt << "\n";
        std::clog << "          no. left aligned                      : " << no_multi_la << "\n";
        std::clog << "          no. right trimmed                     : " << no_multi_rt << "\n";
        std::clog << "\n"; 
    };

 	~Igor()
    {

    };

    private:
};

}

void normalize(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.normalize();
    igor.print_stats();
};
