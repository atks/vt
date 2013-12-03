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

#include "view.h"

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
    std::vector<std::string> samples;
    std::string filter;
    std::string variant;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

	int* imap;
    kstring_t s;	
	
    /////////
    //stats//
    /////////
    uint32_t no_variants;
	uint32_t no_samples; 

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
    	try
    	{
    		std::string desc = "views a VCF file";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
    		TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
    		TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
    		TCLAP::ValueArg<std::string> arg_sample_list("s", "s", "file containing list of sample []", false, "", "file", cmd);
    		TCLAP::ValueArg<std::string> arg_filter("f", "f", "filter expression []", false, "", "exp", cmd);
    		TCLAP::ValueArg<std::string> arg_variant("v", "v", "variant type []", false, "", "exp", cmd);
    		TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
    		TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

    		cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
    		output_vcf_file = arg_output_vcf_file.getValue();
    		parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
    		read_sample_list(samples, arg_sample_list.getValue());
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
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->set_hdr(odr->hdr);
        bcf_hdr_fmt_text(odr->hdr);
        
        //std::vector<int32_t> imap(samples.size(), -1);
        //bcf_hdr_t * h = bcf_hdr_subset_samples(odr->hdr, samples, imap);
        odw->set_hdr(odr->hdr);
        odw->write_hdr();

        s.s = 0;
        s.l = s.m = 0;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
    	no_samples = 0;
	}

    void view()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);
            bcf_get_pos1(v);
            //bcf_set_variant_types(v);
            
            //filter
            
            //subset if necessary
			//int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap)
            
            
            odw->write(v);
            v = odw->get_bcf1_from_pool();
        }

        odw->flush();
    };

    void print_options()
    {
        std::clog << "view v" << version << "\n\n";

		std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
		std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [i] intervals             " << intervals.size() << " intervals" << "\n";
    	std::clog << "\n"; 
	}

    void print_stats()
    {
        std::clog << "stats: no. variants    : " << no_variants << "\n";
        std::clog << "       no. samples     : " << no_samples << "\n";
        std::clog << "\n"; 
    };
 
 	~Igor() {};

    private:
};

}

void view(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.view();
    igor.print_stats();
};
