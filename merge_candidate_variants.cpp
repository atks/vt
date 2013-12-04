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
        
//        boost::logging::core::get()->set_filter
//        (
//            boost::logging::trivial::severity >= logging::trivial::info
//        );
        
        //////////////////////////
        //options initialization//
        //////////////////////////
    	try
    	{
    		std::string desc =
"Merge candidate variants across samples.\n\
Each VCF file is required to have the FORMAT flags E and N and should have exactly one sample.\n\
e.g. vt merge_candidate_variants -o - NA19130.vcf.gz HG00096.vcf.gz\n\n";

       		version = "0.5";
    		TCLAP::CmdLine cmd(desc, ' ', version);
    		TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "File containing list of input VCF files", false, "", "str", cmd);
    		TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "Output VCF file [-]", false, "-", "", cmd);
    		TCLAP::ValueArg<std::string> arg_intervals("i", "i", "Interval (e.g. 20:1000-2000) [all]", false, "all", "str", cmd);
    		TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "File containing list of VCF files", false, "", "str", cmd);
    		TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("input-vcf-files", "Input VCF Files", true, "str", cmd);
		
    		cmd.parse(argc, argv);

    		input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
    		parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            
            ///////////////////////
    		//parse input VCF files
    		///////////////////////
    		if (input_vcf_files.size()==0)
		    {       
		        //reads in file list
                htsFile *file = hts_open(input_vcf_file_list.c_str(), "r");
                if (file==NULL)
                {
                    std::cerr << "cannot open file\n";
                    exit(1);
                }    
                kstring_t *s = &file->line;
                while (hts_getline(file, KS_SEP_LINE, s) >= 0) 
                {
                    input_vcf_files.push_back(std::string(s->s));    
                }
                hts_close(file);  
	        }
		    else
	        {
	            input_vcf_files = arg_input_vcf_files.getValue();
	        }
	        
	        if (input_vcf_files.size()==1)
            {
               // BOOST_LOG_TRIVIAL(warning) << "Only 1 input VCF file specified";
	        }
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
        odw->hdr_append_metainfo("##fileformat=VCFv4.1");
        odw->hdr_append_metainfo("##FORMAT=<ID=RL,Number=1,Type=String,Description=\"Length of each read\">");
        odw->hdr_append_metainfo("##INFO=<ID=LR,Number=1,Type=String,Description=\"Likelihood Ratio Statistic\">");
                
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
	    no_candidate_indels = 0;
	}

    void merge_candidate_variants()
    {
        //list of variants with the same start position
        std::vector<bcf1_t*> bs;
        
        //map to contain the variant type
        std::map<std::string, std::pair<int32_t, int32_t> > bhs;
        std::vector<bcf_hdr_t *> ivcf_hdrs; 
        
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
    
        std::clog << "options:     input VCF file        " << input_vcf_files.size() << "VCF files\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";    
    	    
    	    
	}

    void print_stats()
    {
        std::cerr << "\nStats: Total Number of Candidate SNPs     " << no_candidate_snps << "\n";
        std::cerr << "         Total Number of Candidate Indels   " << no_candidate_indels << "\n\n";
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
   
