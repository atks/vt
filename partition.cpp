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

#include "partition.h"

namespace
{
    
class OverlapStats
{
    public:
        
    uint32_t a,ab,b;
    
    OverlapStats()
    {
        a = 0;
        ab = 0;
        b = 0;
    };
};

class Igor : Program
{
    public:

    std::string version; 

    ///////////
    //options//
    ///////////
    std::string filters;  
    std::vector<std::string> input_vcf_files;   
   	
    std::string build;
    std::vector<GenomeInterval> intervals;
    
    std::string interval_list;
   
    std::vector<OverlapStats> stats;
        
    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    bcf1_t *v;
    Filter *filter;
    kstring_t line;
    
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
"Paritition variants.  Check the overlap of variants between 2 data sets.\n";

       		version = "0.5";
    		TCLAP::CmdLine cmd(desc, ' ', version);
			VTOutput my;
			cmd.setOutput(&my);    		
			TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
			TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
    		TCLAP::ValueArg<std::string> arg_filters("f", "filters", "Filter (e.g. AF>0.3) ", false, "", "str", cmd);
    		TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in.vcf>", "input VCF file", true, "file", cmd);
		
    		cmd.parse(argc, argv);
    	
   			filters = arg_filters.getValue();
			parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
			input_vcf_files = arg_input_vcf_files.getValue();
    	    
    	    filter=NULL;		            
    		if (filters!="")
		    {
		        filter = new Filter();
		        filter->parse(filters);
		    }
		    
		    filter = new Filter("AF", 0, 0.005);
    		
	        ///////////////////////
    		//parse input VCF files
    		///////////////////////
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
        line = {0,0,0};
                
        //input vcfs
        sr = new BCFSyncedReader(input_vcf_files, intervals); 
            
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
	    no_candidate_indels = 0;
	}

    void partition()
    {        
        //map to contain the variant type
        std::vector<bcf_hdr_t *> ivcf_hdrs; 
        
        //for combining the alleles
        std::vector<bcfptr> current_recs;
        std::map<std::string, std::map<std::string, int32_t> > variants;    
        std::stringstream ss;
        std::stringstream centers;
        
        std::vector<Interval*> intervals;
            
        while(sr->read_next_position(current_recs))
        {
            variants.clear();
            ss.str("");
            centers.str("");
                
            //for each file that contains the next record
            char *ref = const_cast<char*>("N");
            char *alt = const_cast<char*>("N");
            bool in_ref = false;    
            for (uint32_t i=0; i<current_recs.size(); ++i)
            {
                int32_t d = current_recs[i].file_index;
                bcf1_t *v = current_recs[i].v;
                //bcf_set_variant_types(v);

                if (d==0)
                {
                    bcf_hdr_t *h = sr->hdrs[0];
                    if (!bcf_is_passed(h, v))
                    {    
                        continue;
                    }
                    
                    if (filter!=NULL && !filter->apply(h, v))
                    {   
                        continue;
                    }
                }
                
                if (bcf_get_var_type(v)==VCF_SNP || bcf_get_var_type(v)==VCF_MNP || bcf_get_n_allele(v)!=2)
                {
                   continue;
                }    
                
                if (d==0)
                {
                    ref = bcf_get_alt(v, 0);
                    alt = bcf_get_alt(v, 1);
                    
                    int32_t pos1 = bcf_get_pos1(v);
                    int32_t ref_len = strlen(ref);
                    
                    bcf_hdr_t *h = sr->hdrs[0];
                    
                    in_ref = true;
    
                    for (uint32_t j=1; j<2; ++j)
                    {
                        ++stats[j].a;
                    }
                }
                else
                {
                    update_stats(d, v, in_ref, ref, alt);
                     
                    if (d==7)
                    {   
                        //rare
                        bcf_hdr_t *h = sr->hdrs[7];
                        if (filter->apply(h, v))
                        {   
                           update_stats(8, v, in_ref, ref, alt);
                        }
                        else
                        {
                           update_stats(9, v, in_ref, ref, alt); 
                        }
                    }    
                }
            }
        }
    };

    void update_stats(int32_t d, bcf1_t *v, bool in_ref, char* ref, char* alt)
    {
        char* r = bcf_get_alt(v, 0);
        char* a = bcf_get_alt(v, 1);
        
        int32_t ins = 0;
        int32_t del = 0;
        
        if (strlen(r) > strlen(a))
        {
            ++del;
        }
        else
        {
            ++ins;
        }
            
        if (in_ref)
        {
            if (!strcmp(ref,r) && !strcmp(alt,a))
            {
                --stats[d].a;
                ++stats[d].ab;
            }
            else
            {
                ++stats[d].b;
            }
        }
        else
        {
            ++stats[d].b; 
        }
    }

    void print_options()
    {
   		std::clog << "partition v" << version << "\n\n";
    
    	std::clog << "Options:     Input VCF File 1  " << input_vcf_files[0] << "\n";
    	std::clog << "             Input VCF File 2  " << input_vcf_files[1] << "\n";
    	std::clog << "         [f] Filters           " << filters << "\n";
   }

    void print_stats()
    {      
//        for (uint32_t j=1; j<dataset_labels.size(); ++j)
//        {
//            double insdel_a_ab =  (stats[j].a_del+stats[j].ab_del)==0 ? -1 : ((double)(stats[j].a_ins+stats[j].ab_ins)/(double)(stats[j].a_del+stats[j].ab_del)); 
//            double insdel_a =  stats[j].a_del==0 ? -1 : ((double)stats[j].a_ins/(double)stats[j].a_del); 
//            double insdel_ab =  stats[j].ab_del==0 ? -1 : ((double)stats[j].ab_ins/(double)stats[j].ab_del); 
//            double insdel_b =  stats[j].b_del==0 ? -1 : ((double)stats[j].b_ins/(double)stats[j].b_del); 
//            double insdel_b_ab =  (stats[j].b_del+stats[j].ab_del)==0 ? -1 : ((double)(stats[j].b_ins+stats[j].ab_ins)/(double)(stats[j].b_del+stats[j].ab_del)); 
//         
//            uint32_t totala = stats[j].a+stats[j].ab;
//            double ab_of_a = totala==0 ? -1 : ((double)stats[j].ab/(double)(totala));
//             
//            uint32_t totalb = stats[j].b+stats[j].ab;
//            double ab_of_b = totalb==0 ? -1 : ((double)stats[j].ab/(double)(totalb));            
//        
//            if (j==1)
//            {
//                std::cout << std::setprecision(3);
//                std::cout << "\n";
//                std::cout << dataset_labels[0] << "\n";
//                std::cout << "variants: " << totala << "\n";
//                std::cout << "ins/del ratio: " << insdel_a_ab << "\n"; 
//                std::cout << "FS Proportion: " << (float)fs/(float)(fs+nfs) << " (" << fs << "," << nfs << ")\n\n";
//                std::cout << "FS Proportion (Rare): " << (float)rare_fs/(float)(rare_fs+rare_nfs) << " (" << rare_fs << "," << rare_nfs << ")\n\n";
//                std::cout << "FS Proportion (Common): " << (float)common_fs/(float)(common_fs+common_nfs) << " (" << common_fs << "," << common_nfs << ")\n\n";
//            
//            }    
//            
//            std::cout << std::setprecision(3);
//            std::cout << dataset_labels[j] << " (" << totalb << ") " << "[" << insdel_b_ab << "]\n";
//            std::cout << "TP\t" << ab_of_b << "(" << stats[j].ab << "/" << totalb << ") " << "[" << insdel_ab << "," << insdel_b << "] \n";
//            std::cout << "FP\t" << ab_of_a << "(" << stats[j].ab << "/" << totala << ") " << "[" << insdel_ab << "," << insdel_a << "] \n\n";
//        }
    };
	
 	~Igor()
    {
    };
    
    private:
 
};

}

void partition(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.partition();
    igor.print_stats();
}
 