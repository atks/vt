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

#include "profile_indels.h"

namespace
{
class BEDRecord: public Interval
{
    public:
    std::string chrom;
    
    BEDRecord(std::string& _chrom, uint32_t _start, uint32_t _end)
    {
        chrom = _chrom;
        start = _start;
        end = _end;
    };
    
    void print()
    {
        std::cerr << "chrom   : " << chrom << "\n";
        std::cerr << "[" << start << "," << end << "]\n";
    };
     
    private:    
};

class GTFRecord: public Interval
{
    public:
    std::string gene;
    std::string feature;
    std::string chrom;
    char strand;
    int32_t frame;
    int32_t exonNo;
    bool fivePrimeConservedEssentialSpliceSite;
    bool threePrimeConservedEssentialSpliceSite; 
    bool containsStartCodon;
    bool containsStopCodon;
    uint32_t level;
    std::string attrib;
    
    GTFRecord(std::string& _chrom, uint32_t _start, uint32_t _end, char _strand,
              std::string& _gene, std::string& _feature, int32_t _frame, int32_t _exonNo, 
              bool _fivePrimeConservedEssentialSpliceSite, bool _threePrimeConservedEssentialSpliceSite, 
              bool _containsStartCodon, bool _containsStopCodon,
              uint32_t _level, std::string& _attrib)
    {
        chrom = _chrom;
        start = _start;
        end = _end;
        strand = _strand;
        gene = _gene;
        feature = _feature;
        frame = _frame;
        exonNo = _exonNo;
        fivePrimeConservedEssentialSpliceSite = _fivePrimeConservedEssentialSpliceSite;
        threePrimeConservedEssentialSpliceSite = _threePrimeConservedEssentialSpliceSite;
        containsStartCodon = _containsStartCodon;
        containsStopCodon = _containsStopCodon;
        level = _level;
        attrib = _attrib;
    };
    
    void print()
    {
        std::cerr << "chrom   : " << chrom << "\n";
        std::cerr << "[" << start << "," << end << "]\n";
        std::cerr << "strand                    : " << strand << "\n";
        std::cerr << "address                   : " << this << "\n";
        std::cerr << "gene                      : " << gene << "\n";
        std::cerr << "feature                   : " << feature << "\n";
        std::cerr << "frame                     : " << frame << "\n";
        std::cerr << "exon number               : " << exonNo << "\n";
        std::cerr << "5' conserved splice site  : " << fivePrimeConservedEssentialSpliceSite << "\n";
        std::cerr << "3' conserved splice site  : " << threePrimeConservedEssentialSpliceSite << "\n";
        std::cerr << "contains start codon      : " << containsStartCodon << "\n";
        std::cerr << "contains stop codon       : " << containsStopCodon << "\n";
        std::cerr << "level                     : " << level << "\n";
        std::cerr << "attrib                    : " << attrib << "\n";
    };
     
    private:    
};    
    
    
class OverlapStats
{
    public:
        
    uint32_t a,ab,b,a_ins,ab_ins,b_ins,a_del,ab_del,b_del;
    
    OverlapStats()
    {
        a = 0;
        ab = 0;
        b = 0;
        
        a_ins=0;
        a_del=0;
        ab_ins=0;
        ab_del=0;
        b_ins=0;
        b_del=0;
    };
};

class Igor : Program
{
    public:

    std::string version; 

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;    
    std::vector<std::string> dataset_labels;
    std::string filters;  
    std::vector<std::string> input_ref_files;    
   	std::vector<std::string> input_vcf_files;   
   	
    std::string build;
    std::vector<GenomeInterval> intervals;
    
    std::string interval_list;
   
    std::vector<OverlapStats> stats;
        
    std::string gencode_gtf_file;
    std::map<std::string, IntervalTree*> GENCODE; 
    std::vector<Interval*> exons;
    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    bcf1_t *v;
    Filter *filter;
    kstring_t line;
    Filter *rare_filter;
    
    /////////
    //stats//
    /////////
    uint32_t no_candidate_snps;
	uint32_t no_candidate_indels;
	uint32_t nfs;
	uint32_t fs;
    uint32_t rare_nfs;
	uint32_t rare_fs;
    uint32_t common_nfs;
	uint32_t common_fs;
    
    Igor(int argc, char ** argv)
    {   
        //////////////////////////
        //options initialization//
        //////////////////////////
    	try
    	{
    		std::string desc =   		    
"Profile Indels.\n\
Each VCF file is required to have the FORMAT flags E and N and should have exactly one sample.\n\
e.g. vt profile_snps_variants -o - NA19130.vcf.gz HG00096.vcf.gz\n\n";

       		version = "0.5";
    		TCLAP::CmdLine cmd(desc, ' ', version);
			VTOutput my;
			cmd.setOutput(&my);    		
			TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
			TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
    		TCLAP::ValueArg<std::string> arg_filters("f", "filters", "Filter (e.g. AF>0.3) ", false, "", "str", cmd);
    		TCLAP::ValueArg<std::string> arg_ref_dataset_list("f", "filters", "Filter (e.g. AF>0.3) ", false, "", "str", cmd);
    		TCLAP::ValueArg<std::string> arg_dataset_labels("l", "data-set-labels", "List of names of Data sets, should be the same number as the number of input files", false, "", "str", cmd);
    		TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
		
    		cmd.parse(argc, argv);
    	
   			filters = arg_filters.getValue();
			parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
			input_vcf_file = arg_input_vcf_file.getValue();
    	    input_vcf_files.push_back(input_vcf_file);
    	    
    	    filter=NULL;		            
    		if (filters!="")
		    {
		        filter = new Filter();
		        filter->parse(filters);
		    }
		    
		    rare_filter = new Filter("AF", 0, 0.005);
    		
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
        line.s=0;
        line.l=line.m=0;
        
        gencode_gtf_file = "/net/fantasia/home/atks/ref/encode/gencode.v15.annotation.gtf.gz";
        
        //input_vcf_files.push_back(input_vcf_file);
        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/1000g_phase1.indels.cplxsubs.sites.vcf.gz");
        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/mills.doublehit.gatk.qc.indels.cplxsubs.sites.vcf.gz");
        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/mills.chip.1000g_phase3.poly.indels.cplxsubs.sites.vcf.gz");
        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/affymetrix.exome.chip.1000g_phase3.poly.indels.cplxsubs.sites.vcf.gz");
        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/cg.mendel.errors.indels.cplxsubs.sites.vcf.gz");  
        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/affymetrix.exome.chip.mono.indels.cplxsubs.sites.vcf.gz");
        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/cg.1000g_phase3.poly.indels.cplxsubs.sites.vcf.gz");        
        
        //input vcfs
        sr = new BCFSyncedReader(input_vcf_files, intervals); 
        
        bcf_hdr_append(sr->hdrs[0], "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
        bcf_hdr_append(sr->hdrs[0], "##INFO=<ID=FR,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
            
        std::stringstream ss;   
        for (uint32_t i=0; i<intervals.size(); ++i)
        {
             ss.str("");  
            ss << "chr" << intervals[i].to_string();
            populate_gencode_tree(ss.str().c_str());
        }       
        fs = 0;
        nfs = 0;
        rare_fs = 0;
        rare_nfs = 0;
        common_fs = 0;
        common_nfs = 0;
        
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
	    no_candidate_indels = 0;
	}

    void profile_indels()
    {        
        //map to contain the variant type
        std::vector<bcf_hdr_t *> ivcf_hdrs; 
        
        //for combining the alleles
        std::vector<bcfptr> current_recs;
        std::map<std::string, std::map<std::string, int32_t> > variants;    
        std::stringstream ss;
        std::stringstream centers;
        
        dataset_labels = {"data", "1000g_phase1", "mills_doublehit", "mills_chip", "affymetrix", "cg_mendel", "affy_mono", "cg", "cg_rare" , "cg_common"};
        stats.resize(dataset_labels.size());
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
//                        std::cerr << "FAIL PASS\n";
//                        vcf_format1(h, v, &line);    
//                        std::cerr << line.s << "\n";    
                        continue;
                    }
                    
                    if (filter!=NULL && !filter->apply(h, v))
                    {   
//                        std::cerr << "FAIL AF filter\n"; 
//                        vcf_format1(h, v, &line);    
//                        std::cerr << line.s << "\n";   
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
                    //std::cerr << bcf_get_chrom(h, v) << ":" << (pos1+1) << ":" << (pos1+ref_len-1) << "\n";
                                         
                    GENCODE[std::string(bcf_get_chrom(h, v))]->search(pos1+1, pos1+ref_len-1, exons);
//                    
//                    if ($start<=$F[$END]-1 && $end>=$F[$START])
//					{
//						my @alts = split(",", $alt);
//						my $refLength = length($ref); 
//						my %diff = ();
//						$diff{($refLength-length($_))%3} = 1 foreach @alts;
//																		
//						if (exists($diff{1})||exists($diff{2}))
//						{
//							$GENCODE{"CDS_FS"} = 1;
//						}
//						else
//						{
//							$GENCODE{"CDS_NFS"} = 1;
//						}
//					}
                    
                    if (exons.size()!=0)
                    {  
                        int32_t alt_len = strlen(alt);
                        bool not_frame_shift = ((abs(ref_len-alt_len)%3)==0);
                        if (not_frame_shift)
                        {
                            nfs += exons.size();
                        }
                        else
                        {
                            fs += exons.size();
                        }
                        
                        if (rare_filter->apply(h, v))
                        {   
                            if (not_frame_shift)
                            {
                                rare_nfs += exons.size();
                            }
                            else
                            {
                                rare_fs += exons.size();
                            }
                        }
                        else
                        {
                            if (not_frame_shift)
                            {
                                common_nfs += exons.size();
                            }
                            else
                            {
                                common_fs += exons.size();
                            }
                        }   
                    }  
                    
                    in_ref = true;
                    int32_t ins = 0;
                    int32_t del = 0;
                    
                    if (strlen(ref) > strlen(alt))
                    {
                        ++del;
                    }
                    else
                    {
                        ++ins;
                    }
                        
                    for (uint32_t j=1; j<dataset_labels.size(); ++j)
                    {
                        ++stats[j].a;
                        stats[j].a_ins += ins;
                        stats[j].a_del += del;
                    }
                }
                else
                {
                    update_stats(d, v, in_ref, ref, alt);
                     
                    if (d==7)
                    {   
                        //rare
                        bcf_hdr_t *h = sr->hdrs[7];
                        if (rare_filter->apply(h, v))
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
                stats[d].a_ins -= ins;
                stats[d].a_del -= del;
                ++stats[d].ab;
                stats[d].ab_ins += ins;
                stats[d].ab_del += del;
            }
            else
            {
                ++stats[d].b;
                stats[d].b_ins += ins;
                stats[d].b_del += del; 
            }
        }
        else
        {
            ++stats[d].b; 
            stats[d].b_ins += ins;
            stats[d].b_del += del;                       
        }
    }

    void print_options()
    {
   		std::clog << "profile_indels v" << version << "\n\n";
    
    	std::clog << "Options: [_] Input VCF File    " << input_vcf_files.size() << "\n";
    	std::clog << "         [i] Filters           " << filters << "\n";
    	std::clog << "         [i] intervals         " << intervals.size() << " intervals\n";
   }

    void print_stats()
    {      
        for (uint32_t j=1; j<dataset_labels.size(); ++j)
        {
            double insdel_a_ab =  (stats[j].a_del+stats[j].ab_del)==0 ? -1 : ((double)(stats[j].a_ins+stats[j].ab_ins)/(double)(stats[j].a_del+stats[j].ab_del)); 
            double insdel_a =  stats[j].a_del==0 ? -1 : ((double)stats[j].a_ins/(double)stats[j].a_del); 
            double insdel_ab =  stats[j].ab_del==0 ? -1 : ((double)stats[j].ab_ins/(double)stats[j].ab_del); 
            double insdel_b =  stats[j].b_del==0 ? -1 : ((double)stats[j].b_ins/(double)stats[j].b_del); 
            double insdel_b_ab =  (stats[j].b_del+stats[j].ab_del)==0 ? -1 : ((double)(stats[j].b_ins+stats[j].ab_ins)/(double)(stats[j].b_del+stats[j].ab_del)); 
         
            uint32_t totala = stats[j].a+stats[j].ab;
            double ab_of_a = totala==0 ? -1 : ((double)stats[j].ab/(double)(totala));
             
            uint32_t totalb = stats[j].b+stats[j].ab;
            double ab_of_b = totalb==0 ? -1 : ((double)stats[j].ab/(double)(totalb));            
        
            if (j==1)
            {
                std::cout << std::setprecision(3);
                std::cout << "\n";
                std::cout << dataset_labels[0] << "\n";
                std::cout << "variants: " << totala << "\n";
                std::cout << "ins/del ratio: " << insdel_a_ab << "\n"; 
                std::cout << "FS Proportion: " << (float)fs/(float)(fs+nfs) << " (" << fs << "," << nfs << ")\n\n";
                std::cout << "FS Proportion (Rare): " << (float)rare_fs/(float)(rare_fs+rare_nfs) << " (" << rare_fs << "," << rare_nfs << ")\n\n";
                std::cout << "FS Proportion (Common): " << (float)common_fs/(float)(common_fs+common_nfs) << " (" << common_fs << "," << common_nfs << ")\n\n";
            
            }    
            
            std::cout << std::setprecision(3);
            std::cout << dataset_labels[j] << " (" << totalb << ") " << "[" << insdel_b_ab << "]\n";
            std::cout << "TP\t" << ab_of_b << "(" << stats[j].ab << "/" << totalb << ") " << "[" << insdel_ab << "," << insdel_b << "] \n";
            std::cout << "FP\t" << ab_of_a << "(" << stats[j].ab << "/" << totala << ") " << "[" << insdel_ab << "," << insdel_a << "] \n\n";
        }
    };
	
 	~Igor()
    {
    };
    
    private:
        
    /** 
     * Populate GENCODE tree 
     */
    void populate_gencode_tree(const char* interval)
    {
        std::string line;
    	std::vector<std::string> fields;
    	std::vector<std::string> alleles;
       
        std::clog << "Populating GENCODE tree ... \n";
        
        std::vector<GTFRecord*> exons; 
        std::vector<Interval*> intervals;
        uint32_t noRecords = 0;
    	
        tbx_t *tbx;
	    hts_itr_t *iter;
	    htsFile* gffgz = hts_open(gencode_gtf_file.c_str(), "r");
        if (!(tbx=tbx_index_load(gencode_gtf_file.c_str())))
        {
    		std::cerr << "Fail to load the VCF GZ index\n";
    		abort();
    	}
   	    	   
   	    std::stringstream region;
	    kstring_t s;
        s.s = 0; s.l = s.m = 0;
      
        if (!(iter = tbx_itr_querys(tbx, interval))) 
        {
            std::cerr << "Cannot access interval " << interval << " in GENCODE file\n";
                exit(1);
		}
    
	    while (tbx_itr_next(gffgz, tbx, iter, &s) >= 0)
        {
            
            //populate interval trees with reference sets
            //gene transfer format
            //1	chr1
            //2	HAVANA - source
            //3	exon   - feature
            //4	13221  - start
            //5	14409  - end
            //6	.      - score
            //7	+      - strand
            //8	.      - frame
            //9	gene_id "ENSG - attribute
            
            //std::cerr << s.s << "\n";
            
            std::string line = std::string(s.s);
        	
    		if(line.at(0)!='#')
        	{   
        	 	split(fields, "\t", line);
                std::string chrom = fields[0].substr(3);	
    		    
    		    
                std::string feature = fields[2];
                    
                //if (feature!="CDS" && feature!="exon")
                if (feature!="CDS")    
                {
                    continue;
                }
                        
                //create tree for chromosome
                if(!exists(GENCODE, chrom))
                {
                    std::clog << "Creating tree for chromosome " << chrom << "\n";
                    GENCODE[chrom] = new IntervalTree();
                }
                
                //process fields			
                std::map<std::string, std::string> attrib_map;
                splitGTFAttributeFields(attrib_map, fields[8]);
    
//                   	std::cerr << line << "\n"; 
//                	for (std::map<std::string, std::string>::iterator i=attrib_map.begin(); i!=attrib_map.end() ; ++i)
//                	{
//                	    std::cerr << i->first << " => " << i->second << "\n";
//                	}
                //exit(1);    	
            	uint32_t start = atoi(fields[3].c_str());
    			uint32_t end = atoi(fields[4].c_str());
    			char strand = fields[6].at(0);	
                std::string gene = attrib_map["gene_id"];
    		    int32_t frame = fields[7]=="." ? -1 : atoi(fields[7].c_str());
                int32_t exonNo = attrib_map.end()==attrib_map.find("exon_number") ? -1 : atoi(attrib_map["exon_number"].c_str());
            	uint32_t level = atoi(attrib_map["level"].c_str());
            	std::string attrib = fields[8];
            	                  
    			if (feature=="CDS")
    			{
    			    if (attrib_map["gene_type"]!="protein_coding")
                    {   
    			        continue;
    			    }
    			   
                    GTFRecord* gtfRecord = new GTFRecord(chrom, start, end, strand,
                                                         gene, feature, frame, exonNo, 
                                                         true, true, 
                                                         true, true,
                                                         level, attrib);
                    
                    GENCODE[chrom]->insert(gtfRecord);
    			}
                ++noRecords;
            }
        }
        
        
        	
    	std::clog << " ... completed\n";
    	if (GENCODE.size()==0)
    	{
    	    std::cerr << "No reference GENCODE features!\n";
            exit(1); 
    	}
    	
    	//std::cerr << "start validation\n"; 
        //GENCODE["X"]->validate();
        //G = GENCODE["chr20"];
        //G->validate();
        //std::cerr << "end validation\n";
        //std::cerr << "height : " << GENCODE["X"]->height << "\n";     
        //std::cerr << "size : " << GENCODE["X"]->size() << "\n";     
//        std::cerr << "height : " << G->height << "\n";     
//        std::cerr << "size : " << G->size() << "\n";  
        std::cerr << "size : " << noRecords << "\n";     


        if(s.s) free(s.s); 
    }
    bool exists(std::map<std::string, IntervalTree*>& map, const std::string& key)
    {
        return map.end()!=map.find(key);
    }
    
    /**
     *Splits a line into a map - PERL style
     */
    void splitGTFAttributeFields(std::map<std::string, std::string>& map, std::string& str)
    {
    	map.clear();	
    	const char* tempStr = str.c_str();
    	int32_t i=0, lastIndex = str.size()-1;
    	std::string key;
        std::string val;
    	std::stringstream token;
    	
    	if (lastIndex<0) return;
    	
    	while (i<=lastIndex)
    	{
    	    //read next character
    		if(tempStr[i]!=';' && tempStr[i]!=' ') 
    		{
    			token << tempStr[i];
    		}
    
            //store key-value pair 
    		if (i==lastIndex || (tempStr[i]==';' && tempStr[i+1]==' ')) 
    		{
    		    val = token.str();
    		    if (val.at(0)=='"')
    		    {    
        			val = token.str();
        			val = val.substr(1,val.size()-2);
        		}
        		map[key] = val;
    		    
    			token.str("");
    			key.clear();
    			
    			if (i!=lastIndex)
    			    ++i;
    		}		
    		
    		//store key
    		if (tempStr[i]==' ') 
    		{
    			key = token.str();
    			token.str("");
    		}
    		
    		++i;
    	} 
    };      
};

}

void profile_indels(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_indels();
    igor.print_stats();
}
 