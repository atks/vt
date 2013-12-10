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

#include "genotype.h"

namespace 
{
    
class PileupRecord
{
    public:
    bcf_hdr_t* ivcf_hdr;
    bcf1_t* ivcf_rec;
    std::map<std::string, int> read_ids;
    uint32_t read_no; 
     
    std::stringstream rqs;
 	std::stringstream aqs;
 	std::stringstream mqs;
 	std::stringstream tps;
 	std::stringstream ss;    
    
    std::stringstream log;
    
    PileupRecord()
    {
        ivcf_rec = bcf_init1();
        read_ids.clear();
        rqs.str("");
		aqs.str("");
		mqs.str("");
        tps.str("");
        ss.str("");
        log.str("");
    }
    
    ~PileupRecord()
    {
        bcf_destroy1(ivcf_rec);
    }

    void print(std::ofstream *OUT_GLF)
    {
        const char* chrom = bcf_get_chrom(ivcf_hdr, ivcf_rec);
        uint32_t pos0 = bcf_get_pos0(ivcf_rec);
        char* ref = bcf_get_ref(ivcf_rec);
        char* alt = bcf_get_alt(ivcf_rec, 1);
       
     	ss << chrom << "\t" 
		   << (pos0+1) << "\t" 
		   << ref  << "\t" 
		   << alt << "\t";

		if (read_no)
		{				
			ss << read_no << "\t"
		       << rqs.str() << aqs.str() << "\t"
		       << mqs.str() << "\t"
		       << tps.str();
		}
		else
		{
			ss << "0\t.\t.\t.";
		}
		
		*OUT_GLF << ss.str() << "\n";
		
		if (log.str()!="")
		{
		    std::clog << log.str();  
        }
    }

    void clear()
    {
        read_ids.clear();
        read_no = 0; 
        rqs.str("");
		aqs.str("");
		mqs.str("");
        tps.str("");
        ss.str("");
        log.str("");
    }
}; 
  
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
    std::vector<GenomeInterval> intervals;
    bool debug;    
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *vodr;
    BAMOrderedReader *sodr;
    BCFOrderedWriter *vodw;
    bcf1_t *v;
    bam1_t *s;

    /////////
    //stats//
    /////////
    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;
    
    /////////
    //tools//
    /////////
    LogTool lt;
    LHMM lhmm_ref, lhmm_alt;
    
    //////////////////////////////////////////
    //share objects for genotype_<snp|indel>//
    //////////////////////////////////////////
    std::stringstream region;
    hts_itr_t *iter, *bam_iter, *vcf_iter;
    std::map<std::string, int> read_ids;	
         
    double log_p_rr;
    double log_p_ra;
    double log_p_aa;        
    
    std::stringstream rqs;
 	std::stringstream aqs;
 	std::stringstream mqs;
 	std::stringstream tps;

    //////////////////////////////////////////
    //share objects for genotype_<snp|indel>//
    //////////////////////////////////////////
    kstring_t str;
    kstring_t readseq;
    kstring_t readqual;
    kstring_t refProbe; 
    kstring_t altProbe; 
    kstring_t plenstring;
    
    Igor(int argc, char ** argv)
    :ivcf(NULL),
    ivcfgz(NULL)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
    	try 
    	{
    		std::string desc = 
"Generates SNP and Indel pileup for each sample\n\
e.g. vt pileup -i probes.sites.vcf -o out.pileup -b NA19130.bam -s NA19130\n\n";
       			
       		version = "0.57";
       		TCLAP::CmdLine cmd(desc, ' ', version);
  		    VTOutput my; cmd.setOutput(&my);
    		TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file("i", "i", "Input Candidate VCF file", true, "", "file", cmd);
    		TCLAP::ValueArg<std::string> arg_input_sam_file("b", "i", "Input BAM file", true, "", "str", cmd);
    		TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "Output PILEUP file", false, "-", "file", cmd);
    		TCLAP::ValueArg<std::string> arg_sample_id("s", "s", "Sample ID", true, "", "str", cmd);
    		TCLAP::SwitchArg arg_debug("d", "d", "Debug alignments", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
 		                
            cmd.parse(argc, argv);
        
            input_vcf_file = arg_input_vcf_file.getValue();
    		input_sam_file = arg_input_sam_file.getValue();
    		output_vcf_file = arg_output_vcf_file.getValue();
    		sample_id = arg_sample_id.getValue();	
    		parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
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
        vodr = new BCFOrderedReader(input_vcf_file, intervals);       
        sodr = new BAMOrderedReader(input_sam_file, intervals);       
        vodw = new BCFOrderedWriter(output_vcf_file, intervals);
        
        ////////////////////////
        //stats initialization//
        ////////////////////////    
        no_snps_genotyped = 0;
        no_indels_genotyped = 0;    
        
        str = {0,0,0};
        readseq = {0,0,0};
        readqual = {0,0,0};	
        refProbe = {0,0,0};		
        altProbe = {0,0,0};	
        plenstring = {0,0,0};
    };
 	
 	~Igor()
    {
    };

    void print_options()
    {
        std::clog << "pileup v" << version << "\n\n";
    	
	    std::clog << "Options: Input VCF File      " << input_vcf_file << "\n";
	    std::clog << "         Input BAM File      " << input_sam_file << "\n";
	    std::clog << "         Output VCF File  " << output_vcf_file << "\n";
	    std::clog << "         Sample ID           " << sample_id << "\n\n";
    }

    void print_stats()
    {
	    std::clog << "Stats: SNPs genotyped     " << no_snps_genotyped << "\n";
	    std::clog << "       Indels genotyped   " << no_indels_genotyped << "\n\n";
    }

    /**
     *The most basic model for GLs.
     */
    double log10Emission(char readBase, char probeBase, uint32_t qual)
    {
        double e = lt.pl2prob(qual);
        
    	if (readBase=='N' || probeBase=='N')
    	{
    		return 0;
    	}
    	
   		return readBase!=probeBase ? log10(e/3) : log10(1-e);
    };
    
    ////////////////////////////////////////
    //Buffer for processing pileup records//
    ////////////////////////////////////////
    //pileup records to be processed
    std::list<PileupRecord*> precs;
    //pool of unused pileup records
    std::list<PileupRecord*> pool;
    
    uint32_t bref, vref;
    uint32_t bstart, bend, vpos;
    
    std::stringstream ss;
    /**
     * Adds pileup records till the first record after end. 
     * Coordinates all 1 based.
     */
    bool add_prec(int32_t epos1)
    {   
        //add records only if the new record overlaps with read
        if (precs.size()!=0)
        {
            bcf1_t *ivcf_rec = precs.back()->ivcf_rec;
            int32_t vpos1 = bcf_get_pos1(ivcf_rec);
            
            if (vpos1>epos1)
            {
                return false;
            }    
        }
        
        kstring_t s;
        s.s = 0; s.l = s.m = 0;
        bool added_record = false;
        while (tbx_itr_next(ivcfgz, ivcf_tbx, vcf_iter, &s) >= 0)
        {   
            PileupRecord *p = NULL;
            if (pool.size())
            {
                p = pool.front();
                pool.pop_front();
            }
            else
            {
                p = new PileupRecord();
            }
            
            p->clear();
            p->ivcf_hdr = ivcf_hdr;
            
            vcf_parse1(&s, ivcf_hdr, p->ivcf_rec);
            bcf_unpack(p->ivcf_rec, BCF_UN_STR);
            bcf_set_variant_types(p->ivcf_rec);
            
            precs.push_back(p); 
            added_record = true;
        
            if ((bcf_get_pos1(p->ivcf_rec))>epos1)
                break;
        }  
        
        return added_record;      
    }

    int32_t read_is_before_first_vcf_record(bam1_t *isam_rec)
    {
        bcf1_t *ivcf_rec = (*(precs.begin()))->ivcf_rec;
        int32_t vpos1 = bcf_get_pos1(ivcf_rec);
        int32_t epos1 = bam_get_end_pos1(isam_rec);
  
        return (epos1+100000)<vpos1 ? vpos1 : 0;
    }

    /**
     * Genotype all variants in the VCF file.
     */
    void pileup1()
    {
        int32_t nseqs;
        int32_t* lens;
        const char** seqs;
        bcf_hdr_get_seqs_and_lens(ivcf_hdr, seqs, lens, &nseqs);
        
        //iterate through chromosomes in VCF file
		for (int32_t i = 0; i < nseqs; ++i) 
		{
            std::cerr << "VCF id:" << i << " SQ:" << seqs[i]
			          << " LN:" << lens[i] << "\n";

            //choose records from chromosome i
            if (vcf_ftype==IS_BCF)
            {
    			if (!(vcf_iter = bcf_itr_queryi(ivcf_idx, i, 1, lens[i]))) 
    			{
    				fprintf(stderr, "[E::] fail to parse region VCF idx\n");
    				abort();
    			}
    		}
    		else if (vcf_ftype==IS_VCF_GZ)
		    {
		        std::stringstream region;
		        region << seqs[i] << ":" << 1 << "-" << lens[i];
		        if (!(vcf_iter = tbx_itr_querys(ivcf_tbx, region.str().c_str()))) 
		        {
    				continue;
    			}
		    }

		    int32_t chromid = bam_name2id(isam_hdr,  seqs[i]);
		    if (chromid==-1)
		    {
		        //print out all VCF records as empty
		        continue;
		    }    
		    int32_t clen = isam_hdr->target_len[chromid];
		    
		    std::cerr << "BAM id:" << chromid << " SQ:" << seqs[i] << " LN:" << clen << "\n";
		    
            //choose reads from chromosome i
    	    if (!(bam_iter = bam_itr_queryi(isam_idx, chromid, 1, clen)))
            {
               std::cerr << "fail to parse regions\n";
               abort();    
            }
            
            int32_t count = 0;
            //WARP:
            //std::cerr << "warped\n";    
            while(bam_itr_next((BGZF*)isam->fp, bam_iter, srec)>=0)
            {               
                //read in vcf records till it passes the end of this read
                int32_t spos1 = bam_get_pos1(srec);
                int32_t epos1 = bam_get_end_pos1(srec);
                add_prec(epos1);
                
                //std::cerr << "EReadPos1: " << epos1 << "\n";   
                                                
//                int32_t npos1 = 0;
//                if ( (npos1 = read_is_before_first_vcf_record(srec)))
//                {   
//                    
//                   // std::cerr << "\tEReadPos1: " << epos1 << " SVCFPos1: "  << npos1 << "\n";   
//                    
//                    //continue;
//                    if ((bam_iter = bam_itr_queryi(isam_idx, chromid, npos1+1, clen)))
//                    {
//                        //should we do something here?
//                        //::cerr << "\t\t jumping forward " << chromid << ":" << npos1 << "(" << clen << ")\n";   
//                        goto WARP;
//                        continue;
//                    }
//                }
////             
                //iterate through VCF records
                std::list<PileupRecord*>::iterator prec_iter = precs.begin();
                //int32_t count = 0;
                while (prec_iter!=precs.end())
                {
                    //std::cerr << count++ << "\n";
                    
                    PileupRecord* prec = *prec_iter;
                    bcf1_t *ivcf_rec = prec->ivcf_rec;
                    int32_t vpos1 = bcf_get_pos1(ivcf_rec);
                    
                    //print out records and remove
                    if (vpos1 < spos1)
                    {
                        if (bcf_get_var_type(ivcf_rec) == VCF_SNP) 
                        {
                            prec->print(OUT_GLF);
                            ++no_snps_genotyped;
                        }
                        else if (bcf_get_var_type(ivcf_rec) == VCF_INDEL || bcf_get_var_type(ivcf_rec) == VCF_OTHER) 
                        {
                            prec->print(OUT_GLF);
                            ++no_indels_genotyped;
                        }
                        
                        //move record to pool
                        prec->clear();
                        pool.push_front(prec);
                        prec_iter = precs.erase(prec_iter);
                        
                    }
                    else if (vpos1>=spos1 && vpos1<=epos1)
                    {
                        if (bcf_get_var_type(ivcf_rec) == VCF_SNP) 
                        {
                            genotype_snp1(*prec, srec);
                        }
                        else if (bcf_get_var_type(ivcf_rec) == VCF_INDEL || bcf_get_var_type(ivcf_rec) == VCF_OTHER) 
                        {
                            genotype_indel1(*prec, srec);
                        }
                        else
                        {
                            std::cerr << "Only SNPs and Indels/Complex Substitutions supported for this pileup, skipping this variant\n";
                            std::cerr << "no alleles: " << ivcf_rec->n_allele << "\n";          
                            for (uint32_t i=0; i<ivcf_rec->n_allele; ++i)
                            {
                                std::cerr << "\t" << ivcf_rec->d.allele[i] << "\n";
                            } 
                        }

                        ++prec_iter;
                    }
                    else
                    {
                        break;
                        //till we meet again
                    }
                }
                
                ++count;
            }//end bam iteration for this chromosome
            
            //add_prec(clen);
            
            //flush out remaining records
            std::list<PileupRecord*>::iterator prec_iter = precs.begin();
            while (prec_iter!=precs.end())
            {
                PileupRecord* prec = *prec_iter;
                bcf1_t *ivcf_rec = prec->ivcf_rec;
                
                if (bcf_get_var_type(ivcf_rec) == VCF_SNP) 
                {
                    prec->print(OUT_GLF);
                    ++no_snps_genotyped;
                }
                else if (bcf_get_var_type(ivcf_rec) == VCF_INDEL || bcf_get_var_type(ivcf_rec) == VCF_OTHER) 
                {
                    prec->print(OUT_GLF);
                    ++no_indels_genotyped;
                }
                
                //move record to pool
                prec->clear();
                pool.push_front(prec);
                prec_iter = precs.erase(prec_iter);                    
            }
            
            kstring_t s;
            s.s = 0; s.l = s.m = 0;
            while (tbx_itr_next(ivcfgz, ivcf_tbx, vcf_iter, &s) >= 0)
            {   
                vcf_parse1(&s, ivcf_hdr, ivcf_rec);
                bcf_unpack(ivcf_rec, BCF_UN_STR);
                bcf_set_variant_types(ivcf_rec);
    
                *OUT_GLF << bcf_get_chrom(ivcf_hdr, ivcf_rec) << "\t" 
        		         << bcf_get_pos1(ivcf_rec) << "\t" 
        		         << bcf_get_ref(ivcf_rec)  << "\t" 
        		         << bcf_get_alt(ivcf_rec, 1) << "\t"
        	             << "0\t.\t.\t.";
            }  
                        
            std::cerr << "chrom " << i << ":" << count << " reads\n";            
        }//end chromosome
    }

    void genotype_snp1(PileupRecord &prec, bam1_t* srec)
    {    
        //maximum depth cap
        if (prec.read_no>255) 
        {
            return;
        }
        
        //has secondary alignment or fail QC or is duplicate or is unmapped
        if (bam_get_flag(srec) & 0x0704) 
        {
            return;
        }
    
        //make this setable
        //ignore poor quality mappings
        if(bam_get_mapq(srec)<20)
        {
            return;
        }
        
        bcf1_t *ivcf_rec = prec.ivcf_rec;
        const char* chrom = bcf_get_chrom(ivcf_hdr, ivcf_rec);
        uint32_t pos0 = bcf_get_pos0(ivcf_rec);
        char ref = bcf_get_snp_ref(ivcf_rec);
        char alt = bcf_get_snp_alt(ivcf_rec);
        char* read_name = bam_get_qname(srec);
       
        //std::cerr << chrom << ":" << pos << ":" << ref << ":" << alt << "\n";
        
        //handle mate pairs
		if(prec.read_ids.find(read_name)==prec.read_ids.end())
		{
		    prec.read_ids[read_name] = 1;
		}
		else
	    {
	        prec.read_ids.erase(read_name);
	        return;
	    }
	    
	    //get base and qual
		char base, qual; int32_t rpos;
        bam_get_base_and_qual(srec, ivcf_rec->pos, base, qual, rpos);         
        
        //strand
        char strand = bam_is_rev(srec) ? 'R' : 'F';

        //map qual
        uint32_t mapQual = bam_get_mapq(srec);
       
        //fail to find the expected mapped base on the read
        if (rpos==BAM_READ_INDEX_NA)
        {
            return;
        }

        //ignore ambiguous bases
        if (base=='N')
        {
            return;
        }
                           
		//////////////////////////////////////////////////
		//perform GL computation here, uses map alignments
		//////////////////////////////////////////////////
        std::string refCigar = "";
        std::string altCigar = "";
            
        //compute genotype likelihood using alignment coordinates
        double refllk = log10Emission(ref, base, qual);
        double altllk = log10Emission(alt, base, qual);
        
        ////////////////////////////////
        //aggregate genotype likelihoods
        ////////////////////////////////
        uint32_t baseqr = 0, baseqa = 0;
        
        if (refllk>altllk)
        {
            baseqa = round(-10*(altllk-refllk));
        }
        else if (refllk<altllk)
        {
            baseqr = round(-10*(refllk-altllk));
        }

        if (debug)
        {
            kstring_t cigar;
            cigar.l = cigar.m = 0; cigar.s = 0;
            bam_get_cigar_string(srec, &cigar);

            kstring_t aligned_read;
        	aligned_read.l = aligned_read.m = 0;
        	aligned_read.s = 0;

        	kstring_t aligned_qual;
        	aligned_qual.l = aligned_qual.m = 0;
        	aligned_qual.s = 0;

        	kstring_t expanded_cigar;
        	expanded_cigar.l = expanded_cigar.m = 0;
        	expanded_cigar.s = 0;

        	kstring_t annotations;
        	annotations.l = annotations.m = 0;
        	annotations.s = 0;

            bam_get_seq_string(srec, &readseq);
            bam_get_qual_string(srec, &readqual, readseq.s);
            
            char allele = strand=='F' ? 'N' : 'n';
            if (refllk>altllk)
            {
            	allele = strand=='F' ? 'R' : 'r';
            }
            else if (refllk<altllk)
            {
            	allele = strand=='F' ? 'A' : 'a';
            }   
            
            std::string pad = "\t";
	        prec.log << pad << "==================\n";
	        prec.log << pad <<  (prec.read_no+1) << ") " << chrom << ":" << (pos0+1) << ":" << ref << ":" << alt << "\n";
            prec.log << pad << bam_get_qname(srec) << "\n";
            prec.log << pad << "==================\n";
            print_read_alignment(readseq, readqual, cigar, aligned_read, aligned_qual, expanded_cigar, annotations, rpos);

            prec.log << pad << "read  " << aligned_read.s  << "\n";
            prec.log << pad << "qual  " << aligned_qual.s  << "\n";
            prec.log << pad << "cigar " << expanded_cigar.s  << "\n";
            prec.log << pad << "anno  " << annotations.s  << "\n";
            prec.log << pad << "==================\n";
            prec.log << pad << "base   " << base  << "\n";
            prec.log << pad << "qual   " << (int32_t)(qual-33)  << "\n";
            prec.log << pad << "rpos   " << rpos  << "\n";
			prec.log << pad << "refllk " << refllk << "\n";
			prec.log << pad << "altllk " << altllk << "\n";
			prec.log << pad << "allele " << allele << "\n\n";
		}
        
        prec.rqs << (uint8_t) (baseqr>93 ? 126 : baseqr+33);
		prec.aqs << (uint8_t) (baseqa>93 ? 126 : baseqa+33);
		prec.mqs << (uint8_t) (mapQual>93 ? 126 : mapQual+33);
		int32_t tp =  pack_tp(!bam_is_rev(srec), bam_is_fread1(srec), rpos, bam_get_l_qseq(srec));
	    prec.tps << std::hex << (tp<16?"0":"") << tp << std::dec ;
		
		++prec.read_no;         		
    }  
    
    void genotype_indel1(PileupRecord &prec, bam1_t* srec)
    {
        //maximum depth cap
        if (prec.read_no>255) 
        {
            return;
        }
        
        //has secondary alignment or fail QC or is duplicate or is unmapped
        if (bam_get_flag(srec) & 0x0704) 
        {
            return;
        }
    
        //ignore poor quality mappings
        if (bam_get_mapq(srec)<13)
        {
            return;
        }
        
        char* read_name = bam_get_qname(srec);
        if(prec.read_ids.find(read_name)==prec.read_ids.end())
		{
		    //assign id to reads to ease checking of overlapping reads
			prec.read_ids[read_name] = 1;
		}
		else
	    {
	        prec.read_ids.erase(read_name);
	        return;
	    }

        bcf1_t *ivcf_rec = prec.ivcf_rec;
        
        const char* chrom = bcf_get_chrom(ivcf_hdr, ivcf_rec);
        uint32_t pos0 = bcf_get_pos0(ivcf_rec);
        char* ref = bcf_get_indel_ref(ivcf_rec);
        char* alt = bcf_get_indel_alt(ivcf_rec);
        
        //std::cerr << chrom << ":" << pos << ":" << ref << ":" << alt << "\n";
        
        //get probe
        std::vector<std::string> candidate_alleles;
        candidate_alleles.push_back(ref); 
        candidate_alleles.push_back(alt);
        std::vector<std::string> probes;
        	
        char* refProbe = bcf_get_info1(ivcf_hdr, ivcf_rec, "REFPROBE");
        assert(refProbe!=NULL);
        char* altProbe = bcf_get_info1(ivcf_hdr, ivcf_rec, "ALTPROBE");
		assert(altProbe!=NULL);
		char* plenstring = bcf_get_info1(ivcf_hdr, ivcf_rec, "PLEN");
		assert(plenstring!=NULL);
		uint32_t plen = boost::lexical_cast<uint32_t>(plenstring);
		
//		std::cerr << "\t" << refProbe << "\n";
//		std::cerr << "\t" << altProbe << "\n";
		
		int32_t probeLength = strlen(refProbe);
        int32_t variantLengthDifference = (int32_t)strlen(alt)-(int32_t)strlen(ref);
       
        //read name
        
        //sequence
        readseq.l=0;
        bam_get_seq_string(srec, &readseq);

        //qual
        bam_get_qual_string(srec, &readqual, readseq.s);
        
        //map qual
        uint32_t mapQual = bam_get_mapq(srec);
 
 		std::string refCigar = "";
        std::string altCigar = "";
            
        double refllk = 0, altllk = 0;
        
        if (1)
        {
    	   	lhmm_ref.align(refllk, refProbe, readseq.s, readqual.s);
            lhmm_ref.computeLogLikelihood(refllk, lhmm_ref.getPath(), readqual.s);
        }
        else
        {
            //if all Ms
            //if (!attempt_to_quick_align(lhmm_ref, refProbe, readseq.s, refllk, readqual.s))
            {
                lhmm_ref.align(refllk, refProbe, readseq.s, readqual.s);
                lhmm_ref.computeLogLikelihood(refllk, lhmm_ref.getPath(), readqual.s);
            }
    	}
       
        lhmm_alt.align(altllk, altProbe, readseq.s, readqual.s);
        lhmm_alt.computeLogLikelihood(altllk, lhmm_alt.getPath(), readqual.s);
        
        std::string pad = "\t";
        if (debug)
        {
	        prec.log << pad << "==================\n"; 
	        prec.log << pad <<  (prec.read_no+1) << ") " << chrom << ":" << (pos0+1) << ":" << ref << ":" << alt << ":" << variantLengthDifference << "\n";
            prec.log << pad <<  read_name << "\n";
            prec.log << pad << "==================\n";
            prec.log << pad << "ref probe     " << refProbe << " (" << plen << "/" << probeLength << ")\n";
            prec.log << pad << "read sequence " << readseq.s  << "\n";    
            prec.log << pad << "==================\n";
            lhmm_ref.printAlignment(pad, prec.log);
            prec.log << pad << "ref llk: " << refllk << "\n";
	        prec.log << pad << "expected indel location: " << plen+1 << "\n";
	        prec.log << pad << "==================\n";
	        prec.log << pad << "==================\n";
            prec.log << pad << "alt probe     " << altProbe << " (" << plen << ")\n";
            prec.log << pad << "read sequence " << readseq.s  << "\n";    
            prec.log << pad << "==================\n";
            lhmm_alt.printAlignment(pad, prec.log);
            prec.log << pad << "alt llk: " << altllk << "\n";
            prec.log << pad << "expected indel location: " << plen+1 << "\n";
            prec.log << pad << "==================\n\n";
        }
				           
        //check if the Insertions and Deletions are at the expected places.
        //deletion
        uint32_t ref_rpos=0, alt_rpos=0;
        if (variantLengthDifference<0)
        {
            //deletion
            if (!deletion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
                 insertion_start_exists(lhmm_alt, plen+1, alt_rpos))
            {
                if (refllk<altllk)
                {
                    swap(refllk, altllk);
                }
                
                //ref allele
                if (debug)
                    std::cerr << pad << "DEL: REF ALLELE\n"; 
            }
            else if (deletion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
                    !insertion_start_exists(lhmm_alt, plen+1, alt_rpos)) 
            {
                if (refllk>altllk)
                {
                    swap(refllk, altllk);
                }    
                
                //alt allele
                if (debug)
                    std::cerr << pad << "DEL: ALT ALLELE\n"; 
            }
            else
            {
                refllk = 0;
                altllk = 0;
                if (debug)
                    std::cerr << pad << "Deletion not at expected location, set as ambiguous\n"; 
            }
        }
        else if (variantLengthDifference>0)
        { 
            //insertion
            if (!insertion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
                 deletion_start_exists(lhmm_alt, plen+1, alt_rpos))
            {
                if (refllk<altllk)
                {
                    swap(refllk, altllk);
                }    
                
                //ref allele
                if (debug)
                    std::cerr << pad << "INS: REF ALLELE\n"; 
            }
            else if (insertion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
                    !deletion_start_exists(lhmm_alt, plen+1, alt_rpos)) 
            {
                if (refllk>altllk)
                {
                    swap(refllk, altllk);
                } 
                
                //alt allele
                if (debug)
                    std::cerr << pad << "INS: ALT ALLELE\n"; 
            }
            else
            {
                refllk = 0;
                altllk = 0;
                if (debug)
                    std::cerr << pad << "Insertion not at expected location "  << plen+1 << ", set as ambiguous\n"; 
            }
        }
 
        if (debug)
        {    
	        std::cerr << pad << "++++++++++++++++++\n";
            std::cerr << pad << "reflk " << refllk << "\n";
		    std::cerr << pad << "altlk " << altllk << "\n";
        }

        uint32_t baseqr = 0, baseqa = 0;
        //char allele = 'N';
        uint32_t rpos = 0;
        if (refllk>altllk)
        {
        	//allele = 'R';            
            rpos = alt_rpos;
            baseqa = round(-10*(altllk-refllk));
        }
        else if (refllk<altllk)
        {
        	//allele = 'A';            	    
      	    rpos = ref_rpos;
      	    baseqr = round(-10*(refllk-altllk));
        }
        
        prec.rqs << (uint8_t)(baseqr>93? 126 : baseqr+33);
		prec.aqs << (uint8_t)(baseqa>93? 126 : baseqa+33);
		prec.mqs << (uint8_t)(mapQual>93? 126 : mapQual+33);
		int32_t tp =  pack_tp(!bam_is_rev(srec), bam_is_fread1(srec), rpos, bam_get_l_qseq(srec));
	    prec.tps << std::hex << (tp<16?"0":"") << tp << std::dec ;
					
		++prec.read_no;   

        free(refProbe);
        free(altProbe);
        free(plenstring);
    }
    
//    void pileup()
//    {
//        while (vcf_read1(ivcf, ivcf_hdr, ivcf_rec)==0) 
//        {
//            bcf_unpack(ivcf_rec, BCF_UN_STR); 
//            bcf_set_variant_types(ivcf_rec);
//            
//            if (ivcf_rec->d.var_type == VCF_SNP) 
//            {
//                genotype_snp();
//                ++no_snps_genotyped;
//            }
//            else if (ivcf_rec->d.var_type == VCF_INDEL || ivcf_rec->d.var_type == VCF_OTHER) 
//            {
//                genotype_indel();
//                ++no_indels_genotyped;
//            }
//            else
//            {
//                std::cerr << "Only SNPs and Indels/Complex Substitutions supported for this pileup, skipping this variant\n";
//                std::cerr << "no alleles: " << ivcf_rec->n_allele << "\n";          
//                for (uint32_t i=0; i<ivcf_rec->n_allele; ++i)
//                {
//                    std::cerr << "\t" << ivcf_rec->d.allele[i] << "\n";
//                } 
//            }
//        }
//    }

    //expands cigar string
    void generateCigarString(kstring_t* expanded_cigar, kstring_t* cigar)
    {
        expanded_cigar->l = 0;
        int32_t i=0, lastIndex = cigar->l-1;
        std::stringstream token;

        if (lastIndex<0)
        {
            return;
        }
        char c;
        bool seenM = false;

        while (i<=lastIndex)
        {
        	c = cigar->s[i];

        	//captures the count
            if (c<'A')
            {
            	token << c;
            }

            if (c>'A' ||
                i==lastIndex)
            {
                uint32_t count = boost::lexical_cast<uint32_t>(token.str());

                //it is possible for I's to be observed before the first M's in the cigar string
                //in this case, we treat them as 'S'
                if (!seenM)
                {
                    if (c=='I')
                    {
                        c = 'S';
                    }
                    else if (c=='M')
                    {
                        seenM = true;
                    }
                }

                for (uint32_t j=0; j<count; ++j)
                    kputc_(c, expanded_cigar);
                token.str("");
            }

            ++i;
        }
        
        kputc_(0, expanded_cigar);
    };
        
    void print_read_alignment(kstring_t& read, kstring_t& qual, kstring_t& cigar, 
                              kstring_t& aligned_read, kstring_t& aligned_qual, kstring_t& expanded_cigar, kstring_t& annotations,
                              uint32_t rpos)
    {
    	aligned_read.l = 0;
    	aligned_qual.l = 0;
     	expanded_cigar.l = 0;
    	generateCigarString(&expanded_cigar, &cigar);
     	annotations.l = 0;
    	
    	uint32_t j=0;
    	for (uint32_t i=0; i<expanded_cigar.l; ++i)
    	{
    	    char state = expanded_cigar.s[i];
    	    if (state=='M' || state=='S' || state=='I')
    	    {
    	        kputc(read.s[j], &aligned_read);
    	        kputc(qual.s[j], &aligned_qual);
    	        
    	        if (j==rpos)
    	        {    
    	            kputc('^', &annotations);
    	        }
    	        else
    	        {    
    	            kputc(' ', &annotations);
    	        }
    	        ++j;
    	    }
    	    else if (state=='D')
    	    {
                kputc(' ', &aligned_read);
    	        kputc(' ', &aligned_qual);
    	    }
    	}
    }

//    void genotype_snp()
//    {
//        if (ivcf_rec->n_allele==2)
//        {
//            read_ids.clear();	
//                 
//            log_p_rr = 0;
//            log_p_ra = 0;
//            log_p_aa = 0;        
//                            
//        	rqs.str("");
//    		aqs.str("");
//    		mqs.str("");
//            tps.str("");
//    		
//            const char* chrom = bcf_get_chrom(ivcf_hdr, ivcf_rec);
//            uint32_t pos0 = bcf_get_pos0(ivcf_rec);
//            char ref = bcf_get_snp_ref(ivcf_rec);
//            char alt = bcf_get_snp_alt(ivcf_rec);
//            uint32_t read_no = 0;     
//            
//            if (!(iter = bam_itr_queryi(isam_idx, bam_name2id(isam_hdr, bcf_get_chrom(ivcf_hdr, ivcf_rec)), ivcf_rec->pos, ivcf_rec->pos+1)))
//            {
//                std::cerr << "fail to parse regions\n";
//                abort();    
//            }
//            
//            while(bam_itr_next((BGZF*)isam->fp, iter, srec)>=0  && read_no<255)
//            {
//                //has secondary alignment or fail QC or is duplicate or is unmapped
//                if (bam_get_flag(srec) & 0x0704) 
//                {
//                    //std::cerr << "fail flag\n";
//                    continue;
//                }
//            
//                //make this setable
//                if(bam_get_mapq(srec)<20)
//    	        {
//    	            //std::cerr << "low mapq\n";
//    	            //ignore poor quality mappings
//    	            continue;
//    	        }
//               
//                char* read_name = bam_get_qname(srec);
//               
//    	        //std::map<std::string, int> readIDs;	
//    			if(read_ids.find(read_name)==read_ids.end())
//    			{
//    			    //assign id to reads to ease checking of overlapping reads
//    				read_ids[read_name] = read_ids.size();
//    			}
//    			//ignore overlapping paired end
//    			else
//    		    {
//    		        //std::cerr << "Mate pair\n";
//    		        read_ids.erase(read_name);
//    		        continue;
//    		    }
//    		        		        
//                //get base and qual
//        		char base, qual; int32_t rpos;
//                
//    			if (1)
//    			{
//    			    bam_get_base_and_qual(srec, ivcf_rec->pos, base, qual, rpos);         
//                }
//                else
//                {
//                    readseq.l = 0;
//                    readqual.l = 0;
//                    bam_get_base_and_qual_and_read_and_qual(srec, ivcf_rec->pos, base, qual, rpos, &readseq, &readqual);
//                }
//                //strand
//                char strand = bam_is_rev(srec) ? 'R' : 'F';
//
//                //map qual
//                uint32_t mapQual = bam_get_mapq(srec);
//               
//                //fail to find a mapped base on the read
//                if (rpos==BAM_READ_INDEX_NA)
//                {
//                    //std::cerr << "rpos NA\n";
//                    continue;
//                } 
//                
//                //ignore ambiguous bases
//                if (base=='N')
//                {
//                    //std::cerr << "N allele\n";
//                    continue;
//                }
//                                            
//    			//////////////////////////////////////////////////
//    			//perform GL computation here, uses map alignments
//    			//////////////////////////////////////////////////
//                std::string refCigar = "";
//                std::string altCigar = "";
//                    
//     	        //compute genotype likelihood using alignment coordinates
//     	        double refllk = log10Emission(ref, base, qual);
//    	        double altllk = log10Emission(alt, base, qual);
//                
//                ////////////////////////////////
//                //aggregate genotype likelihoods
//                ////////////////////////////////
//                log_p_rr = lt.elog10prod(log_p_rr, refllk);
//                log_p_ra = lt.elog10prod(log_p_ra, lt.elog10prod(-log10(2), lt.elog10sum(refllk, altllk)));
//                log_p_aa = lt.elog10prod(log_p_aa, altllk);    
//    
//                uint32_t baseqr = 0, baseqa = 0;
//                char allele = strand=='F' ? 'N' : 'n';
//                if (refllk>altllk)
//                {
//                	allele = strand=='F' ? 'R' : 'r';
//
//                    baseqa = round(-10*(altllk-refllk));
//                }
//                else if (refllk<altllk)
//                {
//                	allele = strand=='F' ? 'A' : 'a';
//
//              	    baseqr = round(-10*(refllk-altllk));
//                }
//
//                if (debug)
//                {
//                    kstring_t cigar;
//                    cigar.l = cigar.m = 0; cigar.s = 0;
//                    bam_get_cigar_string(srec, &cigar);
//
//                    kstring_t aligned_read;
//                	aligned_read.l = aligned_read.m = 0;
//                	aligned_read.s = 0;
//
//                	kstring_t aligned_qual;
//                	aligned_qual.l = aligned_qual.m = 0;
//                	aligned_qual.s = 0;
//
//                	kstring_t expanded_cigar;
//                	expanded_cigar.l = expanded_cigar.m = 0;
//                	expanded_cigar.s = 0;
//
//                	kstring_t annotations;
//                	annotations.l = annotations.m = 0;
//                	annotations.s = 0;
//
//                    std::string pad = "\t";
//    		        std::cerr << pad << "==================\n";
//    		        std::cerr << pad <<  (read_no+1) << ") " << chrom << ":" << (pos0+1) << ":" << ref << ":" << alt << "\n";
//                    std::cerr << pad << bam_get_qname(srec) << "\n";
//                    std::cerr << pad << "==================\n";
//                    print_read_alignment(readseq, readqual, cigar, aligned_read, aligned_qual, expanded_cigar, annotations, rpos);
//
//                    std::cerr << pad << "read  " << aligned_read.s  << "\n";
//                    std::cerr << pad << "qual  " << aligned_qual.s  << "\n";
//                    std::cerr << pad << "cigar " << expanded_cigar.s  << "\n";
//                    std::cerr << pad << "anno  " << annotations.s  << "\n";
//                    std::cerr << pad << "==================\n";
//                    std::cerr << pad << "base   " << base  << "\n";
//                    std::cerr << pad << "qual   " << (int32_t)(qual-33)  << "\n";
//                    std::cerr << pad << "rpos   " << rpos  << "\n";
//        			std::cerr << pad << "refllk " << refllk << "\n";
//        			std::cerr << pad << "altllk " << altllk << "\n";
//        			std::cerr << pad << "allele " << allele << "\n\n";
//        		}
//            
//                rqs << (uint8_t)(baseqr>93? 126 : baseqr+33);
//    			aqs << (uint8_t)(baseqa>93? 126 : baseqa+33);
//    			mqs << (uint8_t)(mapQual>93? 126 : mapQual+33);
//    			int32_t tp =  pack_tp(!bam_is_rev(srec), bam_is_fread1(srec), rpos+1, bam_get_l_qseq(srec));
//    		    tps << std::hex << (tp<16?"0":"") << tp << std::dec ;
//    				
//    			++read_no;     
//            }
//        
//            //////////////////////////////////////////
//            //compute PHRED scores and assign genotype
//            //////////////////////////////////////////
//    		uint32_t pl_rr = (uint32_t) round(-10*log_p_rr);
//    		uint32_t pl_ra = (uint32_t) round(-10*log_p_ra);
//    		uint32_t pl_aa = (uint32_t) round(-10*log_p_aa);
//    		
//    		uint32_t min = pl_rr;
//    		std::string bestGenotype = "0/0";
//    		if (pl_ra < min)
//    		{    
//    		    min = pl_ra;
//    		    bestGenotype = "0/1";
//    		}
//    		if (pl_aa < min)
//    		{    
//    		    min = pl_aa;
//    		    bestGenotype = "1/1";
//    		}
//    			
//    		pl_rr -= min;
//    		pl_ra -= min;
//    		pl_aa -= min;
//    		
//    		if (pl_rr+pl_ra+pl_aa==0)
//    		{
//    		    bestGenotype = "./.";
//    	    }
//    	    		
//    		std::stringstream ss;
//    		ss << chrom << "\t" 
//    		   << (pos0+1) << "\t" 
//    		   << ref  << "\t" 
//    		   << alt << "\t";
//    
//    		if (read_no)
//    		{				
//    			ss << read_no << "\t"
//    		       << rqs.str() << aqs.str() << "\t"
//    		       << mqs.str() << "\t"
//    		       << tps.str();
//    		}
//    		else
//    		{
//    			ss << "0\t.\t.\t.";
//    		}
//    		
//    		//std::cout << ss.str() << "\n";
//    	
//    		//*OUT_GLF << ss.str() << "\n";
//    	}
//        else
//        {
//            //multiallelics to be implemented
//            std::cerr << "Only Biallelic SNPs supported for this pileup, skipping this variant\n";
//            std::cerr << "no alleles: " << ivcf_rec->n_allele << "\n";          
//            for (uint32_t i=0; i<ivcf_rec->n_allele; ++i)
//            {
//                std::cerr << "\t" << ivcf_rec->d.allele[i] << "\n";
//            } 
//        }
//    }    
//    
//    void genotype_indel()
//    {
//        if (ivcf_rec->n_allele==2)
//        {
//            read_ids.clear();	
//                 
//            log_p_rr = 0;
//            log_p_ra = 0;
//            log_p_aa = 0;        
//                            
//        	rqs.str("");
//    		aqs.str("");
//    		mqs.str("");
//            tps.str("");
//            
//            const char* chrom = bcf_get_chrom(ivcf_hdr, ivcf_rec);
//            uint32_t pos = bcf_get_pos0(ivcf_rec);
//            char* ref = bcf_get_indel_ref(ivcf_rec);
//            char* alt = bcf_get_indel_alt(ivcf_rec);
//            
//            //get probe
//            int32_t plen = 0;
//            bcf_get_info_string(ivcf_hdr, ivcf_rec, "REFPROBE", &refProbe);
//            bcf_get_info_string(ivcf_hdr, ivcf_rec, "ALTPROBE", &altProbe);
//            bcf_get_info_int(ivcf_hdr, ivcf_rec, "PLEN", plen);
//            
//            vcf_format1(ivcf_hdr, ivcf_rec, &str);
//            //std::cerr << str.s << "\n";
//            //std::cerr << refProbe.s << ":" << altProbe.s << ":" << plen << "\n";
//            
//           // uint32_t plen = boost::lexical_cast<uint32_t>(plenstring.s);
//    		
//    		int32_t probeLength = refProbe.l;
//                                		 	        
//            int32_t variantLengthDifference = (int32_t)strlen(alt)-(int32_t)strlen(ref);
//            	
//            if (!(iter = bam_itr_queryi(isam_idx, bam_name2id(isam_hdr, bcf_get_chrom(ivcf_hdr, ivcf_rec)), ivcf_rec->pos, ivcf_rec->pos+1)))
//            {
//                std::cerr << "fail to parse regions\n";
//            }
//    
//            uint32_t read_no = 0;
//            while(bam_itr_next((BGZF*)isam->fp, iter, srec)>=0 && read_no<255)
//            {
//                //has secondary alignment or fail QC or is duplicate or is unmapped
//                if (bam_get_flag(srec) & 0x0704) 
//                {
//                    //std::cerr << "fail flag\n";
//                    continue;
//                }
//            
//                //make this setable
//                if (bam_get_mapq(srec)<13)
//    	        {
//    	            //std::cerr << "low mapq\n";
//                    //++readFilteredNo;
//    	            //ignore poor quality mappings
//    	            continue;
//    	        }
//                         
//    	        //std::map<std::string, int> readIDs;	
//    			if(read_ids.find(bam_get_qname(srec))==read_ids.end())
//    			{
//    			    //assign id to reads to ease checking of overlapping reads
//    				read_ids[bam_get_qname(srec)] = read_ids.size();
//    			}
//    			//ignore overlapping paired end
//    			else
//    		    {
//    		        
//    		        read_ids.erase(bam_get_qname(srec));
//    		       // std::cerr << "Mate pair\n";
//    		        continue;
//    		    }
//      
//                //read name
//                char* read_name = bam_get_qname(srec);
//                
//                //sequence
//                readseq.l=0;
//                bam_get_seq_string(srec, &readseq);
//    
//                //qual
//                bam_get_qual_string(srec, &readqual, readseq.s);
//                
//                //map qual
//                uint32_t mapQual = bam_get_mapq(srec);
//         
//         		std::string refCigar = "";
//                std::string altCigar = "";
//                    
//                double refllk = 0, altllk = 0;
//                
//    		   	lhmm_ref.align(refllk, refProbe.s, readseq.s, readqual.s);
//    	        lhmm_ref.computeLogLikelihood(refllk, lhmm_ref.getPath(), readqual.s);
//    	       
//    	        lhmm_alt.align(altllk, altProbe.s, readseq.s, readqual.s);
//                lhmm_alt.computeLogLikelihood(altllk, lhmm_alt.getPath(), readqual.s);
//                
//                std::string pad = "\t";
//    	               
//    	        if (debug)
//    	        {
//    		        std::cerr << pad << "==================\n"; 
//    		        std::cerr << pad <<  (read_no+1) << ") " << chrom << ":" << (pos+1) << ":" << ref << ":" << alt << ":" << variantLengthDifference << "\n";
//                    std::cerr << pad <<  read_name << "\n";
//                    std::cerr << pad << "==================\n";
//                    std::cerr << pad << "ref probe     " << refProbe.s << " (" << plen << "/" << probeLength << ")\n";
//                    std::cerr << pad << "read sequence " << readseq.s  << "\n";    
//                    std::cerr << pad << "==================\n";
//                    lhmm_ref.printAlignment(pad);
//                    std::cerr << pad << "ref llk: " << refllk << "\n";
//    		        std::cerr << pad << "expected indel location: " << plen+1 << "\n";
//    		        std::cerr << pad << "==================\n";
//    		        std::cerr << pad << "==================\n";
//                    std::cerr << pad << "alt probe     " << altProbe.s << " (" << plen << ")\n";
//                    std::cerr << pad << "read sequence " << readseq.s  << "\n";    
//                    std::cerr << pad << "==================\n";
//                    lhmm_alt.printAlignment(pad);
//                    std::cerr << pad << "alt llk: " << altllk << "\n";
//                    std::cerr << pad << "expected indel location: " << plen+1 << "\n";
//                    std::cerr << pad << "==================\n\n";
//    	        }
//    					           
//                //check if the Insertions and Deletions are at the expected places.
//                //deletion
//                uint32_t ref_rpos=0, alt_rpos=0;
//                if (variantLengthDifference<0)
//                {
//                    //deletion
//                    if (!deletion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                         insertion_start_exists(lhmm_alt, plen+1, alt_rpos))
//                    {
//                        if (refllk<altllk)
//                        {
//                            swap(refllk, altllk);
//                        }
//                        
//                        //ref allele
//                        if (debug)
//                            std::cerr << pad << "DEL: REF ALLELE\n"; 
//                    }
//                    else if (deletion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                            !insertion_start_exists(lhmm_alt, plen+1, alt_rpos)) 
//                    {
//                        if (refllk>altllk)
//                        {
//                            swap(refllk, altllk);
//                        }    
//                        
//                        //alt allele
//                        if (debug)
//                            std::cerr << pad << "DEL: ALT ALLELE\n"; 
//                    }
//                    else
//                    {
//                        refllk = 0;
//                        altllk = 0;
//                        if (debug)
//                            std::cerr << pad << "Deletion not at expected location, set as ambiguous\n"; 
//                    }
//                }
//                else if (variantLengthDifference>0)
//                { 
//                    //insertion
//                    if (!insertion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                         deletion_start_exists(lhmm_alt, plen+1, alt_rpos))
//                    {
//                        if (refllk<altllk)
//                        {
//                            swap(refllk, altllk);
//                        }    
//                        
//                        //ref allele
//                        if (debug)
//                            std::cerr << pad << "INS: REF ALLELE\n"; 
//                    }
//                    else if (insertion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                            !deletion_start_exists(lhmm_alt, plen+1, alt_rpos)) 
//                    {
//                        if (refllk>altllk)
//                        {
//                            swap(refllk, altllk);
//                        } 
//                        
//                        //alt allele
//                        if (debug)
//                            std::cerr << pad << "INS: ALT ALLELE\n"; 
//                    }
//                    else
//                    {
//                        refllk = 0;
//                        altllk = 0;
//                        if (debug)
//                            std::cerr << pad << "Insertion not at expected location "  << plen+1 << ", set as ambiguous\n"; 
//                    }
//                }
//         
//    	        if (debug)
//    	        {    
//    		        std::cerr << pad << "++++++++++++++++++\n";
//    	            std::cerr << pad << "reflk " << refllk << "\n";
//    			    std::cerr << pad << "altlk " << altllk << "\n";
//                }
//    
//                log_p_rr = lt.elog10prod(log_p_rr, refllk);
//                log_p_ra = lt.elog10prod(log_p_ra, lt.elog10prod(-log10(2), lt.elog10sum(refllk, altllk)));
//                log_p_aa = lt.elog10prod(log_p_aa, altllk);    
//    
//                uint32_t baseqr = 0, baseqa = 0;
//                //char allele = 'N';
//                uint32_t rpos = 0;
//                if (refllk>altllk)
//                {
//                	//allele = 'R';            
//                    rpos = alt_rpos;
//                    baseqa = round(-10*(altllk-refllk));
//                }
//                else if (refllk<altllk)
//                {
//                	//allele = 'A';            	    
//              	    rpos = ref_rpos;
//              	    baseqr = round(-10*(refllk-altllk));
//                }
//                
//                rqs << (uint8_t)(baseqr>93? 126 : baseqr+33);
//    			aqs << (uint8_t)(baseqa>93? 126 : baseqa+33);
//       			mqs << (uint8_t)(mapQual>93? 126 : mapQual+33);
//      		    int32_t tp =  pack_tp(!bam_is_rev(srec), bam_is_fread1(srec), rpos, bam_get_l_qseq(srec));
//    		    tps << std::hex << (tp<16?"0":"") << tp << std::dec ;
//    			++read_no;     
//            }
//    
//    		//////////////////////////////////////////
//            //compute PHRED scores and assign genotype
//            //////////////////////////////////////////
//    		uint32_t pl_rr = (uint32_t) round(-10*log_p_rr);
//    		uint32_t pl_ra = (uint32_t) round(-10*log_p_ra);
//    		uint32_t pl_aa = (uint32_t) round(-10*log_p_aa);
//    		
//    		uint32_t min = pl_rr;
//    		std::string bestGenotype = "0/0";
//    		if (pl_ra < min)
//    		{    
//    		    min = pl_ra;
//    		    bestGenotype = "0/1";
//    		}
//    		if (pl_aa < min)
//    		{    
//    		    min = pl_aa;
//    		    bestGenotype = "1/1";
//    		}
//    			
//    		pl_rr -= min;
//    		pl_ra -= min;
//    		pl_aa -= min;
//    		
//    		if (pl_rr+pl_ra+pl_aa==0)
//    		{
//    		    bestGenotype = "./.";
//    	    }
//       		
//    		ss.str("");
//    		ss << chrom << "\t" 
//    		   << (pos+1) << "\t" 
//    		   << ref  << "\t" 
//    		   << alt << "\t";
//    
//    		if (read_no)
//    		{				
//    			ss << read_no << "\t"
//    		       << rqs.str() << aqs.str() << "\t"
//    		       << mqs.str() << "\t"
//    		       << tps.str();
//    		}
//    		else
//    		{
//    			ss << "0\t.\t.\t.";
//    		}
//    		
//    		*OUT_GLF << ss.str() << "\n";
//    			
//            
//        }
//        else
//        {
//            //multiallelics to be implemented
//            std::cerr << "Only Biallelic Indels supported for this pileup, skipping this variant\n";
//            std::cerr << "no alleles: " << ivcf_rec->n_allele << "\n";          
//            for (uint32_t i=0; i<ivcf_rec->n_allele; ++i)
//            {
//                std::cerr << "\t" << ivcf_rec->d.allele[i] << "\n";
//            }   
//        }
//    }
    
    private:
        
    bool deletion_start_exists(LHMM& lhmm, uint32_t pos, uint32_t& rpos)
    {
        rpos = 0;
        for (uint32_t i=0; i<lhmm.indelStatusInPath.size(); ++i) 
        {
            if (lhmm.indelStatusInPath[i]=='D' &&
                lhmm.indelStartsInX[i]==pos)
            {
                rpos = lhmm.indelStartsInY[i];
                return true;
            }    
        }
        
        return false;
    }
    
    bool insertion_start_exists(LHMM& lhmm, uint32_t pos, uint32_t& rpos)
    {
        rpos = 0;
        for (uint32_t i=0; i<lhmm.indelStatusInPath.size(); ++i) 
        {
            if (lhmm.indelStatusInPath[i]=='I' &&
                lhmm.indelStartsInX[i]==pos)
            {
                rpos = lhmm.indelStartsInY[i];
                return true;
            }    
        }
        
        return false;
    }
    
    void swap(double& a, double& b)
    {
        b = (a=a+b) - b;
        a -= b; 
    }   

    int32_t pack_tp(bool is_fwd, bool is_read1, uint32_t rpos, uint32_t rlen)
    {    
        if (debug)
            std::cerr << "\tSTR:" << is_fwd << " READMATE:" << is_read1 << " RPOS:" << rpos << " RLEN:"<< rlen << "\n"; 
        
        int32_t t = 0;
        if (debug)
            std::cerr << "\t1 " << std::hex << t << "\n";
        
        t |= is_fwd ? 0 : 128;
        if (debug)
            std::cerr << "\t2 " << std::hex << t << "\n";
        
        t |= is_read1 ? 0 : 64;
        if (debug)
            std::cerr << "\t3 " << std::hex << t << "\n";
        
        int32_t e = rpos<rlen/2 ? rpos : (rlen-rpos);
        t |= (e<64 ? e : 63); 
        if (debug)
            std::cerr << "\t4 " << std::hex << t << "\n" << std::dec;
        
        return t;
    }
};
 
}

void program(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.pileup();    
    igor.print_stats();
}
