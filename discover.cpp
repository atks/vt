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

#include "discover.h"

typedef struct 
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

/**
 * Class for mining candidate variants.
 *
 * Processes an align read and records the variants observed against the reference
 */
class VariantHunter
{		
    public:
	/**
	 * Constructor
	 * baseq_cutoff - q value cutoff to select candidate SNPs
	 */
    VariantHunter(uint32_t vtype, 
                  uint32_t evidence_allele_count_cutoff, 
                  double fractional_evidence_allele_count_cutoff, 
                  uint32_t baseq_cutoff, 
                  faidx_t *fai, 
                  BCFOrderedWriter *odw)
    : vtype(vtype),
      evidence_allele_count_cutoff(evidence_allele_count_cutoff),
      fractional_evidence_allele_count_cutoff(fractional_evidence_allele_count_cutoff),
      baseq_cutoff(baseq_cutoff),
      fai(fai),
      odw(odw),
      bufferSize(800),
      X(bufferSize),
      Y(bufferSize),
      I(bufferSize), 
      D(bufferSize),
      N(bufferSize,0),
      REF(bufferSize),
      start(0), end(0),
      emptyBufferSpace(0),
      minEmptyBufferSize(400),
      startGenomePos(0),
      maxUsedBufferSizeThreshold(bufferSize-minEmptyBufferSize),
      maxINDELLength(50),
      debug(false)
    {	
    };
	
	/**
	 * Reads a read and moves it to the buffer
	 */
	void processRead(bam1_t *b)
    {	
        //extract relevant information from Sam Record
        uint32_t pos = bam_get_pos1(b);
    	bam_get_seq_string(b, read_seq);
    	const char* readSeq = read_seq->s;
    	bam_get_cigar_string(b, cigar_string);
    	const char* cigar = cigar_string->s;
    	bam_get_qual_string(b, qual_string, readSeq);
    	const char* qual = qual_string->s;
    	char a;
    	bool BAQexists = false;
    	const char* bq = NULL;
    	int32_t ref_len;
    	
//    	if (record.getStringTag("BQ")!=NULL)
//    	{    
//    	    bq = record.getStringTag("BQ")->c_str();
//          BAQexists = true;
//      }
        
        uint32_t q;
        						
    	//basically equivalent to emptying the buffer
    	extractCandidateVariants(chrom, pos);
    	
    	std::string genomeSeqString;
    	const char* genomeSeq = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos-1, pos+cigar_string->l-1, &ref_len);
        			
    	//myRefSeq->getString(genomeSeqString, chromNo, (genomeIndex_t) pos-1, cigarString.size()+2);
        //const char* genomeSeq = genomeSeqString.c_str();
        
    	uint32_t genomeSeqPos = 1;
    	uint32_t readSeqPos = 0;
    	uint32_t curPos = getCurPos(pos); //current buffer index
    	bool lastPositionHadSnp = false;
    	bool mnpAlleleConstructionInProgress = false;
    	bool insInit = true;
    	bool delInit = true;
    	uint32_t lastSnpPos = 0;
    	uint32_t mnpInitPos = 0;
    	char mnpInitBase = 'N';
    	uint32_t insInitPos = 0;
    	uint32_t delInitPos = 0;
    	
    //	if (((pos>239300 && pos<239500) && debug))
    //	{
    //		std::cout << "===============\n";	
    //		std::cout << "ADD READ\n";	
    //		std::cout << "startGenomePos      : " << startGenomePos << "\n";
    //		std::cout << "endGenomePos        : " << startGenomePos + diff(end,start)  << "\n";	
    //		std::cout << "start               : " << start << "\n";
    //		std::cout << "end                 : " << end << "\n";
    //		std::cout << "curPos              : " << curPos << "\n";
    //	}
    	
    //			if(pos<239411 && pos+strlen(readSeq)>239411)
    //			{
    //			    std::cerr << "=====================\n";
    //			    std::cerr << "CIGAR " <<  cigar << "\n";
    //			    std::cerr << "POS " <<  pos << "\n";    
    //		        std::cerr << "SEQ " <<  readSeq << "\n";    
    //	        
    //			        
    //			}
    
    //    if (0 && pos>239400 && pos<239421)
    //    {
    //        std::cerr << "=====================\n";
    //        std::cerr << "POS   " <<  pos << "\n";
    //        std::cerr << "CIGAR " <<  cigar << "\n";
    //        std::cerr << "CIGAR STRING " <<  record.getCigar() << "\n";
    //        std::cerr << "SEQ   " <<  readSeq << "\n";    
    //    }
    	
    	if (isEmpty())
    		startGenomePos = pos;
    	
    	//cycle through the cigar string, this cigar string has essentially been expanded
    	for (uint32_t cigarPos=0; cigarPos<cigar_string->l; ++cigarPos)
    	{
    		a = cigar[cigarPos];
    						
    		if (a=='M')
    		{
    			REF[curPos] = genomeSeq[genomeSeqPos];
    			//get quality
    			q = (qual[readSeqPos]-33) - (BAQexists?(bq[readSeqPos]-64):0);
                
    //          if ((bq[readSeqPos]-64)>(qual[readSeqPos]-33))
    //          {   
    //              std::cerr << "===========\n"; 
    //    			std::cerr << "baq " << q << "\n";
    //    			std::cerr << "q "<< qual[readSeqPos]-33 << "\n";
    //    			std::cerr << "bq "<< bq[readSeqPos]-64 << "\n";
    //			}
    			    
                if (genomeSeq[genomeSeqPos]!=readSeq[readSeqPos] && q>=baseq_cutoff)
                {
                    X[curPos].push_back(readSeq[readSeqPos]);
                    
                    //initialize mnp
                    if (lastPositionHadSnp && !mnpAlleleConstructionInProgress)
                    {
                        mnpAlleleConstructionInProgress  = true;
                        mnpInitPos = lastSnpPos;
                        lastPositionHadSnp = false;
                        
                        Y[mnpInitPos].push_back("");
    			        Y[mnpInitPos].back().append(1, mnpInitBase);
    				} 
                                                                   
                    if (mnpAlleleConstructionInProgress)
                    {
                       Y[mnpInitPos].back().append(1, readSeq[readSeqPos]);
                    }
                    
                    lastPositionHadSnp = true;
                    lastSnpPos = curPos;
                    mnpInitBase = readSeq[readSeqPos];
                }
                else
                {
                    lastPositionHadSnp = false;
                    mnpAlleleConstructionInProgress = false;
                }
                
    			insInit = true;
    			delInit = true;
    
    			++N[curPos];					
    			add(curPos);
    			++genomeSeqPos;
    			++readSeqPos;
    		}
    		else if (a=='I')
    		{
    			if (insInit)
    			{   
    			    insInitPos = minus(curPos,1);
    			    I[insInitPos].push_back("");
    			    //place anchor
    			    I[insInitPos].back().append(1, (readSeqPos!=0?readSeq[readSeqPos-1]:genomeSeq[genomeSeqPos-1]));
    								
    				//if (readSeqPos==0)
    				//    ++N[curPos];
    			}
    			
    			lastPositionHadSnp =false;
                mnpAlleleConstructionInProgress = false;
    			I[insInitPos].back().append(1, readSeq[readSeqPos]);
    			insInit = false;
    			delInit = true;
    			++readSeqPos;					
    		}
    		else if (a=='D')
    		{				
    			REF[curPos] = genomeSeq[genomeSeqPos];	
    			if (delInit)
    			{
    			    D[curPos].push_back("");
    			    D[curPos].back().append(1, (readSeqPos!=0?readSeq[readSeqPos-1]:genomeSeq[genomeSeqPos-1]));
    			    ++N[curPos];
    				delInitPos = curPos;
    				
    				//if (readSeqPos==0)
    				//    ++N[curPos];
    			}
    			D[delInitPos].back().append(1, genomeSeq[genomeSeqPos]);
    			
    			lastPositionHadSnp =false;
                mnpAlleleConstructionInProgress = false;
    			delInit = false;
    			insInit = true;
    			add(curPos);
    			++genomeSeqPos;
    		}
    		else //S, H and others
    		{
    			++readSeqPos;
    		}
    	}
    	
    //	if (debug && (pos>239300 && pos<239500))
    //	{
    //		std::cout << "final curPos        : " << curPos << "\n";
    //		std::cout << "final startGenomePos: " << startGenomePos << "\n";
    //		std::cout << "final endGenomePos  : " << startGenomePos + diff(end,start)  << "\n";	
    //		std::cout << "final start         : " << start << "\n";
    //		std::cout << "final end           : " << end << "\n";
    //		std::cout << "====================\n";	
    //		printBuffer();
    //	}
    			
    };

	/**
	 *Processes buffer to pick up variants
	 */	
	void extractCandidateVariants()
	{
        extractCandidateVariants(chrom, 0, true);
    };

	/**
	 * Reset buffer for new chromosome
	 */
	void reset(const char* chrom)
	{
		extractCandidateVariants(chrom, UINT_MAX);
		this->chrom = chrom;
	};
    	
	private:

	uint32_t bufferSize;
	std::vector<std::vector<char> > X; // contains read bases that differ from the genome
	std::vector<std::vector<std::string> > Y; // contains multiple consecutive read bases that differ from the genome
	std::vector<std::vector<std::string> > I; //contains inserted bases
	std::vector<std::vector<std::string> > D; //contains reference bases that are deleted
	std::vector<uint32_t> N; // number of evidences observed here - combination of X, I and D
	std::vector<char> REF;
	std::vector<std::string> ALT;
	const char* chrom;
	uint32_t chromNo;
	uint32_t start, end;
	uint32_t emptyBufferSpace;
	uint32_t minEmptyBufferSize;
	uint32_t startGenomePos;
	uint32_t maxUsedBufferSizeThreshold;
	uint32_t maxINDELLength;
	uint32_t baseq_cutoff;
	uint32_t evidence_allele_count_cutoff;
	double fractional_evidence_allele_count_cutoff;
	faidx_t *fai;
	uint32_t vtype;
	kstring_t *s;
	kstring_t *alleles;
	kstring_t *read_seq;
	kstring_t *qual_string;
	kstring_t *cigar_string;
	
	bcf1_t *v;
	BCFOrderedWriter *odw;
	bool debug;
	
	/**
	 * Processes buffer to pick up variants
	 * Empty buffer to recover space.
     */
    void extractCandidateVariants(const char* chrom, uint32_t pos, bool flush=false)
    {
        //variable to tell when to stop flushing
        uint32_t stop = 0;
    	if (flush)
    	{
    	    stop = end;
    	}
    	else if (isEmpty())
    	{
    		return;				
    	}
    	//extract when separated
    	else if (startGenomePos+(diff(end,start))<pos)
    	{
    		if (pos==160141 && debug)
    		{
    			std::cout << "***********************\n";
    			std::cout << "flush buffer segregated\n";
    				
    		}
    		//flush buffer completely
    		stop = end;
    	}
    	//extract when overlapping
    	else 
    	{
    	    if (pos-startGenomePos>maxUsedBufferSizeThreshold)
    	    {
    	        stop = add(start, pos-startGenomePos);
    
    			if (pos==160141 && debug)
    			{
    				std::cout << "************************\n";
    				std::cout << "flush buffer overlapping\n";
    			}
    
    	    }
    	  	else
    		{	
    			return;
    		}
    	    
    //			    stop = add(start, pos)
    //			    				
    //				if (diff(end,start)>maxUsedBufferSizeThreshold)
    //				{
    //					stop = add(start, (pos>maxINDELLength+startGenomePos ? (pos-maxINDELLength)-startGenomePos : 0));
    //				}
    	}
    				
    //	if (pos==160141 && debug)
    //	{
    //		std::cout << "pos               : " << pos << "\n";
    //		std::cout << "stop              : " << stop << "\n";
    //		std::cout << "usedBufferSize    : " << diff(end,start) << "\n";
    //		std::cout << "bufferSize        : " << bufferSize << "\n";
    //		std::cout << "minEmptyBufferSize: " << minEmptyBufferSize << "\n";
    //		std::cout << "startGenomePos    : " << startGenomePos << "\n";
    //		std::cout << "endGenomePos      : " << startGenomePos + diff(end,start)  << "\n";
    //		std::cout << "start             : " << start << "\n";
    //		std::cout << "end               : " << end << "\n";
    //			
    //		std::cout << pos << "," << diff(end,start) << "," <<  bufferSize-minEmptyBufferSize << "\n";
    //	}
    	
    	std::map<char, uint32_t> snp_alts;
        std::map<std::string, uint32_t> mnp_alts;
    	std::map<std::string, uint32_t> indel_alts;
    	char anchor, ref;
    	int32_t ref_len;
    	//print out candidate variants
    	while (start!=stop)
    	{	    	
    	    //assayed position
    		if (N[start]>=1)
    		{
    		    if (vtype&4)
    	        {
        			//handling insertions
        			indel_alts.clear();        			
                    char* seq = faidx_fetch_seq(fai, chrom, startGenomePos-1, startGenomePos-1, &ref_len);
                    anchor = seq[0];
        			if (ref_len<1) free(seq);
        			
        			if (I[start].size()!=0)
        			{   
        				for (uint32_t i=0; i<I[start].size(); ++i)
        				{
        					if (indel_alts.find(I[start][i])==indel_alts.end())
        					{
        						indel_alts[I[start][i]] = 1;
        					}
        					else
        					{
        						++indel_alts[I[start][i]];
        					}
        				}
        				
        			    for (std::map<std::string, uint32_t>::iterator i =indel_alts.begin(); i!=indel_alts.end(); ++i)
        				{
        				    //make sure that we do not output alleles with N bases.
        				    if (i->second>= evidence_allele_count_cutoff && 
        				        ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff && 
        				        anchor!='N' && (i->first).find_first_of('N')==std::string::npos)
        				    {
        				        bcf_set_chrom(odw->hdr, v, chrom);
        				        bcf_set_pos1(v, startGenomePos);
        				        alleles->l = 0;
        			            kputc(anchor, alleles);
                                kputc(',', alleles);
                                kputs(i->first.c_str(), alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles->s);
        				        bcf_update_info_int32(odw->hdr, v, "E", &i->second, 1);
        				        bcf_update_info_int32(odw->hdr, v, "N", &N[start], 1);
        				        //(*OUT_VCF) << anchor << "\t" << i->first << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";;
        				    }
        				}
        			}
        			
        			//handling deletions
        			indel_alts.clear();
        			if (D[start].size()!=0)
        			{
        				for (uint32_t i=0; i<D[start].size(); ++i)
        				{
        					if (indel_alts.find(D[start][i])==indel_alts.end())
        					{
        						indel_alts[D[start][i]] = 1;
        					}
        					else
        					{
        						++indel_alts[D[start][i]];
        					}
        				}
        				
       				    char* seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), startGenomePos-1, startGenomePos-1, &ref_len);
                        anchor = seq[0];
            			if (ref_len<1) free(seq);
        				//myRefSeq->getString(seq, chromNo, (genomeIndex_t) (startGenomePos-1), 1);
        				//anchor = seq.at(0);
        				
        				for (std::map<std::string, uint32_t>::iterator i = indel_alts.begin(); i!= indel_alts.end(); ++i)
        				{	
        				    //make sure that we do not output alleles with N bases.
        				    if (i->second>= evidence_allele_count_cutoff && 
        				        ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff && 
        				        anchor!='N' && (i->first).find_first_of('N')==std::string::npos)
        				    {
        				        bcf_set_chrom(odw->hdr, v, chrom);
        				        bcf_set_pos1(v, startGenomePos);
                                alleles->l = 0;
        			            kputc(anchor, alleles);
                                kputc(',', alleles);
                                const char* deletedAllele = i->first.c_str();
        				        ++deletedAllele;
        				        kputs(deletedAllele, alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles->s);
        				        bcf_update_info_int32(odw->hdr, v, "E", &i->second, 1);
        				        bcf_update_info_int32(odw->hdr, v, "N", &N[start], 1);
        				        
//        				        (*OUT_VCF) << chrom << "\t" << (startGenomePos-1) << "\t.\t";
//        				        const char* deletedAllele = i->first.c_str();
//        				        ++deletedAllele;
//        				        (*OUT_VCF) << anchor << deletedAllele;
//        					    --deletedAllele;
//        					    (*OUT_VCF) << "\t" << i->first.at(0) << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";
        					}
        				}					
        			}
    		    }
    		    
    		    if (vtype&1)
    		    {
        			//handling SNPs
        			snp_alts.clear();
        			
        			//myRefSeq->getString(seq, chromNo, (genomeIndex_t) (startGenomePos), 1);
        			//ref = seq.at(0);
        			char* seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), startGenomePos-1, startGenomePos-1, &ref_len);
                    ref = seq[0];
            		if (ref_len<1) free(seq);
        				
        			
        			if (X[start].size()!=0)
        			{
        				for (uint32_t i=0; i<X[start].size(); ++i)
        				{
        					if (snp_alts.find(X[start][i])==snp_alts.end())
        					{
        						snp_alts[X[start][i]] = 1;
        					}
        					else
        					{
        						++snp_alts[X[start][i]];
        					}
        				}
        				
        			    for (std::map<char, uint32_t>::iterator i =snp_alts.begin(); i!=snp_alts.end(); ++i)
        				{
        				    //make sure that we do not output alleles with N bases.
        				    if (i->second>= evidence_allele_count_cutoff && 
        				        ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff && 
        				        ref!='N' && (i->first)!='N')
        				    {
        				        bcf_set_chrom(odw->hdr, v, chrom);
        				        bcf_set_pos1(v, startGenomePos);
                                alleles->l = 0;
        			            kputc(anchor, alleles);
                                kputc(',', alleles);
                                kputc(i->first, alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles->s);
        				        bcf_update_info_int32(odw->hdr, v, "E", &i->second, 1);
        				        bcf_update_info_int32(odw->hdr, v, "N", &N[start], 1);
        				        
//        					    (*OUT_VCF) << chrom << "\t" << startGenomePos << "\t.\t";
//        					    (*OUT_VCF) <<ref << "\t" << i->first << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";
        				    }
        				}
        			}
    		    }
    			
    			if (vtype&2)
    		    {
        			//handling MNPs
        			mnp_alts.clear();
        			if (Y[start].size()!=0)
        			{
        				for (uint32_t i=0; i<Y[start].size(); ++i)
        				{
        					if (mnp_alts.find(Y[start][i])==mnp_alts.end())
        					{
        						mnp_alts[Y[start][i]] = 1;
        					}
        					else
        					{
        						++mnp_alts[Y[start][i]];
        					}
        				}
        				
        			    for (std::map<std::string, uint32_t>::iterator i =mnp_alts.begin(); i!=mnp_alts.end(); ++i)
        				{
        				    char* seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), startGenomePos-1, startGenomePos+i->first.size()-1, &ref_len);
//        				    myRefSeq->getString(seq, chromNo, (genomeIndex_t) (startGenomePos), i->first.size());
        				    
        				    //make sure that we do not output alleles with N bases.
        				    if (i->second>= evidence_allele_count_cutoff && 
        				        ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff && 
        				        strchr(seq, 'N') && (i->first).find_first_of('N')==std::string::npos)
        				    {
        				                				        bcf_set_chrom(odw->hdr, v, chrom);
        				        bcf_set_pos1(v, startGenomePos);
                                alleles->l = 0;
        			            kputc(anchor, alleles);
                                kputc(',', alleles);
                                kputs(i->first.c_str(), alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles->s);
        				        bcf_update_info_int32(odw->hdr, v, "E", &i->second, 1);
        				        bcf_update_info_int32(odw->hdr, v, "N", &N[start], 1);
        				        
//        					    (*OUT_VCF) << chrom << "\t" << startGenomePos << "\t.\t";
//        					    (*OUT_VCF) << seq << "\t" << i->first << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";;
        				    }
        				    
        				    if (ref_len<1) free(seq);
        				}
        			}
        		}
    		
    		}				
    		
    		Y[start].clear();
    		X[start].clear();
    		I[start].clear();
    		D[start].clear();
    		N[start] = 0;
    						
    		add(start);
    		++startGenomePos;	
    	}
    	
    //			//clean up final position too
    //			if (isEmpty())
    //			{
    //			    X[start].clear();
    //				I[start].clear();
    //				D[start].clear();
    //				N[start] = 0;
    //			}
    	
    //	if (pos==160141  && debug)
    //	{
    //		std::cout << "final start : " << start << "\n";
    //		std::cout << "final end   : " << end << "\n";
    //		std::cout << "*************\n";
    //		printBuffer();
    //	}
    };
		
	

	/**
	 *Checks if buffer is empty
	 */		
	bool isEmpty()
	{
		return start==end;			
	};
	
	/**
	 *Increments buffer index i by 1.
	 */
	void add(uint32_t& i)
	{
	    if (i>=bufferSize)
	    {
	        std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
	        exit(1);
	    }
	    
		uint32_t temp = (i+1)%bufferSize;
		i = end==i ? (end=temp) : temp;	
	};
	
	/**
	 * Increments buffer index i by j.
	 */
	uint32_t add(uint32_t i, uint32_t j)
	{
	    if (i>=bufferSize)
	    {
	        std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
	        exit(1);
	    }
	    
		return (i+j)%bufferSize;
	};

	/**
	 * Decrements buffer index i by j.
	 */
	uint32_t minus(uint32_t& i, uint32_t j)
	{
	    if (i>=bufferSize)
	    {
	        std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
	        exit(1);
	    }
	    
		return (i>=j ? i-j : bufferSize-(j-i));
	};
	
	/**
	 * Decrements buffer index i by 1.
	 */
	void minus(uint32_t& i)
	{
	    if (i>=bufferSize)
	    {
	        std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
	        exit(1);
	    }
	    
		i = (i>=1 ? i-1 : bufferSize-1);
	};
	
	/**
	 * Returns the difference between 2 buffer positions
	 */
	uint32_t diff(uint32_t i, uint32_t j)
	{
		return (i>=j ? i-j : bufferSize-(j-i));
	};
	
	/**
	 * Gets the position in the buffer that corresponds to 
	 * the genome position indicated by pos.
	 */
	uint32_t getCurPos(uint32_t genomePos)
	{
	    //when buffer is empty
	    if (isEmpty())
		{
			startGenomePos = genomePos;
			return start;
		}
		else
		{
		    if (genomePos-startGenomePos>bufferSize)
		    {
		        std::cerr << "overflow buffer\n" ;
		        //should allow for unbuffering here
		        
		    }
			return (start + (genomePos-startGenomePos))%bufferSize;
		}
	};
	
	/**
	 * Print buffer contents for debugging purpose
	 */
	void printBuffer()
	{
		std::cout << "PRINT BUFFER" << "\n";
		std::cout << "usedBufferSize: " << diff(end,start) << "\n";
		uint32_t curPos = start;
		uint32_t genomePos = startGenomePos;
			
		while (curPos!=end)
		{
			std::cout << genomePos << "\t" << curPos << "\t" << REF[curPos] << "\t";
				
			for (uint32_t j=0; j<I[curPos].size(); ++j)
			{
				std::cout << I[curPos][j] << ","; 
			}
			
			for (uint32_t j=0; j<D[curPos].size(); ++j)
			{
				std::cout << D[curPos][j] << ","; 
			}
			
			std::cout << "\t" <<  N[curPos] << "\n";
			
			add(curPos);
			++genomePos;
		}
	};
		
    /**
     * expands cigar string
     */
    void generateCigarString(std::string& cigar, const char* cigarString)
    {
        cigar.clear();
        int32_t i=0, lastIndex = strlen(cigarString)-1;
        std::stringstream token;
        
        if (lastIndex<0)
        {
            return;
        }
        char c;
        bool seenM = false;
        
        while (i<=lastIndex)
        {
        	c = cigarString[i];
        	
        	//captures the count
            if (c<'A') 
            {
            	token << c;
            }
    
            if (c>'A' ||
                i==lastIndex) 
            {
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
                
                int32_t count;
                std::string s = token.str();
                str2int32(s, count);
                cigar.append(count, c);
                token.str("");
            }
            
            ++i;
        } 
    };
};

namespace
{

class Igor : Program
{
    public:

	///////////
    //options//
    ///////////
    std::vector<GenomeInterval> intervals;
    std::string output_vcf_file;
   	std::string input_bam_file;
   	std::string	ref_fasta_file;
   	std::string sample_id;
	uint32_t mapq_cutoff;
	uint32_t baseq_cutoff;
	//takes on snps, mnps, indels
	std::string variant_type;
	uint32_t evidence_allele_count_cutoff;
	double fractional_evidence_allele_count_cutoff;
	
	uint16_t excludeFlag;  	
	
    ///////
    //i/o//
    ///////
    BAMOrderedReader *odr;
    BCFOrderedWriter *odw;
    bam1_t *s;
    bcf1_t *v;

    /////////
    //stats//
    /////////

    /////////
    //tools//
    /////////
    VariantHunter *variantHunter;
    
    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Discovers variants from reads in a BAM file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
         	TCLAP::ValueArg<std::string> arg_input_bam_file("b", "input-bam-file", "Input BAM file", true, "", "string", cmd);
    		TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            
            TCLAP::ValueArg<std::string> arg_sample_id("s", "sample-id", "Sample ID", true, "", "string", cmd);
    		TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_mapq_cutoff("m", "mapq-cutoff", "MAPQ Cutoff, only alignments with map quality >= mapq are considered (default:20)", false, 20, "integer", cmd);
    		TCLAP::ValueArg<uint32_t> arg_baseq_cutoff("q", "q-cutoff", "BASE Cutoff, only bases with BAQ >= baseq are considered (default:13)", false, 13, "integer", cmd);
    		TCLAP::ValueArg<std::string> arg_variant_type("d", "variant-type", "Variant Types, takes on any combinations of the values snps,mnps,indels comma delimited (default:snps,mnps,indels)", true, "snps,mnps,indels", "string", cmd);
    		TCLAP::ValueArg<uint32_t> arg_evidence_allele_count_cutoff("e", "evidence-count-cutoff", "Sample based evidence filter for candidate allele (default 2)", false, 2, "integer", cmd);
    		TCLAP::ValueArg<double> arg_fractional_evidence_allele_count_cutoff("f", "evidence-fraction-cutoff", "Sample based fractional evidence filter for candidate allele (default:0.1)", false, 0.1, "float", cmd);
    		
            cmd.parse(argc, argv);
            input_bam_file = arg_input_bam_file.getValue();
    		parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            
            sample_id = arg_sample_id.getValue();	
    		ref_fasta_file = arg_ref_fasta_file.getValue();
            mapq_cutoff = arg_mapq_cutoff.getValue();
    	    baseq_cutoff = arg_baseq_cutoff.getValue();	
        	variant_type = arg_variant_type.getValue();
        	evidence_allele_count_cutoff = arg_evidence_allele_count_cutoff.getValue();
    	    fractional_evidence_allele_count_cutoff = arg_fractional_evidence_allele_count_cutoff.getValue();
    	           
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
    	
        excludeFlag = 0x0704;  	                

        odr = new BAMOrderedReader(input_bam_file, intervals);
        s = bam_init1();
   
        odw = new BCFOrderedWriter(output_vcf_file, 1);
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=E,Number=1,Type=Integer,Description=\"Number of reads containing evidence of the alternate allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Total number of reads at a candidate locus with reads that contain evidence of the alternate allele\">");
        bcf_hdr_add_sample(odw->hdr, sample_id.c_str());

    	std::vector<std::string> variant_types;
        split(variant_types, ",", variant_type);
        uint32_t vtype = 0;
        for (uint32_t i = 0; i<variant_types.size(); ++i)
        {
            if (variant_types[i] == "snps")
            {
                vtype |= 1;
            }
            else if (variant_types[i] == "mnps")
            {
                vtype |= 2;
            }
            else if (variant_types[i] == "indels")
            {
                vtype |= 4;
            }
            else if (variant_types[i] == "all")
            {
                vtype = 7;
            }
        }
        
        
        
        
    	
        ////////////////////////
        //stats initialization//
        ////////////////////////
   
        ////////////////////////
        //tools initialization//
        ////////////////////////
        faidx_t *fai = fai_load(ref_fasta_file.c_str());
        
        variantHunter = new VariantHunter(vtype,
    	                            evidence_allele_count_cutoff, 
    	                            fractional_evidence_allele_count_cutoff, 
    	                            baseq_cutoff, 
    	                            fai, 
    	                            odw);
    	
    }

    void discover()
    {
    	uint32_t pos;
    	const char* readSeq;
    	uint32_t len;
    	std::string cigar;
    	const char* qual;
    	uint32_t mapq;
    	std::string genomeSeq;
    	uint32_t count = 0;
    	std::vector<std::string> chromosomes;
    	

        const char* current_chrom = "";
        
        variantHunter->reset(current_chrom);
		
   	    khash_t(rdict) *reads = kh_init(rdict);
        khiter_t k;
        int32_t ret;
		
		std::map<std::string, int > readIDStart;
        std::map<std::string, int > readIDEnd;
        
        while (odr->read(s))
        {
            const char* chrom = bam_get_chrom(odr->hdr, s);
            
            if (!strcmp(chrom, current_chrom))
            {
                current_chrom = chrom;
                variantHunter->reset(current_chrom);
                readIDStart.clear();
                readIDEnd.clear();
            }
            
            uint16_t flag = bam_get_flag(s);
            if(flag & excludeFlag)
            {
                //1. unmapped
                //2. secondary alignment
                //3. not passing QC
                //4. PCR or optical duplicate
                continue;
            }
            
            //remove the downstream overlapping read pair
            
//             khash_t(sdict) *m = kh_init(sdict);
//        khiter_t k;
//        int32_t ret;
//  
//        while (sr->read_next_position(current_recs))
//        {
//            kh_clear(sdict, m);
//			
//			for (uint32_t i=0; i<current_recs.size(); ++i)
//			{
//			    bcf_get_variant(current_recs[i].h, current_recs[i].v, &var);
//		        if (kh_get(sdict, m, var.s)==kh_end(m))
//		        {
//		            odw->write(current_recs[i].v);
//		            kh_put(sdict, m, strdup(var.s), &ret);
//		            ++no_unique_variants;
//		        }
//		        
//		        ++no_total_variants;
//			}
//        }
//
//        odw->close();    
            
            
            char *qname = bam_get_qname(s);
            if((k = kh_get(rdict, reads, qname))==kh_end(reads))
        	{
        	    k = kh_put(rdict, reads, strdup(qname), &ret);
        	    kh_val(reads, k) = {0,0};
                kh_val(reads, k).start1 = bam_get_pos1(s);
                kh_val(reads, k).end1 = bam_get_pos1(s)+bam_get_l_qseq(s);
                
        		//readIDEnd[readName] = record.get1BasedPosition()+record.getReadLength()-1;
        	}
        	else
            {
                int32_t start = bam_get_pos1(s);
        		int32_t end = bam_get_pos1(s)+bam_get_l_qseq(s);
        		if (kh_val(reads, k).end1>=start && kh_val(reads, k).start1<=end)
                {
                    continue;
                }
                
                kh_del(rdict, reads, k);               
            }
        
        	
        	
        	//reads that have too many inconsistencies with the reference
            if (mapq<mapq_cutoff)
            {
                //includes short alignments too
                //close by indels!
                continue;
            }
        
        	if (0)
        	{   
        	    pos = bam_get_pos1(s);
            	kstring_t seq = {0,0,0};
            	bam_get_seq_string(s, &seq);
            	len = bam_get_l_qseq(s);
            	kstring_t qual = {0,0,0};
            	bam_get_qual_string(s, &qual, seq.s);
            	mapq = bam_get_mapq(s);
        	    
        		std::cout << "##################" << "\n";
        		std::cout << "count:" << count << "\n";		
        		std::cout << chrom << ":" << pos << "\n";
        		std::cout << "read :" << seq.s << "\n";
        		//std::cout << "len  :" <<  cigar.size() << "\n";
        		//std::cout << "cigar:" << cigar << "\n";
        		std::cout << "qual  :" << qual.s << "\n";
        		//std::cout << "ref  :" << genomeSeq << "\n";
        		std::cout << "##################" << "\n";        	
        	
        	    if (seq.m) free(seq.s);
      	        if (qual.m) free(qual.s);
        	}
        	
        	variantHunter->processRead(s);
        	
        	++count;
        }
    };

    void print_options()
    {
        std::clog << "discover v" << version << "\n\n";

		std::clog << "options: [b] input BAM File               " << input_bam_file << "\n";
		std::clog << "         [o] output VCF File              " << output_vcf_file << "\n";
        std::clog << "         [s] sample ID                    " << sample_id << "\n";
        std::clog << "         [r] reference FASTA File         " << ref_fasta_file << "\n";
        std::clog << "         [m] MAPQ Cutoff                  " << mapq_cutoff << "\n";
        std::clog << "         [q] base quality Cutoff          " << baseq_cutoff << "\n";
        std::clog << "         [v] variant Type(s)              " << variant_type << "\n";
        std::clog << "         [e] evidence cutoff              " << evidence_allele_count_cutoff << "\n";
        std::clog << "         [f] fractional Evidence Cutoff   " << fractional_evidence_allele_count_cutoff<< "\n\n";    
        if (intervals.size()!=0)
        {
            std::clog << "         [i] intervals             " << intervals.size() <<  " intervals\n";
        }
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\nstats: biallelic\n";
//        std::clog << "          no. left trimmed                      : " << no_lt << "\n";
//        std::clog << "          no. left trimmed and left aligned     : " << no_lt_la << "\n";
//        std::clog << "          no. left trimmed and right trimmed    : " << no_lt_rt << "\n";
//        std::clog << "          no. left aligned                      : " << no_la << "\n";
//        std::clog << "          no. right trimmed                     : " << no_rt << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void discover(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.discover();
    igor.print_stats();
};
