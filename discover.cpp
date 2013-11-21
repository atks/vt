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
		
VariantHunter::VariantHunter(uint32_t _evidenceAlleleCountCutoff, double _fractionalEvidenceAlleleCountCutoff, uint32_t _baseqCutoff, GenomeSequence* _myRefSeq, uint32_t _vtype, IFILE _OUT_VCF)
: bufferSize(800),
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
  baseqCutoff(_baseqCutoff),
  evidenceAlleleCountCutoff(_evidenceAlleleCountCutoff),
  fractionalEvidenceAlleleCountCutoff(_fractionalEvidenceAlleleCountCutoff),
  myRefSeq(_myRefSeq),
  vtype(_vtype),
  OUT_VCF(_OUT_VCF),
  debug(false)
{	
};

void VariantHunter::processRead(SamRecord& record)
{	
    //extract relevant information from Sam Record
    uint32_t pos = record.get1BasedPosition();
	const char* readSeq = record.getSequence();
	std::string cigarString;
	generateCigarString(cigarString, record.getCigar());
	const char* cigar = cigarString.c_str();
	const char* qual = record.getQuality();
	char a;
	bool BAQexists = false;
	const char* bq = NULL;
	
	if (record.getStringTag("BQ")!=NULL)
	{    
	    bq = record.getStringTag("BQ")->c_str();
        BAQexists = true;
    }
    
    uint32_t q;
    						
	//basically equivalent to emptying the buffer
	extractCandidateVariants(chrom, pos);
	
	std::string genomeSeqString;
	myRefSeq->getString(genomeSeqString, chromNo, (genomeIndex_t) pos-1, cigarString.size()+2);
    const char* genomeSeq = genomeSeqString.c_str();
    
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
	for (uint32_t cigarPos=0; cigarPos<cigarString.size(); ++cigarPos)
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
			    
            if (genomeSeq[genomeSeqPos]!=readSeq[readSeqPos] && q>=baseqCutoff)
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


void VariantHunter::extractCandidateVariants()
{
    extractCandidateVariants(chrom, 0, true);
};

/**
Empty buffer to recover space.
*/
void VariantHunter::extractCandidateVariants(const char* chrom, uint32_t pos, bool flush)
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
	std::string seq;
	char anchor, ref;
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
    			myRefSeq->getString(seq, chromNo, (genomeIndex_t) startGenomePos, 1);
    			anchor = seq.at(0);
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
    				    if (i->second>= evidenceAlleleCountCutoff && 
    				        ((double)i->second/(double) N[start]) >= fractionalEvidenceAlleleCountCutoff && 
    				        anchor!='N' && (i->first).find_first_of('N')==std::string::npos)
    				    {
    					    (*OUT_VCF) << chrom << "\t" << (startGenomePos) << "\t.\t";
    					    (*OUT_VCF) << anchor << "\t" << i->first << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";;
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
    				
    				myRefSeq->getString(seq, chromNo, (genomeIndex_t) (startGenomePos-1), 1);
    				anchor = seq.at(0);
    				
    				for (std::map<std::string, uint32_t>::iterator i = indel_alts.begin(); i!= indel_alts.end(); ++i)
    				{	
    				    //make sure that we do not output alleles with N bases.
    				    if (i->second>= evidenceAlleleCountCutoff && 
    				        ((double)i->second/(double) N[start]) >= fractionalEvidenceAlleleCountCutoff && 
    				        anchor!='N' && (i->first).find_first_of('N')==std::string::npos)
    				    {
    				        (*OUT_VCF) << chrom << "\t" << (startGenomePos-1) << "\t.\t";
    				        const char* deletedAllele = i->first.c_str();
    				        ++deletedAllele;
    				        (*OUT_VCF) << anchor << deletedAllele;
    					    --deletedAllele;
    					    (*OUT_VCF) << "\t" << i->first.at(0) << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";
    					}
    				}					
    			}
		    }
		    
		    if (vtype&1)
		    {
    			//handling SNPs
    			snp_alts.clear();
    			myRefSeq->getString(seq, chromNo, (genomeIndex_t) (startGenomePos), 1);
    			ref = seq.at(0);
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
    				    if (i->second>= evidenceAlleleCountCutoff && 
    				        ((double)i->second/(double) N[start]) >= fractionalEvidenceAlleleCountCutoff && 
    				        ref!='N' && (i->first)!='N')
    				    {
    					    (*OUT_VCF) << chrom << "\t" << startGenomePos << "\t.\t";
    					    (*OUT_VCF) <<ref << "\t" << i->first << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";
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
    				    myRefSeq->getString(seq, chromNo, (genomeIndex_t) (startGenomePos), i->first.size());
    				    
    				    //make sure that we do not output alleles with N bases.
    				    if (i->second>= evidenceAlleleCountCutoff && 
    				        ((double)i->second/(double) N[start]) >= fractionalEvidenceAlleleCountCutoff && 
    				        seq.find_first_of('N')==std::string::npos && (i->first).find_first_of('N')==std::string::npos)
    				    {
    					    (*OUT_VCF) << chrom << "\t" << startGenomePos << "\t.\t";
    					    (*OUT_VCF) << seq << "\t" << i->first << "\t.\t.\t.\tE:N\t" << i->second << ":" << N[start] << "\n";;
    				    }
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
		
	
int discover(int argc, char ** argv)
{	
	//options
    std::string outputVCFFileName;
   	std::string inputBAMFileName;
   	std::string	genomeFAFileName;
   	std::string sampleID;
	uint32_t mapqCutoff;
	uint32_t baseqCutoff;
	//takes on
	//snps
	//mnps
	//indels
	std::string variant_type;
	uint32_t evidenceAlleleCountCutoff;
	double fractionalEvidenceAlleleCountCutoff;
			    
	try 
	{
		std::string desc = 
"Discovers variants.\n\
$path = /net/fantasia/home/atks/programs/vt\n\
e.g. $path/vt discover -b $path/in.bam -o - -g ref.fa -v snps \n\
e.g. bam mergeBam --in a.bam --in b.bam -o - | $path/vt discover -b - -o out.sites.vcf -g ref.fa -v all\n\
e.g. samtools view in.bam | $path/vt discover -b $path/in.bam -o - -g ref.fa\n";		    
		
   		std::string version = "0.5";
		TCLAP::CmdLine cmd(desc, ' ', version);
		TCLAP::ValueArg<std::string> argOutputVCFFileName("o", "output-vcf-file", "Output VCF file", true, "-", "string", cmd);
		TCLAP::ValueArg<std::string> argInputBAMFileName("b", "input-bam-file", "Input BAM file", true, "", "string", cmd);
		TCLAP::ValueArg<std::string> argSampleID("s", "sample-id", "Sample ID", true, "", "string", cmd);
		TCLAP::ValueArg<std::string> argGenomeFAFileName("g", "genome-fa-file", "Genome FASTA file", true, "/net/fantasia/home/atks/ref/genome/human.g1k.v37.fa", "string", cmd);
		TCLAP::ValueArg<uint32_t> argMAPQCutoff("m", "mapq-cutoff", "MAPQ Cutoff, only alignments with map quality >= mapq are considered (default:20)", false, 20, "integer", cmd);
		TCLAP::ValueArg<uint32_t> argBASEQCutoff("q", "q-cutoff", "BASE Cutoff, only bases with BAQ >= baseq are considered (default:13)", false, 13, "integer", cmd);
		TCLAP::ValueArg<std::string> argVariantType("d", "variant-type", "Variant Types, takes on any combinations of the values snps,mnps,indels comma delimited (default:snps,mnps,indels)", true, "snps,mnps,indels", "string", cmd);
		TCLAP::ValueArg<uint32_t> argEvidenceAlleleCountCutoff("e", "evidence-count-cutoff", "Sample based evidence filter for candidate allele (default 2)", false, 2, "integer", cmd);
		TCLAP::ValueArg<double> argFractionalEvidenceAlleleCountCutoff("f", "evidence-fraction-cutoff", "Sample based fractional evidence filter for candidate allele (default:0.1)", false, 0.1, "float", cmd);
				
		cmd.parse(argc, argv);

		outputVCFFileName = argOutputVCFFileName.getValue();
		inputBAMFileName = argInputBAMFileName.getValue();
		sampleID = argSampleID.getValue();	
		genomeFAFileName = argGenomeFAFileName.getValue();
		mapqCutoff = argMAPQCutoff.getValue();
	    baseqCutoff = argBASEQCutoff.getValue();	
    	variant_type = argVariantType.getValue();
    	evidenceAlleleCountCutoff = argEvidenceAlleleCountCutoff.getValue();
	    fractionalEvidenceAlleleCountCutoff = argFractionalEvidenceAlleleCountCutoff.getValue();
	    
		std::clog << "discover v0.577\n\n";

		std::clog << "Options: [b] Input BAM File               " << inputBAMFileName << "\n";
		std::clog << "         [o] Output VCF File              " << outputVCFFileName << "\n";
        std::clog << "         [s] Sample ID                    " << sampleID << "\n";
        std::clog << "         [g] Reference FA File            " << genomeFAFileName << "\n";
        std::clog << "         [m] MAPQ Cutoff                  " << mapqCutoff << "\n";
        std::clog << "         [q] Base quality Cutoff          " << baseqCutoff << "\n";
            
        std::clog << "         [v] Variant Type(s)              " << variant_type << "\n";
        std::clog << "         [e] Evidence cutoff              " << evidenceAlleleCountCutoff << "\n";
        std::clog << "         [f] Fractional Evidence Cutoff   " << fractionalEvidenceAlleleCountCutoff<< "\n\n";
        
	}
	catch (TCLAP::ArgException &e) 
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		abort();
	}

	//===============
	//Genome Sequence
	//===============	
	GenomeSequence myRefSeq(genomeFAFileName);    
	//if failed (i.e. the memory mapped file does not exist)
	if (myRefSeq.open()) 
	{
		//check if creation of memory mapped file is successful
		if (myRefSeq.create(false)) 
		{
			std::cerr << "Failed creating memory map file of the genomic reference file.\n";
		}
		
		//attempt opening it again
		if (myRefSeq.open()) 
		{
			std::cerr << "Failed opening memory map file of the genomic reference file.\n";
		}
	} 

	int32_t ftype = fileType(outputVCFFileName, 'r'); 
	IFILE OUT_VCF = ifopen(outputVCFFileName.c_str(), "w", (ftype==0?InputFile::UNCOMPRESSED:(ftype==1?InputFile::GZIP:InputFile::BGZF)));
	
	(*OUT_VCF) << "##fileformat=VCFv4.1\n";
	(*OUT_VCF) << "##FORMAT=<ID=E,Number=1,Type=Integer,Description=\"Number of reads containing evidence of the alternate allele\">\n";
	(*OUT_VCF) << "##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Total number of reads at a candidate locus with reads that contain evidence of the alternate allele\">\n";
	(*OUT_VCF) << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleID << "\n";
	
	std::vector<std::string> variant_types;
    split(variant_types, ',', variant_type);
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
    
	VariantHunter variantHunter(evidenceAlleleCountCutoff, fractionalEvidenceAlleleCountCutoff, baseqCutoff, &myRefSeq, vtype, OUT_VCF);
	
 	SamFile IN_SAM;
	SamFileHeader header;
    if(!IN_SAM.OpenForRead(inputBAMFileName.c_str()))
        fprintf(stderr, "%s\n", IN_SAM.GetStatusMessage());
    if(!IN_SAM.ReadHeader(header))
        fprintf(stderr, "%s\n", IN_SAM.GetStatusMessage());
  	std::string inputBAMIndexFileName = inputBAMFileName;
  	inputBAMIndexFileName.append(".bai");
  	if (inputBAMFileName!="-" && !IN_SAM.ReadBamIndex(inputBAMIndexFileName.c_str()))
  	{    
  		std::cout << "samin status" << IN_SAM.GetStatus() <<"\n";
	}
  	IN_SAM.setSortedValidation(SamFile::COORDINATE);	    	
    uint16_t excludeFlag = 0x0704;  	                
	SamRecord record;

	uint32_t pos;
	const char* readSeq;
	uint32_t len;
	std::string cigar;
	const char* qual;
	uint32_t mapq;
	std::string genomeSeq;
	uint32_t count = 0;
	std::vector<std::string> chromosomes;
	std::string chrom = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y";
    //std::string chrom = "20"; //todo: region lists and a region option - mutex it.
    split(chromosomes, ',', chrom);
	
	if (inputBAMFileName!="-")
	{
    	for (uint32_t i=0; i<chromosomes.size(); ++i)	
    	{	
    		//std::cerr << "chromosome " << chromosomes[i] << "\n";
    		
    		if(!IN_SAM.SetReadSection(header.getReferenceID(chromosomes[i].c_str())))
    		{	
    		    std::cerr << "Cannot read section";			
    			abort();
    		}
    	
    		variantHunter.reset(chromosomes[i].c_str());
    		//uint32_t chromNo = myRefSeq.getChromosome(chromosomes[i].c_str());
            std::map<std::string, int > readIDStart;
            std::map<std::string, int > readIDEnd;
                			
    	    while (IN_SAM.ReadRecord(header, record))
    	    {
    	        uint16_t flag = record.getFlag();
    	        if(flag & excludeFlag)
    	        {
    	            //1. unmapped
    	            //2. secondary alignment
    	            //3. not passing QC
    	            //4. PCR or optical duplicate
    	            continue;
    	        }
    	        
    	        //remove the downstream overlapping read pair
                std::string readName(record.getReadName());
    			if(readIDStart.find(readName)==readIDStart.end())
    			{
    			    readIDStart[readName] = record.get1BasedPosition();
    				readIDEnd[readName] = record.get1BasedPosition()+record.getReadLength()-1;
    			}
    			else
    		    {
    		        int32_t start = record.get1BasedPosition();
    				int32_t end = record.get1BasedPosition()+record.getReadLength()-1;
    				//overlapping paired end reads
    		        if (readIDEnd[readName]>=start && readIDStart[readName]<=end)
    		        {
    		            continue;
    		        }
    		        
    		        readIDStart.erase(readName);
                    readIDEnd.erase(readName);
    		    }
    	
    			pos = record.get1BasedPosition();
    			readSeq = record.getSequence();
    			len = record.getReadLength();
    			qual = record.getQuality();
    			mapq = record.getMapQuality();
    			
    			//reads that have too many inconsistencies with the reference
    	        if (mapq<mapqCutoff)
                {
                    //includes short alignments too
                    //close by indels!
                    continue;
                }
    	
//    			if (0 && pos > 1519200 && pos < 1519450)
//    			{		
//    				std::cout << "##################" << "\n";
//    				std::cout << "count:" << count << "\n";		
//    				std::cout << chrom << ":" << pos << "\n";
//    				std::cout << "read :" << readSeq << "\n";
//    				std::cout << "len  :" <<  cigar.size() << "\n";
//    				std::cout << "cigar:" << cigar << "\n";
//    				std::cout << "qual  :" << qual << "\n";
//    				std::cout << "ref  :" << genomeSeq << "\n";
//    				std::cout << "##################" << "\n";
//    			}
    			
    			variantHunter.processRead(record);
    			++count;
    	    }
    	    
    	    //flush
    	    variantHunter.extractCandidateVariants();
    	}
    }
    else
    {
        std::string currentChrom = "-1";
        
        variantHunter.reset(currentChrom.c_str());
		//uint32_t chromNo = myRefSeq.getChromosome(chromosomes[i].c_str());
        std::map<std::string, int > readIDStart;
        std::map<std::string, int > readIDEnd;
        
        while (IN_SAM.ReadRecord(header, record))
        {
            const char* chrom = record.getReferenceName();
            
            if (chrom!=currentChrom)
            {
                currentChrom = std::string(chrom);
                variantHunter.reset(currentChrom.c_str());
                readIDStart.clear();
                readIDEnd.clear();
            }
            
            uint16_t flag = record.getFlag();
            if(flag & excludeFlag)
            {
                //1. unmapped
                //2. secondary alignment
                //3. not passing QC
                //4. PCR or optical duplicate
                continue;
            }
            
            //remove the downstream overlapping read pair
            std::string readName(record.getReadName());
        	if(readIDStart.find(readName)==readIDStart.end())
        	{
        	    readIDStart[readName] = record.get1BasedPosition();
        		readIDEnd[readName] = record.get1BasedPosition()+record.getReadLength()-1;
        	}
        	else
            {
                int32_t start = record.get1BasedPosition();
        		int32_t end = record.get1BasedPosition()+record.getReadLength()-1;
        		//overlapping paired end reads
                if (readIDEnd[readName]>=start && readIDStart[readName]<=end)
                {
                    continue;
                }
                
                readIDStart.erase(readName);
                readIDEnd.erase(readName);                
            }
        
        	pos = record.get1BasedPosition();
        	readSeq = record.getSequence();
        	len = record.getReadLength();
        	qual = record.getQuality();
        	mapq = record.getMapQuality();
        	
        	//reads that have too many inconsistencies with the reference
            if (mapq<mapqCutoff)
            {
                //includes short alignments too
                //close by indels!
                continue;
            }
        
        	if (0 && pos > 1519200 && pos < 1519450)
        	{		
        		std::cout << "##################" << "\n";
        		std::cout << "count:" << count << "\n";		
        		std::cout << chrom << ":" << pos << "\n";
        		std::cout << "read :" << readSeq << "\n";
        		std::cout << "len  :" <<  cigar.size() << "\n";
        		std::cout << "cigar:" << cigar << "\n";
        		std::cout << "qual  :" << qual << "\n";
        		std::cout << "ref  :" << genomeSeq << "\n";
        		std::cout << "##################" << "\n";        	
        	}
        	
        	variantHunter.processRead(record);
        	
        	++count;
        //	if (count>50000)
        //	{
        //	    break;
        //  }
        }
        
    }
	ifclose(OUT_VCF);

    return 0;
}
