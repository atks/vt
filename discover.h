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

#ifndef DISCOVER_VARIANT_H
#define DISCOVER_VARIANT_H

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>

#include "utils.h"

/**
*Class for mining candidate variants
*/
class VariantHunter
{			
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
	uint32_t baseqCutoff;
	uint32_t evidenceAlleleCountCutoff;
	double fractionalEvidenceAlleleCountCutoff;
	GenomeSequence* myRefSeq;
	uint32_t vtype;
	IFILE OUT_VCF;
	bool debug;
	
	/**
	 *Processes buffer to pick up variants
	 */
	void extractCandidateVariants(const char* chrom, uint32_t pos, bool flush=false);
		
	public:
	      
	/**
	 *Constructor
	 *baseqCutoff - q value cutoff to select candidate SNPs
	 *todo: implement Igor
	 */
	VariantHunter(uint32_t _evidenceAlleleCountCutoff, double _fractionalEvidenceAlleleCountCutoff, uint32_t _baseqCutoff, GenomeSequence* _myRefSeq, uint32_t _vtype, IFILE _OUT_VCF);
	
	/**
	 *Reads a read and moves it to the buffer
	 */
	void processRead(SamRecord& record);

	/**
	 *Processes buffer to pick up variants
	 */	
	void extractCandidateVariants();
	
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
	 *Increments buffer index i by j.
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
	 *Decrements buffer index i by j.
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
	 *Decrements buffer index i by 1.
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
	 *Returns the difference between 2 buffer positions
	 */
	uint32_t diff(uint32_t i, uint32_t j)
	{
		return (i>=j ? i-j : bufferSize-(j-i));
	};
	
	/**
	 *Gets the position in the buffer that corresponds to 
	 *the genome position indicated by pos.
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
	 *Reset buffer for new chromosome
	 */
	void reset(const char* _chrom)
	{
		extractCandidateVariants(_chrom, UINT_MAX);
		chrom = _chrom;
		chromNo = myRefSeq->getChromosome(_chrom);
	};
	
	/**
	 *Print buffer contents for debugging purpose
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
	
	
    //expands cigar string
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
                
                cigar.append(count, c);
                token.str("");
            }
            
            ++i;
        } 
    };
};

int discover(int argc, char ** argv);

#endif