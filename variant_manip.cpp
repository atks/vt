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

#include "variant_manip.h"

VariantManip::VariantManip(std::string ref_fasta_file)
{
    fai = fai_load(ref_fasta_file.c_str());
};

/**
 * Converts VTYPE to string
 */
std::string VariantManip::vtype2string(int32_t VTYPE)
{
	std::string s;

	if (!VTYPE)
	{
		s += (s.size()==0) ? "" : ";";
		s += "REF";
	}
	
	if (VTYPE & VT_SNP)
	{
		s += (s.size()==0) ? "" : ";";
		s += "SNP";
	}
	
	if (VTYPE & VT_MNP)
	{
		s += (s.size()==0) ? "" : ";";
		s += "MNP";
	}	
	
	if (VTYPE & VT_INDEL)
	{
		s += (s.size()==0) ? "" : ";";
		s += "INDEL";
	}	
	
	if (VTYPE & VT_STR)
	{
		s += (s.size()==0) ? "" : ";";
		s += "STR";
	}
	
	if (VTYPE & VT_EXACT_STR)
	{
		s += (s.size()==0) ? "" : ";";
		s += "EXACT_STR";
	}	
	
	if (VTYPE & VT_INEXACT_STR)
	{
		s += (s.size()==0) ? "" : ";";
		s += "INEXACT_STR";
	}	
	
	if (VTYPE & VT_COMPLEX)
	{
		s += (s.size()==0) ? "" : ";";
		s += "COMPLEX";
	}
	
	if (VTYPE & VT_SV)
	{
		s += (s.size()==0) ? "" : ";";
		s += "SV";
	}
	
	if (VTYPE & VT_CR)
	{
		s += (s.size()==0) ? "" : ";";
		s += "CR";
	}
	
	return s;
}

/**
 * Classifies Indels into the following categories:
 * 1. Homopolymers
 * 2. Dinucleotides
 * 3. Trinucleotides
 */
int32_t VariantManip::classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele, std::string& motif, uint32_t& tlen)
{
    int32_t pos0 = pos1-1;
	int32_t VTYPE = 0;
	
    if (n_allele==2)
    {
    	uint32_t len_ref = strlen(allele[0]);
    	uint32_t len_alt = strlen(allele[1]);
    	
        uint32_t len = 0;
        if (len_ref==len_alt)
        {
            if (len_ref==1)
            {
            	if (allele[0][0]!=allele[1][0])
            	{
                	VTYPE |= VT_SNP;
            	}
            }
            else
            {
            	for (uint32_t i=0; i<len; ++i)
            	{
            		if (allele[0][i]==allele[1][i])
            		{
            			VTYPE |= VT_COMPLEX;
            		}	
            	}
            	
                VTYPE |= VT_MNP;
            }
        }
        else
        {
            int32_t ref_len;
            //char *ref = 0;
            //ref = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos0, pos0+50, &ref_len);
			//std::cerr << "REF: " << ref << "\n";
			//if (ref) free(ref);
			
			if (len_ref==1 || len_alt==1)
			{
				VTYPE |= VT_INDEL;
				
				if (len_ref==1) VTYPE |= VT_INSERTION;
				if (len_alt==1) VTYPE |= VT_DELETION;
					
				if (allele[0][0]!=allele[1][0])
				{
					VTYPE |= VT_SNP;
				}
				
				char* ru = 0;
				ru = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos0+1, pos0+1, &ref_len);
				//std::cerr << "first ru: "<< ru << "\n";
				
				int32_t tract_len = 1;
	            int32_t motif_len = 1;
	            		
	            while (1)
	            {
	            	char* next_ru = 0;
	            	next_ru = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos0+tract_len*motif_len+1, pos0+(tract_len+1)*motif_len, &ref_len);
					//std::cerr << "\tnext_ru: "<< next_ru << "\n";
	            	
	            	//motif repeated
	                if (strcmp(ru, next_ru)==0)
	                {
	                	//extend tract length
	                    ++tract_len;
	                    //std::cerr << "\t\ttract_len: "<< tract_len << "\n";
	                }
	                else //try longer motif 
	                {
	                	if (tract_len>1)
	                    {
	                    	motif = std::string(ru);
	                        tlen = tract_len;
	                        VTYPE |= (VT_STR | VT_EXACT_STR);
	                    	free(next_ru);	
	                    	break;
	                    }
	                
	                	//not STR
	                	if (motif_len>10)
	                	{
	                		free(next_ru);	
	                		break;
	                	}	
	                
	                	free(ru);    
	                    ++motif_len;
	                    tract_len=1;
	                    ru = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos0+1, pos0+motif_len, &ref_len);
	                    //std::cerr << "new ru: "<< ru << " " << motif_len << " " << tract_len << "\n";
	                }
	                
	                free(next_ru);
	                
	            }
	            
	            free(ru);
			}
			else
			{
				VTYPE |= VT_COMPLEX;
			}
        }
    }
    else
    {
    }

    return VTYPE;
}

/**
 * Left trims a variant of unnecesary nucleotides.
 */
void VariantManip::left_trim(std::vector<std::string>& alleles, uint32_t& pos1, uint32_t& left_trimmed)
{
    bool may_left_trim =  true;

    for (uint32_t i=0; i<alleles.size(); ++i)
    {
        if (alleles[i].size()==1 || alleles[i].at(0)!=alleles[0].at(0))
        {
            may_left_trim = false;
            break;
        }
    }

    if (may_left_trim)
    {
        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            alleles[i].erase(0, 1);
        }

        ++pos1;
        ++left_trimmed;
        return left_trim(alleles, pos1, left_trimmed);
    }
};

/**
 * Left aligns a variant.
 * @todo: rewrite input as char**
 */
void VariantManip::left_align(std::vector<std::string>& alleles, uint32_t& pos1, const char* chrom, uint32_t& left_aligned, uint32_t& right_trimmed)
{
    bool may_right_trim =  true;
    bool may_left_align = false;
    char lastBase = ' ';

    for (uint32_t i=0; i<alleles.size(); ++i)
    {
        if (!alleles[i].empty())
        {
            lastBase = (lastBase != ' ') ? lastBase : alleles[i].at(alleles[i].size()-1);
            if (lastBase != alleles[i].at(alleles[i].size()-1))
            {
                may_right_trim = false;
            }
        }
        else
        {
            may_left_align = true;
            may_right_trim = false;
            break;
        }
    }

    if(may_left_align)
    {
        std::string base = "";

        --pos1;

        int ref_len = 0;

        char *ref = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos1-1, pos1-1, &ref_len);
        base = std::string(ref);
        free(ref);


        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            alleles[i].insert(0, 1, base.at(0));
        }

        ++left_aligned;
        return left_align(alleles, pos1, chrom, left_aligned, right_trimmed);
    }
    else if (may_right_trim)
    {
        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            alleles[i].erase(alleles[i].size()-1);
        }

        ++right_trimmed;
        return left_align(alleles, pos1, chrom, left_aligned, right_trimmed);
    }
};
