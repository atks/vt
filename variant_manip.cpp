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
 * Classifies Indels into the following categories:
 * 1. Homopolymers
 * 2. Dinucleotides
 * 3. Trinucleotides
 */
int32_t VariantManip::classify_variant(char* chrom, uint32_t pos1, char** allele, int32_t n_allele, std::string& motif, uint32_t& tlen)
{
    int32_t pos0 = pos1-1;

    if (n_allele==2)
    {
        uint32_t len = 0;
        if ((len=strlen(allele[0]))==strlen(allele[1]))
        {
            if (len==1)
            {
                return VT_SNP;
            }
            else
            {
                return VT_MNP;
            }
        }
        else
        {
            int32_t tract_len = 0;
            int32_t motif_len = 1;
            std::string base;

            int32_t ref_len;
            char *ref = faidx_fetch_seq(fai, chrom, pos0, pos0+1, &ref_len);

            std::string ru = base;
            uint32_t i = 2;

            while (1)
            {
                if (base==ru)
                {
                    ++tract_len;

                    i += motif_len;
                }
                else
                {
                    if (motif_len*tract_len>10)
                    {
                        motif = ru;
                        tlen = tract_len;
                        return VT_STR;
                    }
                    else if (motif_len>10)
                    {
                        return VT_INDEL;
                    }

                    tract_len = 1;
                    ++motif_len;
                    
                    int ref_len = 0;
                    char *ref = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos0, pos0+motif_len-1, &ref_len);
                    ru = std::string(ref);

                    free(ref);

                    std::cout << "\t" << ru << " i " << i << "\n";
                    i = motif_len+1;
                }
            }
        }
    }
    else
    {
    }

    return 1;
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
