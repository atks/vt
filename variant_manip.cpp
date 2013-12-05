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
    if (ref_fasta_file=="")
    {
        fai = fai_load(ref_fasta_file.c_str());
    }
};

/**
 * Converts VTYPE to string.
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
 * Classifies variants.
 */
void VariantManip:: void classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele, Variant& v)
{
    int32_t pos0 = pos1-1;
    v.clear();
    
    int32_t len_ref = strlen(allele[0]);

    for (uint32_t i=1; i<n_allele; ++i)
    {
        int32_t len_alt = strlen(allele[i]);

        uint32_t len = 0;
        //SNP or MNP or block substitution
        if (len_ref==len_alt)
        {
            if (len_ref==1)
            {
                if (allele[0][0]!=allele[i][0])
                {
                    v.type |= VT_SNP;
                }
            }
            else
            {
                bool is_mnp = true;
                for (uint32_t j=0; j<len; ++j)
                {
                    //clumped variant treated as complex variant
                    if (allele[0][j]==allele[i][j])
                    {
                        is_mnp = false;
                        v.type |= VT_COMPLEX;
                    }
                }

                v.type |= (is_mnp ? VT_MNP : 0);
            }
        }
        else
        {
            int32_t ref_len;
            
            if (len_ref==1 || len_alt==1)
            {
                v.type |= VT_INDEL;

                if (len_ref==1) v.type |= VT_INSERTION;
                if (len_alt==1) v.type |= VT_DELETION;

                //SNP adjacent to Indel
                if (allele[0][0]!=allele[i][0])
                {
                    v.type |= VT_SNP;
                }

                char* ru = 0;
                ru = faidx_fetch_seq(fai, chrom, pos0+1, pos0+1, &ref_len);
                //std::cerr << "first ru: "<< ru << "\n";

                int32_t tract_len = 1;
                int32_t motif_len = 1;

                while (1)
                {
                    char* next_ru = 0;
                    next_ru = faidx_fetch_seq(fai, chrom, pos0+tract_len*motif_len+1, pos0+(tract_len+1)*motif_len, &ref_len);
                    
                    //motif repeated
                    if (strcmp(ru, next_ru)==0)
                    {
                        //extend tract length
                        ++tract_len;
                    }
                    else //try longer motif
                    {
                        if (tract_len>1)
                        {
                            motif = std::string(ru);
                            tlen = tract_len;
                            v.type |= (VT_STR | VT_EXACT_STR);
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
                        ru = faidx_fetch_seq(fai, chrom, pos0+1, pos0+motif_len, &ref_len);
                    }

                    free(next_ru);
                }

                free(ru);
            }
            else
            {
                v.type |= VT_COMPLEX;
            }
        }
    }
}

/**
 * Left trims a variant with unnecesary nucleotides.
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

/**
 * Generates a probing haplotype with flanks around the variant of interest.
 * Flanks are equal length
 */
void VariantManip::generate_probes(const char* chrom,
                        int32_t pos1, uint32_t probeDiff, //
                        std::vector<std::string>& alleles, //store alleles
                        std::vector<std::string>& probes, //store probes
                        uint32_t min_flank_length,
                        int32_t& preambleLength) //store preamble length
{
    //map to usable number
    probes.resize(alleles.size(), "");

    //check allele lengths
    std::map<uint32_t, uint32_t> alleleLengths;
    for (uint32_t i=0; i<alleles.size(); ++i)
    {
       alleleLengths[alleles[i].size()]=1;
    }

    //for SNPs and MNPs and block substitutions
    if (alleleLengths.size()==1)
    {
        //just get flanking sequences
        //append preamble
        std::map<char, uint32_t> bases;
        std::string preamble;
        std::string postamble;
        char *base;
        uint32_t i = 1;
        int32_t ref_len;
        while (bases.size()<4 || preamble.size()<min_flank_length)
        {
            base = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos1-1, pos1-1, &ref_len);
            preamble.append(1,base[0]);
            bases[base[0]] = 1;
            if (ref_len<1) free(base);
            ++i;
        }

        bases.clear();
        i=0;
        uint32_t alleleLength = alleles[0].size();
        while (bases.size()<4 || postamble.size()<min_flank_length)
        {
            base = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos1+alleleLength+i, pos1+alleleLength+i, &ref_len);
            postamble.append(1,base[0]);
            bases[base[0]] = 1;
            if (ref_len<1) free(base);
            ++i;
        }

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            probes[i] = std::string(preamble.rbegin(), preamble.rend()).append(alleles[i]);
            probes[i] = probes[i].append(postamble);
        }

        preambleLength = preamble.size();
    }
    //for Indels and Complex Substitutions
    else
    {
        //find gald
        uint32_t min_len = alleles[0].size();
        uint32_t max_len = alleles[0].size();
        for (uint32_t i=1; i<alleles.size(); ++i)
        {
            if (alleles[i].size()<min_len)
                min_len = alleles[i].size();
            if (alleles[i].size()>max_len)
                max_len = alleles[i].size();
        }
        uint32_t gald = max_len-min_len;

        uint32_t currentDiff = 0;
        //current length of probe
        uint32_t length = 0;
        //number of point differences for each probe wrt the reference
        std::vector<uint32_t> diff(alleles.size(), 0);
        probes.resize(alleles.size(), "");

        generate_probes(chrom, pos1, min_flank_length, currentDiff, length, gald, diff, alleles, probes);

        //append preamble
        std::map<char, uint32_t> bases;
        std::string preamble;
        char* base;
        uint32_t i = 1;
        int32_t ref_len;
        while (bases.size()<4 || preamble.size()<min_flank_length)
        {
            base = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos1+i-1, pos1+i-1, &ref_len);
            //myRefSeq->getString(base, chromNo, pos-i, 1);
            preamble.append(1,base[0]);
            bases[base[0]] = 1;
            ++i;
            if (base[0]=='N')
            {
                break;
            }
            if (ref_len<1) free(base);
        }

        preambleLength = preamble.size();

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            probes[i] = std::string(preamble.rbegin(), preamble.rend()).append(probes[i]);
        }
    }
}

/**
*Iteratively called function for generating a haplotype.
*/
void VariantManip::generate_probes(const char* chrom,
                        int32_t pos1,
                        uint32_t flankLength,
                        uint32_t& currentDiff,
                        uint32_t& length,
                        uint32_t gald,
                        std::vector<uint32_t>& diff,
                        std::vector<std::string>& alleles,
                        std::vector<std::string>& probes)
{
    //std::cerr << "B:" << currentDiff << " " << alleles.size()-1 << "\n";

    if (currentDiff<alleles.size() || length<=2*gald+flankLength)
    {
        std::map<std::string, uint32_t> probeHash;
        //extend probes
        for (uint32_t i=0; i<alleles.size(); ++i)
        {
//                if (i==0)
//                {
//                    std::string base;
//                    myRefSeq->getString(base, chromNo, (genomeIndex_t) (pos+length), 1);
//                    probes[0].append(1, base.at(0));
//                }
//                else
            {
                //copy from allele
                if (length<alleles[i].size())
                {
                    probes[i].append(1,alleles[i].at(length));
                }
                else//copy from reference
                {
                    int32_t start1 = (pos1+length-alleles[i].size()+alleles[0].size()-1);
                    int32_t ref_len;
                    char* base = faidx_fetch_seq(fai, const_cast<char*>(chrom), start1 , start1, &ref_len);
                    //myRefSeq->getString(base, chromNo, (genomeIndex_t) (pos+length-alleles[i].size()+alleles[0].size()), 1);
                    probes[i].append(1, base[0]);
                    if (ref_len<1) free(base);
                }
            }
            //std::cerr << probes[i] << "\n" ;
            probeHash[probes[i]] = 1;
        }

        currentDiff = probeHash.size();
        ++length;
        //std::cerr << probes[0] << "\n" << probes[1] << "\n";
        generate_probes(chrom, pos1, flankLength, currentDiff, length, gald, diff, alleles, probes);
    }
}