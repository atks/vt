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

/**
 * Constructor.
 *
 * @ref_fasta_file reference sequence FASTA file.
 */
VariantManip::VariantManip(std::string ref_fasta_file)
{
    if (ref_fasta_file!="")
    {
        fai = fai_load(ref_fasta_file.c_str());
        if (fai==NULL)
        {
            fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
            exit(1);
        }
        reference_present = (fai!=NULL);
    }
};

/**
 * Constructor.
 */
VariantManip::VariantManip()
{
    reference_present = false;
}

/**
 * Converts VTYPE to string.
 */
std::string VariantManip::vtype2string(int32_t VTYPE)
{
    std::string s;

    if (!VTYPE)
    {
        s += (s.size()==0) ? "" : "/";
        s += "Reference";
    }

    if (VTYPE & VT_SNP)
    {
        s += (s.size()==0) ? "" : "/";
        s += "SNP";
    }

    if (VTYPE & VT_MNP)
    {
        s += (s.size()==0) ? "" : "/";
        s += "MNP";
    }

    if (VTYPE & VT_INDEL)
    {
        s += (s.size()==0) ? "" : "/";
        s += "Indel";
    }

    if (VTYPE & VT_CLUMPED)
    {
        s += (s.size()==0) ? "" : "/";
        s += "Clumped";
    }

    if (VTYPE & VT_SV)
    {
        s += (s.size()==0) ? "" : "/";
        s += "Structural Variants";
    }

    return s;
}

/**
 * Converts VTYPE to string.
 */
void VariantManip::vtype2string(int32_t vtype, kstring_t *s)
{
    s->l = 0;

    if (vtype & VT_SNP)
    {
        if (s->l) kputc(',', s);
        kputs("SNP", s);
    }

    if (vtype & VT_MNP)
    {
        if (s->l) kputc(',', s);
        kputs("MNP", s);
    }

    if (vtype & VT_INDEL)
    {
        if (s->l) kputc(',', s);
        kputs("INDEL", s);
    }

    if (vtype & VT_CLUMPED)
    {
        if (s->l) kputc(',', s);
        kputs("CLUMPED", s);
    }
}

/**
 * Detects near by STRs.
 */
bool VariantManip::detect_str(bcf_hdr_t *h, bcf1_t *v, Variant& variant)
{
    return detect_str(bcf_get_chrom(h, v), bcf_get_pos1(v), variant);
}

/**
// * Detects near by STRs.
 */
bool VariantManip::detect_str(const char* chrom, uint32_t pos1, Variant& variant)
{
    int32_t ref_len;
    //STR related
    char* ru = 0;
    ru = faidx_fetch_uc_seq(fai, chrom, pos1, pos1, &ref_len);

    int32_t tract_len = 1;
    int32_t motif_len = 1;

    std::string motif = "";
    int32_t tlen = 0;

    while (1)
    {
        char* next_ru = 0;
        next_ru = faidx_fetch_uc_seq(fai, chrom, pos1+tract_len*motif_len, pos1+(tract_len)*motif_len, &ref_len);

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
            ru = faidx_fetch_uc_seq(fai, chrom, pos1, pos1+motif_len-1, &ref_len);
        }

        free(next_ru);
    }

    free(ru);

    return true;
}

/**
 * Checks if a variant is normalized.
 */
bool is_normalized(char** alleles, int32_t n_allele)
{
    char first_base;
    char last_base;
    size_t len;
    bool exists_len_one_allele = false;
    bool first_base_same = true;
    bool last_base_same = true;

    if (n_allele==2)
    {
        //ref
        len = strlen(alleles[0]);
        if (len==1) exists_len_one_allele = true;
        first_base = alleles[0][0];
        last_base = alleles[0][len-1];

        //alt
        len = strlen(alleles[1]);
        if (len==1) exists_len_one_allele = true;
        if (first_base!=alleles[1][0]) first_base_same = false;
        if (last_base!=alleles[1][len-1]) last_base_same = false;

        if (last_base_same || (!exists_len_one_allele && first_base_same))
        {
            return false;
        }

        return true;
    }
    else
    {
        for (size_t i=0; i<n_allele; ++i)
        {
            if (i)
            {
                len = strlen(alleles[i]);
                if (len==1) exists_len_one_allele = true;
                if (first_base!=alleles[i][0]) first_base_same = false;
                if (last_base!=alleles[i][len-1]) last_base_same = false;
            }
            else
            {
                len = strlen(alleles[0]);
                if (len==1) exists_len_one_allele = true;
                first_base = alleles[0][0];
                last_base = alleles[0][len-1];
            }
        }

        if (last_base_same || (!exists_len_one_allele && first_base_same))
        {
            return false;
        }

        return true;
    }
}

/**
 * Classifies variants.
 */
int32_t VariantManip::classify_variant(bcf_hdr_t *h, bcf1_t *v,  Variant& variant, bool in_situ_left_trimming)
{
    bcf_unpack(v, BCF_UN_STR);
    return classify_variant(bcf_get_chrom(h, v), bcf_get_pos1(v), bcf_get_allele(v), bcf_get_n_allele(v), variant, in_situ_left_trimming);
}

/**
 * Classifies variants.
 */
int32_t VariantManip::classify_variant(bcf_hdr_t *h, bcf1_t *v)
{
    Variant variant;
    bcf_unpack(v, BCF_UN_STR);
    return classify_variant(bcf_get_chrom(h, v), bcf_get_pos1(v), bcf_get_allele(v), bcf_get_n_allele(v), variant);
}

/**
 * Classifies variants.
 */
int32_t VariantManip::classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele)
{
    Variant variant;
    return classify_variant(chrom, pos1, allele, n_allele, variant);
}

/**
 * Classifies variants.
 */
int32_t VariantManip::classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele, Variant& v, bool trimming)
{
    int32_t pos0 = pos1-1;
    v.clear(); // this sets the type to VT_REF by default.

    bool homogeneous_length = true;

    char* ref = allele[0];
    int32_t rlen = strlen(ref);

    //if only ref allele, skip this entire for loop
    for (size_t i=1; i<n_allele; ++i)
    {
        int32_t type = VT_REF;

        //check for tags
        if (allele[i][0]=='<')
        {
            type = VT_SV; //support for SVs, might extend to STRs in future which is not really an SV

            v.type |= type;
            std::string sv_type(allele[i]);
            v.alleles.push_back(Allele(type, 0, 0, 0, 0, 0, 0, sv_type));
        }
        else
        {
            kstring_t REF = {0,0,0};
            kstring_t ALT = {0,0,0};

            ref = allele[0];
            char* alt = allele[i];
            int32_t alen = strlen(alt);

            if (rlen!=alen)
            {
                homogeneous_length = false;
            }

            //trimming
            //this is required in particular for the
            //characterization of multiallelics and
            //in general, any unnormalized variant
            int32_t rl = rlen;
            int32_t al = alen;
            if (trimming)
            {
                //trim right
                while (rl!=1 && al!=1)
                {
                    if (ref[rl-1]==alt[al-1])
                    {
                        --rl;
                        --al;
                    }
                    else
                    {
                        break;
                    }
                }

                //trim left
                while (rl !=1 && al!=1)
                {
                    if (ref[0]==alt[0])
                    {
                        ++ref;
                        ++alt;
                        --rl;
                        --al;
                    }
                    else
                    {
                        break;
                    }
                }

                kputsn(ref, rl, &REF);
                kputsn(alt, al, &ALT);

                ref = REF.s;
                alt = ALT.s;
            }


            int32_t mlen = std::min(rl, al);
            int32_t dlen = al-rl;
            int32_t diff = 0;
            int32_t ts = 0;
            int32_t tv = 0;

            if (mlen==1 && dlen)
            {
                char ls, le, ss;

                if (rl>al)
                {
                     ls = ref[0];
                     le = ref[rl-1];
                     ss = alt[0];
                }
                else
                {
                     ls = alt[0];
                     le = alt[al-1];
                     ss = ref[0];
                }

                if (ls!=ss && le!=ss)
                {
                    ++diff;

                    if ((ls=='G' && ss=='A') ||
                        (ls=='A' && ss=='G') ||
                        (ls=='C' && ss=='T') ||
                        (ls=='T' && ss=='C'))
                    {
                        ++ts;
                    }
                    else
                    {
                        ++tv;
                    }
                }
            }
            else
            {
                for (int32_t j=0; j<mlen; ++j)
                {
                    if (ref[j]!=alt[j])
                    {
                        ++diff;

                        if ((ref[j]=='G' && alt[j]=='A') ||
                            (ref[j]=='A' && alt[j]=='G') ||
                            (ref[j]=='C' && alt[j]=='T') ||
                            (ref[j]=='T' && alt[j]=='C'))
                        {
                            ++ts;
                        }
                        else
                        {
                            ++tv;
                        }
                    }
                }
            }

            //substitution variants
            if (mlen==diff)
            {
                type |= mlen==1 ? VT_SNP : VT_MNP;
            }

            //indel variants
            if (dlen)
            {
                type |= VT_INDEL;
            }

            //clumped SNPs and MNPs
            if (diff && diff < mlen) //internal gaps
            {
                type |= VT_CLUMPED;
            }

            v.type |= type;
            v.alleles.push_back(Allele(type, diff, alen, dlen, 0, mlen, ts, tv));

            if (REF.m) free(REF.s);
            if (ALT.m) free(ALT.s);
        }
    }

    //additionally define MNPs by length of all alleles
    if (!(v.type&VT_SV))
    {
        if (homogeneous_length && rlen>1)
        {
            v.type |= VT_MNP;
        }
    }

    return v.type;
}

/**
 * Right trims or left extend a variant.
 */
void VariantManip::right_trim_or_left_extend(std::vector<std::string>& alleles, uint32_t& pos1, const char* chrom, uint32_t& left_extended, uint32_t& right_trimmed)
{
    bool to_right_trim = true;
    bool to_left_extend = false;

    while (to_right_trim || to_left_extend)
    {
        //checks if right trimmable or left extendable
        to_right_trim = true;
        to_left_extend = false;
        for (size_t i=0; i<alleles.size(); ++i)
        {
            if (!alleles[i].empty())
            {
                if (alleles[0].at(alleles[0].size()-1) != alleles[i].at(alleles[i].size()-1))
                {
                    to_right_trim = false;
                    //do not break here!!! you need to check for empty alleles that might exist!!!
                }
                
                if (pos1==1 && alleles[i].size()==1)
                {
                    to_right_trim = false;
                    break;
                }
            }
            else
            {
                to_right_trim = false;
                to_left_extend = true;
                break;
            }
        }
        
        if (to_right_trim)
        {
            for (size_t i=0; i<alleles.size(); ++i)
            {
                alleles[i].erase(alleles[i].size()-1);
            }

            ++right_trimmed;
        }

        if (to_left_extend)
        {
            --pos1;
            int ref_len = 0;

            char *ref = faidx_fetch_uc_seq(fai, chrom, pos1-1, pos1-1, &ref_len);
            if (!ref)
            {
                fprintf(stderr, "[%s:%d %s] failure to extrac base from fasta file: %s:%d: >\n", __FILE__, __LINE__, __FUNCTION__, chrom, pos1-1);
                exit(1);
            }
            char base = ref[0];
            free(ref);

            for (size_t i=0; i<alleles.size(); ++i)
            {
                alleles[i].insert(0, 1, base);
            }

            ++left_extended;
        }
    }
};

/**
 * Left trims a variant.
 */
void VariantManip::left_trim(std::vector<std::string>& alleles, uint32_t& pos1, uint32_t& left_trimmed)
{
    bool to_left_trim =  true;

    while (to_left_trim)
    {
        //checks if left trimmable.
        for (size_t i=0; i<alleles.size(); ++i)
        {
            if (alleles[i].size()==1 || alleles[i].at(0)!=alleles[0].at(0))
            {
                to_left_trim = false;
                break;
            }
        }

        if (to_left_trim)
        {
            for (size_t i=0; i<alleles.size(); ++i)
            {
                alleles[i].erase(0, 1);
            }

            ++pos1;
            ++left_trimmed;
        }
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
            base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), pos1-1, pos1-1, &ref_len);
            preamble.append(1,base[0]);
            bases[base[0]] = 1;
            if (ref_len>0) free(base);
            ++i;
        }

        bases.clear();
        i=0;
        uint32_t alleleLength = alleles[0].size();
        while (bases.size()<4 || postamble.size()<min_flank_length)
        {
            base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), pos1+alleleLength+i, pos1+alleleLength+i, &ref_len);
            postamble.append(1,base[0]);
            bases[base[0]] = 1;
            if (ref_len>0) free(base);
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
        int32_t ref_len = 0;
        while (bases.size()<4 && preamble.size()<min_flank_length)
        {
            base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), pos1-i-1, pos1-i-1, &ref_len);
            preamble.append(1,base[0]);
            bases[base[0]] = 1;
            ++i;
            if (base[0]=='N')
            {
                break;
            }
            if (ref_len>0) free(base);
        }

        preambleLength = preamble.size();

        for (size_t i=0; i<alleles.size(); ++i)
        {
            probes[i] = std::string(preamble.rbegin(), preamble.rend()).append(probes[i]);
        }
    }
}

/**
 * Iteratively called function for generating a haplotype.
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
    if (currentDiff<alleles.size() || length<=2*gald+flankLength)
    {
        std::map<std::string, uint32_t> probeHash;
        //extend probes
        for (uint32_t i=0; i<alleles.size(); ++i)
        {
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
                    char* base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), start1 , start1, &ref_len);
                    probes[i].append(1, base[0]);
                    if (ref_len>0) free(base);
                }
            }
            probeHash[probes[i]] = 1;
        }

        currentDiff = probeHash.size();
        ++length;
        //std::cerr << probes[0] << "\n" << probes[1] << "\n";
        generate_probes(chrom, pos1, flankLength, currentDiff, length, gald, diff, alleles, probes);
    }
}