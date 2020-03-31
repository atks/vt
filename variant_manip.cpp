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
 * Checks if the REF sequence of a VCF entry is consistent.
 *
 * returns
 * 2 - bases are not consistent with unmasked bases
 * 1 - bases are consistent against unmasked sequences
 * 0 - bases are consistent
 */
int32_t VariantManip::is_not_ref_consistent(bcf_hdr_t *h, bcf1_t *v)
{
    const char* chrom = bcf_get_chrom(h, v);
    uint32_t pos0 = bcf_get_pos0(v);
    char* vcf_ref = bcf_get_ref(v);
    uint32_t rlen = strlen(vcf_ref);

    int32_t ref_len = 0;
    char *ref = faidx_fetch_seq(fai, chrom, pos0, pos0+rlen-1, &ref_len);
    if (!ref)
    {
        fprintf(stderr, "[%s:%d %s] failure to extract base from fasta file: %s:%d-%d\n", __FILE__, __LINE__, __FUNCTION__, chrom, pos0, pos0+rlen-1);
        fprintf(stderr, "FAQ: http://genome.sph.umich.edu/wiki/Vt#1._vt_cannot_retrieve_sequences_from_my_reference_sequence_file\n");
        exit(1);
    }

    int32_t is_not_consistent = 0;
    for (uint32_t i=0; i<ref_len; ++i)
    {
        if (toupper(vcf_ref[i]) != toupper(ref[i]))
        {
            if (toupper(ref[i])=='N' || toupper(vcf_ref[i])=='N')
            {
                is_not_consistent = 1;
            }
            else
            {
                is_not_consistent = 2;
                break;
            }
        }
    }
    
    if (is_not_consistent)
    {
       fprintf(stderr, "[%s:%d %s] reference bases not consistent: %s:%d-%d  %s(REF) vs %s(FASTA)\n", __FILE__, __LINE__, __FUNCTION__, 
                                                                                chrom, pos0, pos0+rlen-1, vcf_ref, ref);       
    }    

    if (ref_len) free(ref);
   
    return is_not_consistent;
}

/**
 * Checks if a variant is normalized.
 * Ignores if entry is not a variant.
 */
bool VariantManip::is_normalized(bcf1_t *v)
{
    char** alleles = bcf_get_allele(v);
    int32_t n_allele = bcf_get_n_allele(v);

    if (n_allele==1) return true;

    char first_base;
    char last_base;
    size_t rlen, alen, len;
    bool exists_len_one_allele = false;
    bool first_base_same = true;
    bool last_base_same = true;

    if (n_allele==2)
    {
        rlen = strlen(alleles[0]);
        alen = strlen(alleles[1]);

        if (rlen==1&&alen==1)
        {
            return true;
        }
        else
        {
            //check if variant is reference.
            if (rlen==alen)
            {
                if (strcasecmp(alleles[0], alleles[1])==0)
                {
                    return true;
                }
            }

            //ref
            if (rlen==1) exists_len_one_allele = true;
            first_base = toupper(alleles[0][0]);
            last_base = toupper(alleles[0][rlen-1]);

            //alt
            if (alen==1) exists_len_one_allele = true;
            if (first_base!=toupper(alleles[1][0])) first_base_same = false;
            if (last_base!=toupper(alleles[1][alen-1])) last_base_same = false;

            if (last_base_same || (!exists_len_one_allele && first_base_same))
            {
                return false;
            }

            return true;
        }
    }
    else
    {
        bool same = true;
        for (size_t i=0; i<n_allele; ++i)
        {
            if (i)
            {
                len = strlen(alleles[i]);
                if (len==1) exists_len_one_allele = true;
                if (first_base!=toupper(alleles[i][0])) first_base_same = false;
                if (last_base!=toupper(alleles[i][len-1])) last_base_same = false;

                same = same && strcasecmp(alleles[i],alleles[0])==0;
            }
            else
            {
                len = strlen(alleles[0]);
                if (len==1) exists_len_one_allele = true;
                first_base = toupper(alleles[0][0]);
                last_base = toupper(alleles[0][len-1]);
            }
        }

        //reference entry
        if (same)
        {
            return true;
        }

        if (last_base_same || (!exists_len_one_allele && first_base_same))
        {
            return false;
        }

        return true;
    }
}

/**
 * Checks if a variant contains N bases.
 */
bool VariantManip::contains_N(bcf1_t *v)
{
    char** alleles = bcf_get_allele(v);
    int32_t n_allele = bcf_get_n_allele(v);

    for (uint32_t i=0; i<n_allele; ++i)
    {
        if (strchr(alleles[i], 'N'))
        {
            //symbolic allele
            if (i && alleles[i][0]=='<')
            {
                continue;
            }

            return true;
        }
    }

    return false;
}

/**
 * Classifies variants.
 */
int32_t VariantManip::classify_variant(bcf_hdr_t *h, bcf1_t *v, Variant& var)
{
    var.clear(); // this sets the type to VT_REF by default.

    var.h = h;
    var.v = v;

    bcf_unpack(v, BCF_UN_STR);
    var.chrom.assign(bcf_get_chrom(h, v));
    var.rid = bcf_get_rid(v);
    var.pos1 = bcf_get_pos1(v);
    var.beg1 = var.pos1;
    var.end1 = bcf_get_end1(v);

    char** allele = bcf_get_allele(v);
    int32_t n_allele = bcf_get_n_allele(v);

    uint32_t pos1 = var.pos1;
    int32_t pos0 = pos1-1;

    bool homogeneous_length = true;
    char* ref = allele[0];
    int32_t rlen = strlen(ref);

    if (strchr(ref, 'N'))
    {
        var.contains_N = true;
    }

    //if only ref allele, skip this entire for loop
    for (size_t i=1; i<n_allele; ++i)
    {
        int32_t type = VT_REF;

        //check for symbolic alternative alleles
        if (strchr(allele[i],'<'))
        {
            size_t len = strlen(allele[i]);
            if (len>=5)
            {
                //VN/d+
                if (allele[i][0]=='<' && allele[i][1]=='V' && allele[i][2]=='N' && allele[i][len-1]=='>' )
                {
                    for (size_t j=3; j<len-1; ++j)
                    {
                        if (allele[i][j]<'0' || allele[i][j]>'9')
                        {
                            type = VT_VNTR;
                        }
                    }
                }
                //VNTR
                else if (len==6 &&
                         allele[i][0]=='<' &&
                         allele[i][1]=='V' && allele[i][2]=='N' && allele[i][3]=='T' && allele[i][4]=='R' &&
                         allele[i][5]=='>' )
                {
                     type = VT_VNTR;
                }
                //STR
                else if (len==5 &&
                         allele[i][0]=='<' &&
                         allele[i][1]=='S' && allele[i][2]=='T' && allele[i][3]=='R' &&
                         allele[i][4]=='>' )
                {
                     type = VT_VNTR;
                }
                //ST/d+
                else if (allele[i][0]=='<' && allele[i][1]=='S' && allele[i][2]=='T' && allele[i][len-1]=='>' )
                {
                    type = VT_VNTR;

                    for (size_t j=3; j<len-1; ++j)
                    {
                        if ((allele[i][j]<'0' || allele[i][j]>'9') && allele[i][j]!='.')
                        {
                            type = VT_SV;
                        }
                    }
                }
            }

            if (type==VT_VNTR)
            {
                type = VT_VNTR;
                var.type |= type;
                var.alleles.push_back(Allele(type));
            }
            else
            {
                type = VT_SV;
                var.type |= type;
                std::string sv_type(allele[i]);
                var.alleles.push_back(Allele(type, sv_type));
            }
        }
        //checks for chromosomal breakpoints
        else if (strchr(allele[i],'[')||strchr(allele[i],']'))
        {
            type = VT_SV;
            var.type |= type;
            std::string sv_type("<BND>");
            var.alleles.push_back(Allele(type, sv_type));
        }
        //non variant record
        else if (allele[i][0]=='.' || strcmp(allele[i],allele[0])==0)
        {
            type = VT_REF;
        }
        //explicit sequence of bases
        else
        {
            kstring_t REF = {0,0,0};
            kstring_t ALT = {0,0,0};

            ref = allele[0];
            char* alt = allele[i];
            int32_t alen = strlen(alt);

            if (strchr(alt, 'N'))
            {
                var.contains_N = true;
            }

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

            var.type |= type;
            var.alleles.push_back(Allele(type, diff, alen, dlen, mlen, ts, tv));
            var.ts += ts;
            var.tv += tv;
            var.ins = dlen>0?1:0;
            var.del = dlen<0?1:0;
            var.max_dlen = var.max_dlen<dlen ? dlen : var.max_dlen;
            var.min_dlen = var.min_dlen>dlen ? dlen : var.min_dlen;

            if (REF.m) free(REF.s);
            if (ALT.m) free(ALT.s);
        }
    }

    if (var.type==VT_VNTR)
    {
        var.update_vntr_from_info_fields(h, v);
    }

    //additionally define MNPs by length of all alleles
    if (!(var.type&(VT_VNTR|VT_SV)) && var.type!=VT_REF)
    {
        if (homogeneous_length && rlen>1 && n_allele>1)
        {
            var.type |= VT_MNP;
        }
    }

    return var.type;
}

/**
 * Right trims or left extend a variant.
 */
void VariantManip::right_trim_or_left_extend(std::vector<std::string>& alleles, int32_t& pos1, const char* chrom, int32_t& left_extended, int32_t& right_trimmed)
{
    bool to_right_trim = true;
    bool to_left_extend = false;

    if (alleles.size()==1)
        return;

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
                fprintf(stderr, "[%s:%d %s] failure to extract base from fasta file: %s:%d\n", __FILE__, __LINE__, __FUNCTION__, chrom, pos1-1);
                fprintf(stderr, "FAQ: http://genome.sph.umich.edu/wiki/Vt#1._vt_cannot_retrieve_sequences_from_my_reference_sequence_file\n");
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
void VariantManip::left_trim(std::vector<std::string>& alleles, int32_t& pos1, int32_t& left_trimmed)
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