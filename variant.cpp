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

#include "variant.h"

/**
 * Constructor.
 */
Variant::Variant()
{
    type = VT_REF;
    motif = {0,0,0};
    mlen = 0;
    tlen = 0;
    alleles.clear();
}

/**
 * Destructor.
 */
Variant::~Variant()
{
    if (motif.m) free(motif.s);
}

/**
 * Classifies variants.
 */
int32_t Variant::classify_variant(bcf_hdr_t *h, bcf1_t *v,  Variant& variant, bool in_situ_left_trimming)
{
    bcf_unpack(v, BCF_UN_STR);
    return classify_variant(bcf_get_chrom(h, v), bcf_get_pos1(v), bcf_get_allele(v), bcf_get_n_allele(v), variant, in_situ_left_trimming);
}

/**
 * Classifies variants.
 */
int32_t Variant::classify_variant(bcf_hdr_t *h, bcf1_t *v)
{
    Variant variant;
    bcf_unpack(v, BCF_UN_STR);
    return classify_variant(bcf_get_chrom(h, v), bcf_get_pos1(v), bcf_get_allele(v), bcf_get_n_allele(v), variant);
}

/**
 * Classifies variants.
 */
int32_t Variant::classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele)
{
    Variant variant;
    return classify_variant(chrom, pos1, allele, n_allele, variant);
}

/**
 * Classifies variants.
 */
int32_t Variant::classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele, Variant& v, bool trimming)
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
                    if (REF.s[j]!=ALT.s[j])
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
 * Returns true if variant contains an allele that is potentially frame shifting.
 */
bool Variant::exists_frame_shift()
{
    for (size_t i=0; i<alleles.size(); ++i)
    {
        if (abs(alleles[i].dlen)%3!=0)
        {
            return true;
        }
    }

    return false;
}

/**
 * Prints variant information.
 */
void Variant::print()
{
    std::cerr << "type : " << vtype2string(type) << "\n";
    std::cerr << "rlen : " << rlen << "\n";
    //std::cerr << "motif: " << motif.s << "\n";
    std::cerr << "mlen : " << mlen << "\n";
    std::cerr << "tlen : " << tlen << "\n";
    for (int32_t i=0; i<alleles.size(); ++i)
    {
        std::cerr << "\tallele: " << i << "\n";
        std::cerr << "\t  type: " << vtype2string(alleles[i].type) << "\n";
        std::cerr << "\t  diff: " << alleles[i].diff << "\n";
        std::cerr << "\t  alen: " << alleles[i].alen << "\n";
        std::cerr << "\t  dlen: " << alleles[i].dlen << "\n";
        std::cerr << "\t  tlen: " << alleles[i].tlen << "\n";
    }
};

/**
 * Converts VTYPE to string.
 */
std::string Variant::vtype2string(int32_t VTYPE)
{
    std::string s;

    if (!VTYPE)
    {
        s += (s.size()==0) ? "" : "/";
        s += "REF";
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
        s += "INDEL";
    }

    if (VTYPE & VT_CLUMPED)
    {
        s += (s.size()==0) ? "" : "/";
        s += "CLUMPED";
    }

    return s;
}

/**
 * Clears variant information.
 */
void Variant::clear()
{
    type = VT_REF;
    motif.l = 0;
    mlen = 0;
    tlen = 0;
    alleles.clear();
};