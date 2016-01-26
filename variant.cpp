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
Variant::Variant(bcf_hdr_t* h, bcf1_t* v)
{
    this->v = v;

    type = classify(h, v);

    chrom = bcf_get_chrom(h, v);
    rid = bcf_get_rid(v);
    pos1 = bcf_get_pos1(v);

    no_overlapping_snps = 0;
    no_overlapping_indels = 0;
    no_overlapping_vntrs = 0;

    //attempts to update relevant information on variants
    if (type==VT_SNP)
    {
        beg1 = bcf_get_pos1(v);
        end1 = bcf_get_pos1(v);
    }
    else if (type==VT_INDEL)
    {
        int32_t *flanks = NULL;
        int32_t n = 0;
        if (bcf_get_info_int32(h, v, "FLANKS", &flanks, &n)>0)
        {
            vntr.exact_rbeg1 = flanks[0]+1;
            vntr.exact_rend1 = flanks[1]-1;
            free(flanks);
        }
        else
        {
            vntr.exact_rbeg1 = bcf_get_pos1(v) - 1;
            vntr.exact_rend1 = bcf_get_end_pos1(v) + 1;
        }

        beg1 = vntr.exact_rbeg1-1;
        end1 = vntr.exact_rend1+1;

//        int32_t *fuzzy_flanks = NULL;
//        n = 0;
//        if (bcf_get_info_int32(h, v, "FZ_FLANKS", &fuzzy_flanks, &n)>0)
//        {
//            vntr.fuzzy_rbeg1 = fuzzy_flanks[0];
//            vntr.fuzzy_rend1 = fuzzy_flanks[1];
//            free(fuzzy_flanks);
//        }
//        else
//        {
//            vntr.fuzzy_rbeg1 = bcf_get_pos1(v) - 1;
//            vntr.fuzzy_rend1 = bcf_get_end_pos1(v) + 1;
//        }
//
//        beg1 = std::min(vntr.rbeg1-1, vntr.fuzzy_rbeg1-1);
//        end1 = std::max(vntr.rend1+1, vntr.fuzzy_rend1+1);


    }
    else if (type==VT_VNTR)
    {
        update_vntr_from_info_fields(h, v);

        vs.push_back(v);
        vntr_vs.push_back(v);
    }
}

/**
 * Constructor for a variant object that is to be composed from several variants.
 */
Variant::Variant(Variant* v1, Variant* v2)
{
    type = VT_UNDEFINED;

    chrom = v1->chrom;
    rid = v1->rid;
    pos1 = std::min(v1->pos1, v2->pos1);
    beg1 = std::min(v1->beg1, v2->beg1);
    end1 = std::max(v1->end1, v2->end1);

    no_overlapping_snps = 0;
    no_overlapping_indels = 0;
    no_overlapping_vntrs = 0;

    vs.push_back(v1->v);
    vs.push_back(v2->v);

    if (v1->type==VT_SNP)
    {
        snp_vs.push_back(v1->v);
    }
    else if (v1->type==VT_INDEL)
    {
        indel_vs.push_back(v1->v);
    }
    else if (v1->type==VT_VNTR)
    {
        vntr_vs.push_back(v1->v);
    }

    if (v2->type==VT_SNP)
    {
        snp_vs.push_back(v2->v);
    }
    else if (v2->type==VT_INDEL)
    {
        indel_vs.push_back(v2->v);
    }
    else if (v2->type==VT_VNTR)
    {
        vntr_vs.push_back(v2->v);
    }
}

/**
 * Constructor.
 */
Variant::Variant()
{
    type = VT_REF;
    v = NULL;
    no_overlapping_snps = 0;
    no_overlapping_indels = 0;
    no_overlapping_vntrs = 0;
    alleles.clear();
}

/**
 * Destructor.
 */
Variant::~Variant()
{
//    if (v) bcf_destroy(v);
//    for (uint32_t i=0; i<vs.size(); ++i)
//    {
//        if (vs[i]) bcf_destroy(vs[i]);
//    }
    clear();
};

/**
 * Clears variant information.
 */
void Variant::clear()
{
    type = VT_REF;
    alleles.clear();
    vntr.clear();
    vs.clear();
    snp_vs.clear();
    indel_vs.clear();
    vntr_vs.clear();
};

/**
 * Classifies variants.
 */
int32_t Variant::classify(bcf_hdr_t *h, bcf1_t *v)
{
    this->h = h;
    this->v = v;

    bcf_unpack(v, BCF_UN_STR);
    chrom.assign(bcf_get_chrom(h, v));
    rid = bcf_get_rid(v);
    pos1 = bcf_get_pos1(v);
    char** allele = bcf_get_allele(v);
    end1 = pos1 + strlen(allele[0]) - 1;
    int32_t n_allele = bcf_get_n_allele(v);

    uint32_t pos1 = pos1;
    int32_t pos0 = pos1-1;
    ts = 0;
    tv = 0;
    ins = 0;
    del = 0;

    clear(); // this sets the type to VT_REF by default.

    bool homogeneous_length = true;

    char* ref = allele[0];
    int32_t rlen = strlen(ref);

    //if only ref allele, skip this entire for loop
    for (size_t i=1; i<n_allele; ++i)
    {
        int32_t allele_type = VT_REF;

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
                            allele_type = VT_VNTR;
                        }
                    }
                }
                //VNTR
                else if (len==6 &&
                         allele[i][0]=='<' &&
                         allele[i][1]=='V' && allele[i][2]=='N' && allele[i][3]=='T' && allele[i][4]=='R' &&
                         allele[i][5]=='>' )
                {
                     allele_type = VT_VNTR;
                }
                //ST/d+
                else if (allele[i][0]=='<' && allele[i][1]=='S' && allele[i][2]=='T' && allele[i][len-1]=='>' )
                {
                    for (size_t j=3; j<len-1; ++j)
                    {
                        if (allele[i][j]<'0' || allele[i][j]>'9')
                        {
                            allele_type = VT_VNTR;
                        }
                    }
                }
                //STR
                else if (len==5 &&
                         allele[i][0]=='<' &&
                         allele[i][1]=='S' && allele[i][2]=='T' && allele[i][3]=='R' &&
                         allele[i][4]=='>' )
                {
                     allele_type = VT_VNTR;
                }
            }

            if (allele_type==VT_VNTR)
            {
                allele_type = VT_VNTR;
                type |= allele_type;
                alleles.push_back(Allele(allele_type));
            }
            else
            {
                allele_type = VT_SV;
                type |= allele_type;
                std::string sv_type(allele[i]);
                alleles.push_back(Allele(allele_type, sv_type));
            }
        }
        else if (strchr(allele[i],'[')||strchr(allele[i],']'))
        {
            allele_type = VT_SV;
            type |= allele_type;
            std::string sv_type("<BND>");
            alleles.push_back(Allele(allele_type, sv_type));
        }
        else if (allele[i][0]=='.' || strcmp(allele[i],allele[0])==0)
        {
            type = VT_REF;
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
                allele_type |= mlen==1 ? VT_SNP : VT_MNP;
            }

            //indel variants
            if (dlen)
            {
                allele_type |= VT_INDEL;
            }

            //clumped SNPs and MNPs
            if (diff && diff < mlen) //internal gaps
            {
                allele_type |= VT_CLUMPED;
            }

            type |= allele_type;
            alleles.push_back(Allele(type, diff, alen, dlen, mlen, ts, tv));
            ts += ts;
            tv += tv;
            ins = dlen>0?1:0;
            del = dlen<0?1:0;

            if (REF.m) free(REF.s);
            if (ALT.m) free(ALT.s);
        }
    }

    if (type==VT_VNTR)
    {
        update_vntr_from_info_fields(h, v);
    }

    //additionally define MNPs by length of all alleles
    if (!(type&(VT_VNTR|VT_SV)) && type!=VT_REF)
    {
        if (homogeneous_length && rlen>1 && n_allele>1)
        {
            type |= VT_MNP;
        }
    }

    return type;
}

/**
 * Updates VNTR related information from INFO fields.
 */
void Variant::update_vntr_from_info_fields(bcf_hdr_t *h, bcf1_t *v)
{
    vntr.motif = bcf_get_rid(v);
    char** allele = bcf_get_allele(v);
    vntr.exact_repeat_tract.assign(allele[0]);

    char *motif = NULL;
    int32_t n = 0;
    if (bcf_get_info_string(h, v, "MOTIF", &motif, &n)>0)
    {
        vntr.motif.assign(motif);
        free(motif);

        vntr.basis = vntr.get_basis(vntr.motif);
    }
    else
    {
        vntr.motif = "";
    }

    char *ru = NULL;
    n = 0;
    if (bcf_get_info_string(h, v, "RU", &ru, &n)>0)
    {
        vntr.ru.assign(ru);
        free(ru);
    }
    else
    {
        vntr.ru = "";
    }

    float *exact_motif_concordance = NULL;
    n = 0;
    if (bcf_get_info_float(h, v, "CONCORDANCE", &exact_motif_concordance, &n)>0)
    {
        vntr.exact_motif_concordance = exact_motif_concordance[0];
        free(exact_motif_concordance);
    }
    else
    {
        vntr.exact_motif_concordance = -1;
    }

    float *fuzzy_motif_concordance = NULL;
    n = 0;
    if (bcf_get_info_float(h, v, "FZ_CONCORDANCE", &fuzzy_motif_concordance, &n)>0)
    {
        vntr.fuzzy_motif_concordance = fuzzy_motif_concordance[0];
        free(fuzzy_motif_concordance);
    }
    else
    {
        vntr.fuzzy_motif_concordance = -1;
    }

    int32_t *flanks = NULL;
    n = 0;
    if (bcf_get_info_int32(h, v, "FLANKS", &flanks, &n)>0)
    {
        vntr.exact_rbeg1 = flanks[0];
        vntr.exact_rend1 = flanks[1];
        free(flanks);

        if (bcf_get_pos1(v)==vntr.exact_rbeg1 && bcf_get_end1(v)==vntr.exact_rend1)
        {
            char** allele = bcf_get_allele(v);
            vntr.exact_repeat_tract.assign(allele[0]);
        }
        else
        {
            vntr.exact_repeat_tract = "";
        }
    }
    else
    {
        vntr.exact_rbeg1 = bcf_get_pos1(v) - 1;
        vntr.exact_rend1 = bcf_get_end_pos1(v) + 1;
    }

    int32_t *fuzzy_flanks = NULL;
    n = 0;
    if (bcf_get_info_int32(h, v, "FZ_FLANKS", &fuzzy_flanks, &n)>0)
    {
        vntr.fuzzy_rbeg1 = fuzzy_flanks[0];
        vntr.fuzzy_rend1 = fuzzy_flanks[1];
        free(fuzzy_flanks);

        if (bcf_get_pos1(v)==vntr.fuzzy_rbeg1 && bcf_get_end1(v)==vntr.fuzzy_rend1)
        {
            char** allele = bcf_get_allele(v);
            vntr.fuzzy_repeat_tract.assign(allele[0]);
        }
        else
        {
            vntr.fuzzy_repeat_tract = "";
        }
    }
    else
    {
        vntr.fuzzy_rbeg1 = 0;
        vntr.fuzzy_rend1 = 0;
    }

    beg1 = std::min(vntr.exact_rbeg1-1, vntr.fuzzy_rbeg1-1);
    end1 = std::max(vntr.exact_rend1+1, vntr.fuzzy_rend1+1);
}

/**
 * Prints variant information.
 */
void Variant::print()
{
    std::cerr << "type : " << vtype2string(type) << "\n";
    std::cerr << "motif: " << vntr.motif << "\n";
    std::cerr << "rlen : " << vntr.motif.size() << "\n";

    for (int32_t i=0; i<alleles.size(); ++i)
    {
        std::cerr << "\tallele: " << i << "\n";
        std::cerr << "\t  type: " << vtype2string(alleles[i].type) << "\n";
        std::cerr << "\t  diff: " << alleles[i].diff << "\n";
        std::cerr << "\t  alen: " << alleles[i].alen << "\n";
        std::cerr << "\t  dlen: " << alleles[i].dlen << "\n";
    }
};

/**
 * Gets a string representation of the underlying VNTR by exact alignment.
 */
void Variant::get_vntr_string(kstring_t* s)
{
    s->l = 0;
    kputs(chrom.c_str(), s);
    kputc(':', s);
    kputw(vntr.exact_rbeg1, s);
    kputc(':', s);
    kputs(vntr.exact_repeat_tract.c_str(), s);
    kputc(':', s);
    kputs("<VNTR>", s);
    kputc(':', s);
    kputs(vntr.motif.c_str(), s);
};

/**
 * Gets a string representation of the underlying VNTR by fuzzy alignment.
 */
void Variant::get_fuzzy_vntr_string(kstring_t* s)
{
    s->l = 0;
    kputs(chrom.c_str(), s);
    kputc(':', s);
    kputw(vntr.fuzzy_rbeg1, s);
    kputc(':', s);
    kputs(vntr.fuzzy_repeat_tract.c_str(), s);
    kputc(':', s);
    kputs("<VNTR>", s);
    kputc(':', s);
    kputs(vntr.motif.c_str(), s);
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

    if (VTYPE & VT_VNTR)
    {
        s += (s.size()==0) ? "" : "/";
        s += "VNTR";
    }

    if (VTYPE & VT_SV)
    {
        s += (s.size()==0) ? "" : "/";
        s += "SV";
    }

    return s;
}