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
    this->h = h;
    this->v = v;

    type = classify(h, v);

    chrom = bcf_get_chrom(h, v);
    rid = bcf_get_rid(v);
    pos1 = bcf_get_pos1(v);

    no_overlapping_snps = 0;
    no_overlapping_indels = 0;
    no_overlapping_vntrs = 0;

    is_new_multiallelic =  false;

    //attempts to update relevant information on variants
    if (type==VT_SNP)
    {
        beg1 = bcf_get_pos1(v);
        end1 = bcf_get_pos1(v);
    }
    else if (type==VT_INDEL)
    {
        beg1 = bcf_get_pos1(v);
        end1 = bcf_get_end1(v);
    }
    else if (type & (VT_SNP|VT_MNP|VT_INDEL|VT_CLUMPED))
    {
        beg1 = bcf_get_pos1(v);       
        end1 = bcf_get_info_int(h, v, "END", 0);
        if (!end1) end1 = bcf_get_end1(v);
    }
    else if (type==VT_VNTR)
    {
        beg1 = bcf_get_pos1(v);       
        end1 = bcf_get_info_int(h, v, "END", 0);
        if (!end1) end1 = bcf_get_end1(v);
        
        update_vntr_from_info_fields(h, v);

        vs.push_back(v);
        vntr_vs.push_back(v);
    }
    else if (type==VT_SV)
    {
        beg1 = bcf_get_pos1(v);
        end1 = bcf_get_info_int(h, v, "END", 0);
        if (!end1) end1 = strlen(bcf_get_allele(v)[0]);
    }
    else
    {
        std::cerr <<  "unexpected type in variant construction\n";
        print();
        exit(1);
    }
      
}

/**
 * Constructor for a variant object that is to be composed from several variants.
 * ??????? WHY DO THIS????????
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
    clear();
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
    
    h = NULL;
    v = NULL;
    
    is_new_multiallelic = false;
    is_involved_in_a_multiallelic = false;
    associated_new_multiallelic = NULL;
    updated_multiallelic = false;
    
    chrom.clear();
    rid = 0;
    pos1 = 0;
    beg1 = 0;
    end1 = 0;
    ts = 0;
    tv = 0;
    ins = 0;
    del = 0;
    max_dlen = 0;
    min_dlen = 0;

    no_overlapping_snps = 0;
    no_overlapping_indels = 0;
    no_overlapping_vntrs = 0;
    
    contains_N = false;
    alleles.clear();
    vntr.clear();
    vs.clear();
    snp_vs.clear();
    indel_vs.clear();
    vntr_vs.clear();
    consolidated_vntr_vs.clear();
};

/**
 * Classifies variants.
 */
int32_t Variant::classify(bcf_hdr_t *h, bcf1_t *v)
{
    clear();
    
    this->h = h;
    this->v = v;

    bcf_unpack(v, BCF_UN_STR);
    chrom.assign(bcf_get_chrom(h, v));
    rid = bcf_get_rid(v);
    pos1 = bcf_get_pos1(v);
    end1 = bcf_get_end1(v);
    char** allele = bcf_get_allele(v);
    int32_t n_allele = bcf_get_n_allele(v);
    int32_t pos0 = pos1-1;

    bool homogeneous_length = true;
    char* ref = allele[0];
    int32_t rlen = strlen(ref);

    if (strchr(ref, 'N'))
    {
        contains_N = true;
    }

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
                //STR
                else if (len==5 &&
                         allele[i][0]=='<' &&
                         allele[i][1]=='S' && allele[i][2]=='T' && allele[i][3]=='R' &&
                         allele[i][4]=='>' )
                {
                     allele_type = VT_VNTR;
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
        //checks for chromosomal breakpoints
        else if (strchr(allele[i],'[')||strchr(allele[i],']'))
        {
            allele_type = VT_SV;
            type |= allele_type;
            std::string sv_type("<BND>");
            alleles.push_back(Allele(allele_type, sv_type));
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
                contains_N = true;
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
void Variant::update_vntr_from_info_fields()
{    
    vntr.motif = bcf_get_rid(v);
    char** allele = bcf_get_allele(v);
//    vntr.exact_repeat_tract.assign(allele[0]);
//   std::string tags[16] = {"MOTIF", "RU", "BASIS", "MLEN", "BLEN", "REPEAT_TRACT", "COMP", "ENTROPY", "ENTROPY2", "KL_DIVERGENCE", "KL_DIVERGENCE2", "RL", "LL", "RU_COUNTS", "SCORE", "TRF_SCORE"};
    
    vntr.motif = bcf_get_info_str(h, v, "MOTIF");
    vntr.ru = bcf_get_info_str(h, v, "RU");
    vntr.basis = bcf_get_info_str(h, v, "BASIS");
    if (vntr.basis=="") vntr.basis = VNTR::get_basis(vntr.motif);
    vntr.mlen = vntr.motif.size();    
    vntr.blen = (int32_t) vntr.basis.size();
    std::vector<int32_t> i_vec = bcf_get_info_int_vec(h, v, "REPEAT_TRACT", 2, 0);
    vntr.beg1 = i_vec[0];
    vntr.end1 = i_vec[1];
    i_vec = bcf_get_info_int_vec(h, v, "COMP", 4, 0);
    vntr.comp[0] = i_vec[0];
    vntr.comp[1] = i_vec[1];
    vntr.comp[2] = i_vec[2];
    vntr.comp[3] = i_vec[3];
    vntr.entropy = bcf_get_info_flt(h, v, "ENTROPY");
    vntr.entropy2 = bcf_get_info_flt(h, v, "ENTROPY2");
    vntr.kl_divergence = bcf_get_info_flt(h, v, "KL_DIVERGENCE");
    vntr.kl_divergence2 = bcf_get_info_flt(h, v, "KL_DIVERGENCE2");
    vntr.rl = bcf_get_info_int(h, v, "RL");
    vntr.ll = bcf_get_info_int(h, v, "LL");
    i_vec = bcf_get_info_int_vec(h, v, "RU_COUNTS", 2, 0);
    vntr.no_perfect_ru = i_vec[0];
    vntr.no_ru = i_vec[1];
    vntr.score = bcf_get_info_flt(h, v, "SCORE");
    vntr.trf_score = bcf_get_info_int(h, v, "TRF_SCORE");

    vntr.exact_motif = bcf_get_info_str(h, v, "EX_MOTIF");
    vntr.exact_ru = bcf_get_info_str(h, v, "EX_RU");
    vntr.exact_basis = bcf_get_info_str(h, v, "EX_BASIS");
    vntr.exact_mlen = (int32_t) vntr.exact_motif.size();
    vntr.exact_blen = (int32_t) vntr.exact_basis.size();
    i_vec = bcf_get_info_int_vec(h, v, "EX_REPEAT_TRACT", 2, 0);
    vntr.exact_beg1 = i_vec[0];
    vntr.exact_end1 = i_vec[1];
    i_vec = bcf_get_info_int_vec(h, v, "EX_COMP", 4, 0);
    vntr.exact_comp[0] = i_vec[0];
    vntr.exact_comp[1] = i_vec[1];
    vntr.exact_comp[2] = i_vec[2];
    vntr.exact_comp[3] = i_vec[3];
    vntr.exact_entropy = bcf_get_info_flt(h, v, "EX_ENTROPY");
    vntr.exact_entropy2 = bcf_get_info_flt(h, v, "EX_ENTROPY2");
    vntr.exact_kl_divergence = bcf_get_info_flt(h, v, "EX_KL_DIVERGENCE");
    vntr.exact_kl_divergence2 = bcf_get_info_flt(h, v, "EX_KL_DIVERGENCE2");
    vntr.exact_rl = bcf_get_info_int(h, v, "EX_RL");
    vntr.exact_ll = bcf_get_info_int(h, v, "EX_LL");
    i_vec = bcf_get_info_int_vec(h, v, "EX_RU_COUNTS", 2, 0);
    vntr.exact_no_perfect_ru = i_vec[0];
    vntr.exact_no_ru = i_vec[1];
    vntr.exact_score = bcf_get_info_flt(h, v, "EX_SCORE");
    vntr.exact_trf_score = bcf_get_info_int(h, v, "EX_TRF_SCORE");   
    
    vntr.fuzzy_motif = bcf_get_info_str(h, v, "FZ_MOTIF");
    vntr.fuzzy_ru = bcf_get_info_str(h, v, "FZ_RU");
    vntr.fuzzy_basis = bcf_get_info_str(h, v, "FZ_BASIS");
    vntr.fuzzy_mlen = (int32_t) vntr.fuzzy_motif.size();
    vntr.fuzzy_blen = (int32_t) vntr.fuzzy_basis.size();
    i_vec = bcf_get_info_int_vec(h, v, "FZ_REPEAT_TRACT", 2, 0);
    vntr.fuzzy_beg1 = i_vec[0];
    vntr.fuzzy_end1 = i_vec[1];
    i_vec = bcf_get_info_int_vec(h, v, "FZ_COMP", 4, 0);
    vntr.fuzzy_comp[0] = i_vec[0];
    vntr.fuzzy_comp[1] = i_vec[1];
    vntr.fuzzy_comp[2] = i_vec[2];
    vntr.fuzzy_comp[3] = i_vec[3];
    vntr.fuzzy_entropy = bcf_get_info_flt(h, v, "FZ_ENTROPY");
    vntr.fuzzy_entropy2 = bcf_get_info_flt(h, v, "FZ_ENTROPY2");
    vntr.fuzzy_kl_divergence = bcf_get_info_flt(h, v, "FZ_KL_DIVERGENCE");
    vntr.fuzzy_kl_divergence2 = bcf_get_info_flt(h, v, "FZ_KL_DIVERGENCE2");
    vntr.fuzzy_rl = bcf_get_info_int(h, v, "FZ_RL");
    vntr.fuzzy_ll = bcf_get_info_int(h, v, "FZ_LL");
    i_vec = bcf_get_info_int_vec(h, v, "FZ_RU_COUNTS", 2, 0);
    vntr.fuzzy_no_perfect_ru = i_vec[0];
    vntr.fuzzy_no_ru = i_vec[1];
    vntr.fuzzy_score = bcf_get_info_flt(h, v, "FZ_SCORE");
    vntr.fuzzy_trf_score = bcf_get_info_int(h, v, "FZ_TRF_SCORE");
}

/**
 * Updates VNTR related information from INFO fields.
 */
void Variant::update_vntr_from_info_fields(bcf_hdr_t *h, bcf1_t *v)
{
    this->h = h;
    this->v = v;

    update_vntr_from_info_fields();
}

/**
 * Updates a bcf1_t object with VNTR information.
 */
void Variant::update_bcf_with_vntr_info(bcf_hdr_t *h, bcf1_t *v)
{
    bcf_update_info_flag(h, v, "FUZZY", NULL, 1);

    //VNTR position and sequences
    bcf_set_pos1(v, vntr.fuzzy_beg1);
    kstring_t s = {0,0,0};
    s.l = 0;
    kputs(vntr.fuzzy_repeat_tract.c_str(), &s);
    kputc(',', &s);
    kputs("<VNTR>", &s);
    bcf_update_alleles_str(h, v, s.s);

    //VNTR motif
    bcf_update_info_string(h, v, "MOTIF", vntr.motif.c_str());
    bcf_update_info_string(h, v, "RU", vntr.ru.c_str());

    //VNTR characteristics
    bcf_update_info_float(h, v, "FZ_CONCORDANCE", &vntr.fuzzy_score, 1);
    bcf_update_info_float(h, v, "FZ_RL", &vntr.fuzzy_rl, 1);
    bcf_update_info_float(h, v, "FZ_LL", &vntr.fuzzy_ll, 1);
    int32_t flank_pos1[2] = {vntr.exact_beg1-1, vntr.exact_end1+1};
    bcf_update_info_int32(h, v, "FLANKS", &flank_pos1, 2);

    //flank positions
    int32_t fuzzy_flank_pos1[2] = {vntr.fuzzy_beg1-1, vntr.fuzzy_end1+1};
    bcf_update_info_int32(h, v, "FZ_FLANKS", &fuzzy_flank_pos1, 2);
    int32_t ru_count[2] = {vntr.fuzzy_no_perfect_ru, vntr.fuzzy_no_ru};
    bcf_update_info_int32(h, v, "FZ_RU_COUNTS", &ru_count, 2);
    
    if (vntr.is_large_repeat_tract) bcf_update_info_flag(h, v, "LARGE_REPEAT_REGION", NULL, 1);

    if (s.m) free(s.s);
}

/**
 * Prints variant information.
 */
void Variant::print()
{
    std::cerr << "type  : " << vtype2string(type) << "\n";
    std::cerr << "\n";
    if (h!=NULL && v!=NULL) bcf_print_liten(h, v); 
    std::cerr << "chrom : " << chrom << "\n";
    std::cerr << "rid   : " << rid << "\n";
    std::cerr << "beg1  : " << beg1 << "\n";
    std::cerr << "end1  : " << end1 << "\n";

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
 * Gets a string representation of the variant.
 */
std::string Variant::get_variant_string()
{
    kstring_t var = {0,0,0};
    bcf_unpack(v, BCF_UN_STR);
    var.l = 0;
    kputs(bcf_get_chrom(h, v), &var);
    kputc(':', &var);
    kputw(bcf_get_pos1(v), &var);
    kputc(':', &var);
    for (size_t i=0; i<bcf_get_n_allele(v); ++i)
    {
        if (i) kputc('/', &var);
        kputs(bcf_get_alt(v, i), &var);
    }
            
    std::string str(var.s);
    
    if (var.m) free(var.s);
    
    return str;
} 

/**
 * Gets a string representation of the underlying VNTR by exact alignment.
 */
void Variant::get_vntr_string(kstring_t* s)
{
    s->l = 0;
    kputs(chrom.c_str(), s);
    kputc(':', s);
    kputw(vntr.exact_beg1, s);
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
    kputw(vntr.fuzzy_beg1, s);
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