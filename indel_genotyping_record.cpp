/* The MIT License

   Copyright (c) 2016 Hyun Min Kang <hmkang.umich.edu> and Adrian Tan <atks@umich.edu>

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

#include "indel_genotyping_record.h"

/**
 * Constructor.
 * @v - VCF record.
 */
IndelGenotypingRecord::IndelGenotypingRecord(bcf_hdr_t *h, bcf1_t *v, int32_t nsamples, int32_t ploidy, Estimator* est)
{
    clear();

    this->h = h;
    this->rid = bcf_get_rid(v);
    this->pos1 = bcf_get_pos1(v);
    this->nsamples = nsamples;

    this->alleles = {0,0,0};
    char** tmp_alleles = bcf_get_allele(v);
    for (size_t i=0; i< bcf_get_n_allele(v); ++i) 
    {
        if (i) kputc(',', &this->alleles);
        kputs(tmp_alleles[i], &this->alleles);
        v_alleles.push_back(tmp_alleles[i]);
    }

    n_filter = 0;

    //rid = bcf_get_rid(v);
    dlen = strlen(tmp_alleles[1])-strlen(tmp_alleles[0]);
    len = abs(dlen);
    
    int32_t *flanks_pos1 = NULL;
    int32_t n = 0;
    
    if (bcf_get_info_int32(h, v, "FLANKS", &flanks_pos1, &n)>0) 
    {
        this->beg1 = flanks_pos1[0];
        this->end1 = flanks_pos1[1];
        free(flanks_pos1);
    }
    else 
    {
        this->beg1 = bcf_get_pos1(v) - 3;
        this->end1 = bcf_get_end1(v) + 3;
    }

    if (dlen>0) 
    {
        indel.append(&tmp_alleles[1][1]);
    }
    else 
    {
        indel.append(&tmp_alleles[0][1]);
    }

    if (bcf_has_filter(h, v, const_cast<char*>("overlap_snp"))==1)
        n_filter |= FILTER_MASK_OVERLAP_SNP;
    
    if (bcf_has_filter(h, v, const_cast<char*>("overlap_indel"))==1)
        n_filter |= FILTER_MASK_OVERLAP_INDEL;
    
    if (bcf_has_filter(h, v, const_cast<char*>("overlap_vntr"))==1)
        n_filter |= FILTER_MASK_OVERLAP_VNTR;

    
    pls = (uint8_t*)calloc( nsamples*3, sizeof(uint8_t) );
    ads = (uint8_t*)calloc( nsamples*3, sizeof(uint8_t) );
}

/**
 * Clears this record.
 */
void IndelGenotypingRecord::clearTemp()
{
    tmp_dp_q20 = 0;
    tmp_dp_ra = 0;
    tmp_bq_s1 = tmp_bq_s2 = 0;
    tmp_mq_s1 = tmp_mq_s2 = 0;
    tmp_cy_s1 = tmp_cy_s2 = 0;
    tmp_st_s1 = tmp_st_s2 = 0;
    tmp_al_s1 = tmp_bq_al = tmp_mq_al = tmp_cy_al = tmp_st_al = tmp_nm_al = 0;
    tmp_nm_s1 = tmp_nm_s2 = 0;
    tmp_oth_exp_q20 = tmp_oth_obs_q20 = 0;
    tmp_pls[0] = tmp_pls[1] = tmp_pls[2] = 1.;
    tmp_ads[0] = tmp_ads[1] = tmp_ads[2] = 0;
}

void IndelGenotypingRecord::clear()
{
    vtype = -1;

    pls = ads = NULL;

    bqr_num = bqr_den = 0;
    mqr_num = mqr_den = 0;
    cyr_num = cyr_den = 0;
    str_num = str_den = 0;
    nmr_num = nmr_den = 0;
    ior_num = ior_den = 0;
    nm0_num = nm0_den = 0;
    nm1_num = nm1_den = 0;
    abe_num = abe_den = 0;
    abz_num = abz_den = 0;
    ns_nref = dp_sum = max_gq = 0;
    clearTemp();
}

/**
 * Destructor.
 */
IndelGenotypingRecord::~IndelGenotypingRecord()
{
  //if (v) bcf_destroy(v);
    if ( pls ) free(pls);
    if ( ads ) free(ads);
    if ( alleles.l > 0 ) free(alleles.s);
    if ( est ) delete est;
}

/**
 * Destructor.
 */
bcf1_t* IndelGenotypingRecord::flush_variant(bcf_hdr_t* hdr)
{
    bcf1_t *nv = bcf_init();
    bcf_clear(nv);
    bcf_set_n_sample(nv, nsamples);

    bcf_set_rid(nv, rid);
    bcf_set_pos1(nv, pos1);
    bcf_update_alleles_str(hdr, nv, alleles.s);

    int32_t* gt = (int32_t*) calloc ( nsamples * 2, sizeof(int32_t) );
    int32_t* pl = (int32_t*) calloc ( nsamples * 3, sizeof(int32_t) );
    int32_t* ad = (int32_t*) calloc ( nsamples * 2, sizeof(int32_t) );
    int32_t* td = (int32_t*) calloc ( nsamples * 1, sizeof(int32_t) );
    int32_t* gq = (int32_t*) calloc ( nsamples * 1, sizeof(int32_t) );
    float MLE_HWE_AF[2];
    float MLE_HWE_GF[3];
    double gp, gp_sum, max_gp;
    int32_t best_gt;
    int32_t best_a1, best_a2;
    int32_t* pls_i;
    int32_t an = 0;
    int32_t acs[2] = {0,0};
    int32_t gcs[3] = {0,0,0};
    float afs[3];
    int32_t max_gq = 0;
    int64_t dp_sum = 0;
    int32_t ploidy = 2;

    std::copy( pls, pls + (nsamples*3), pl );
    for(int32_t i=0; i < nsamples; ++i) 
    {
        dp_sum += ( td[i] = ( (ad[2*i] = ads[3*i]) + (ad[2*i+1] = ads[3*i+1]) + ads[3*i+2] ) );
    }

    int32_t n = 0;
    int32_t adSumHet[2] = {0,0};
    est->compute_gl_af_hwe(pl, nsamples, ploidy, 2, MLE_HWE_AF, MLE_HWE_GF, n, 1e-20);

    for(int32_t i=0; i < nsamples; ++i) 
    {
        int32_t* pli = &pl[ i * 3 ];
        max_gp = gp_sum = gp = ( est->lt->pl2prob(pli[0]) * MLE_HWE_AF[0] * MLE_HWE_AF[0] );
        best_gt = 0; best_a1 = 0; best_a2 = 0;
        for(size_t l=1; l < 2; ++l) 
        {
            for(size_t m=0; m <= l; ++m) 
            {
                gp = ( est->lt->pl2prob(pli[ l*(l+1)/2 + m]) * MLE_HWE_AF[l] * MLE_HWE_AF[m] * (l == m ? 1 : 2) );
                gp_sum += gp;
                if ( max_gp < gp ) 
                {
                    max_gp = gp;
                    best_gt = l*(l+1)/2 + m;
                    best_a1 = m;
                    best_a2 = l;
                }
            }
        }

        if ( best_gt == 1 ) 
        {
            adSumHet[0] += ad[2*i];
            adSumHet[1] += ad[2*i+1];
        }

        double prob = 1.-max_gp/gp_sum;  // to calculate GQ
        if ( prob <= 3.162278e-26 )
        {
            prob = 3.162278e-26;
        }
       
        if ( prob > 1 )
        {
            prob = 1;
        }
        
        gq[i] = (int32_t)est->lt->prob2pl(prob);
    
        if ( ( best_gt > 0 ) && ( max_gq < gq[i] ) )
        {
            max_gq = gq[i];
        }
        
        gt[2*i]   = ((best_a1 + 1) << 1);
        gt[2*i+1] = ((best_a2 + 1) << 1);
        an += 2;             // still use diploid representation of chrX for now.
        ++acs[best_a1];
        ++acs[best_a2];
        ++gcs[best_gt];
    }

    for(size_t i=0; i < 2; ++i) 
    {
        afs[i] = acs[i]/(float)an;
    }

    bcf_update_format_int32(hdr, nv, "GT", gt, nsamples * 2);
    bcf_update_format_int32(hdr, nv, "GQ", gq, nsamples );
    bcf_update_format_int32(hdr, nv, "AD", ad, nsamples * 2);
    bcf_update_format_int32(hdr, nv, "DP", td, nsamples );
    bcf_update_format_int32(hdr, nv, "PL", pl, nsamples * 3);

    float avgdp = (float)dp_sum / (float)nsamples;

    nv->qual = (float) max_gq;
    bcf_update_info_float(hdr, nv, "AVGDP", &avgdp, 1);
    bcf_update_info_int32(hdr, nv, "AC", &acs[1], 1);
    bcf_update_info_int32(hdr, nv, "AN", &an, 1);
    bcf_update_info_float(hdr, nv, "AF", &afs[1], 1);
    bcf_update_info_int32(hdr, nv, "GC", gcs, 3);
    bcf_update_info_int32(hdr, nv, "GN", &nsamples, 1);

    if (n) 
    {
        float* MLE_HWE_AF_PTR = &MLE_HWE_AF[1];
        bcf_update_info_float(hdr, nv, "HWEAF", MLE_HWE_AF_PTR, 1);
    }

    // calculate the allele frequencies under HWD
    float MLE_AF[2];
    float MLE_GF[3];
    n = 0;
    est->compute_gl_af(pl, nsamples, ploidy, 2, MLE_AF, MLE_GF,  n, 1e-20);
    if (n) 
    {
        //float* MLE_AF_PTR = &MLE_AF[1];
        //bcf_update_info_float(odw->hdr, nv, "HWDAF", MLE_AF_PTR, n_alleles-1);
        bcf_update_info_float(hdr, nv, "HWDGF", &MLE_GF, 3);
    }

    float fic = 0;
    n = 0;
    est->compute_gl_fic(pl, nsamples, ploidy, MLE_HWE_AF, 2, MLE_GF, fic, n);
    if ( std::isnan((double)fic) ) fic = 0;
    if (n) 
    {
        bcf_update_info_float(hdr, nv, "IBC", &fic, 1);
    }
    
    // calculate the LRT statistics related to HWE
    float lrts;
    float logp;
    int32_t df;
    n = 0;
    est->compute_hwe_lrt(pl, nsamples, ploidy, 2, MLE_HWE_GF, MLE_GF, n, lrts, logp, df);
    if (n) 
    {
        if ( fic < 0 ) logp = 0-logp;
        {
            bcf_update_info_float(hdr, nv, "HWE_SLP", &logp, 1);
        }
    }

    float abe = (adSumHet[0] + 0.5)/(adSumHet[0] + adSumHet[1] + 1.0);
    bcf_update_info_float(hdr, nv, "ABE", &abe, 1);

    // update filter
    if ( n_filter & FILTER_MASK_OVERLAP_SNP )
    bcf_add_filter(hdr, nv, bcf_hdr_id2int(hdr, BCF_DT_ID, "overlap_snp"));
    if ( n_filter & FILTER_MASK_OVERLAP_INDEL )
    bcf_add_filter(hdr, nv, bcf_hdr_id2int(hdr, BCF_DT_ID, "overlap_indel"));
    if ( n_filter & FILTER_MASK_OVERLAP_VNTR )
    bcf_add_filter(hdr, nv, bcf_hdr_id2int(hdr, BCF_DT_ID, "overlap_vntr"));

    //bcf_update_info_int32(odw->hdr, nv, "NS_NREF", &v_ns_nrefs[k], 1);
    abe_num /= (abe_den+1e-6); bcf_update_info_float(hdr, nv, "ABE",  &abe_num, 1);
    abz_num /= sqrt(abz_den+1e-6); bcf_update_info_float(hdr, nv, "ABZ",  &abz_num, 1);
    bqr_num /= sqrt(bqr_den+1e-6); bcf_update_info_float(hdr, nv, "BQZ", &bqr_num, 1);
    mqr_num /= sqrt(mqr_den+1e-6); bcf_update_info_float(hdr, nv, "MQZ", &mqr_num, 1);
    cyr_num /= sqrt(cyr_den+1e-6); bcf_update_info_float(hdr, nv, "CYZ", &cyr_num, 1);
    str_num /= sqrt(str_den+1e-6); bcf_update_info_float(hdr, nv, "STZ", &str_num, 1);
    nmr_num /= sqrt(nmr_den+1e-6); bcf_update_info_float(hdr, nv, "NMZ", &nmr_num, 1);
    ior_num = log(ior_num/ior_den+1e-6)/log(10.); bcf_update_info_float(hdr, nv, "IOR", &ior_num, 1);
    nm1_num /= (nm1_den+1e-6); bcf_update_info_float(hdr, nv, "NM1", &nm1_num, 1);
    nm0_num /= (nm0_den+1e-6); bcf_update_info_float(hdr, nv, "NM0", &nm0_num, 1);

    //odw->write(nv);
    //bcf_destroy(nv);
    free(gt);
    free(gq);
    free(pl);
    free(ad);
    free(td);

    free(pls); pls = NULL;
    free(ads); ads = NULL;

    delete est;
    free(alleles.s);

    return nv;
}

void IndelGenotypingRecord::flush_sample(int32_t sampleIndex)
{
    uint8_t* p_pls = &pls[sampleIndex*3];
    uint8_t* p_ads = &ads[sampleIndex*3];

    int32_t imax = ( tmp_pls[0] > tmp_pls[1] ) ? ( tmp_pls[0] > tmp_pls[2] ? 0 : 2 ) : ( tmp_pls[1] > tmp_pls[2] ? 1 : 2);
    for(int32_t i=0; i < 3; ++i) 
    {
        uint32_t l = est->lt->prob2pl(tmp_pls[i]/tmp_pls[imax]);
        p_pls[i] = ((l > 255) ? 255 : l);
        p_ads[i] = ((tmp_ads[i] > 255) ? 255 : (uint8_t)tmp_ads[i]);
    }

    float sqrt_dp_ra = sqrt((float)tmp_dp_ra);
    float ior = (float)(tmp_oth_obs_q20 / (tmp_oth_exp_q20 + 1e-6));
    float nm1 = tmp_al_s1 == 0 ? 0 : tmp_nm_al / (float)tmp_al_s1;
    float nm0 = (tmp_dp_ra - tmp_al_s1) == 0 ? 0 : (tmp_nm_s1-tmp_nm_al) / (float)(tmp_dp_ra - tmp_al_s1);
    float w_dp_ra  = log(tmp_dp_ra+1.); //sqrt(dp_ra);
    float w_dp_q20 = log(tmp_dp_q20+1.); //sqrt(dp_q20);
    float w_al_s1  = log(tmp_al_s1+1.); //sqrt(al_s1);
    float w_ref_s1 = log(tmp_dp_ra - tmp_al_s1+1.);

    if ( p_pls[1] == 0 ) 
    { // het genotypes
        abe_num += (w_dp_ra * (tmp_dp_ra - tmp_al_s1 + 0.05) / (double)(tmp_dp_ra + 0.1));
        abe_den += w_dp_ra;
    
        // E(r) = 0.5(r+a) V(r) = 0.25(r+a)
        abz_num += w_dp_ra * (tmp_dp_ra - tmp_al_s1 - tmp_dp_ra*0.5)/sqrt(0.25 * tmp_dp_ra + 1e-3);
        abz_den += (w_dp_ra * w_dp_ra);
    
        float bqr = sqrt_dp_ra * Estimator::compute_correlation( tmp_dp_ra, tmp_bq_al, tmp_bq_s1, tmp_bq_s2, tmp_al_s1, tmp_al_s1, .1 );
        float mqr = sqrt_dp_ra * Estimator::compute_correlation( tmp_dp_ra, tmp_mq_al, tmp_mq_s1, tmp_mq_s2, tmp_al_s1, tmp_al_s1, .1 );
        float cyr = sqrt_dp_ra * Estimator::compute_correlation_f( tmp_dp_ra, tmp_cy_al, tmp_cy_s1, tmp_cy_s2, (float)tmp_al_s1, (float)tmp_al_s1, .1 );
        float str = sqrt_dp_ra * Estimator::compute_correlation( tmp_dp_ra, tmp_st_al, tmp_st_s1, tmp_st_s1, tmp_al_s1, tmp_al_s1, .1 );
        float nmr = sqrt_dp_ra * Estimator::compute_correlation( tmp_dp_ra, tmp_nm_al, tmp_nm_s1, tmp_nm_s2, tmp_al_s1, tmp_al_s1, .1 );
    
        // Use Stouffer's method to combine the z-scores, but weighted by log of sample size
        bqr_num += (bqr * w_dp_ra); bqr_den += (w_dp_ra * w_dp_ra);
        mqr_num += (mqr * w_dp_ra); mqr_den += (w_dp_ra * w_dp_ra);
        cyr_num += (cyr * w_dp_ra); cyr_den += (w_dp_ra * w_dp_ra);
        str_num += (str * w_dp_ra); str_den += (w_dp_ra * w_dp_ra);
        nmr_num += (nmr * w_dp_ra); nmr_den += (w_dp_ra * w_dp_ra);
    }

    ior_num += (ior * w_dp_q20); ior_den += w_dp_q20;
    nm1_num += (nm1 * w_al_s1);  nm1_den += w_al_s1;
    nm0_num += (nm0 * w_ref_s1); nm0_den += w_ref_s1;

    clearTemp();
}

/**
 * Adds an allele based with collected sufficient statistics.
 */
void IndelGenotypingRecord::add_allele(double contam, int32_t allele, uint8_t mapq, bool fwd, uint32_t q, int32_t cycle, uint32_t nm)
{
    double pe = est->lt->pl2prob(q);
    double pm = 1 - pe;

    if (q>40)
    {
        q = 40;
    }

    if ( allele == 0 )
    {
        ++tmp_ads[0];
        tmp_pls[0] *= ( pm * (1-contam) + pe * contam / 3 );
        tmp_pls[1] *= ( pm / 2 + pe / 6 );
        tmp_pls[2] *= ( pm * contam + pe * (1-contam) / 3 );
    }
    else if ( allele > 0 ) // currently, bi-allelic only
    {
        ++tmp_ads[1];
        tmp_pls[0] *= ( pm * contam + pe * (1-contam) / 3 );
        tmp_pls[1] *= ( pm / 2 + pe / 6 );
        tmp_pls[2] *= ( pm * (1-contam) + pe * contam / 3 );
    }
    else
    {
        ++tmp_ads[2];
    }
    double sump = tmp_pls[0] + tmp_pls[1] + tmp_pls[2] + 1e-300;
    tmp_pls[0] /= sump;
    tmp_pls[1] /= sump;
    tmp_pls[2] /= sump;

    if ( allele >= 0 )
    {
        if ( q > 20 )
        {
          tmp_oth_exp_q20 += (est->lt->pl2prob(q) * 2. / 3.);
          ++tmp_dp_q20;
        }

        float log_td = (cycle > 0) ? 0-logf((float)cycle) : 0;

        ++tmp_dp_ra;
        tmp_bq_s1 += q;
        tmp_bq_s2 += (q*q);
        tmp_mq_s1 += mapq;
        tmp_mq_s2 += (mapq * mapq);
        tmp_cy_s1 += log_td;
        tmp_cy_s2 += (log_td * log_td);
        tmp_st_s1 += fwd;

        if ( allele > 0 )
        {
            ++tmp_al_s1;
            tmp_bq_al += q;
            tmp_mq_al += mapq;
            tmp_cy_al += log_td;
            tmp_st_al += fwd;
            tmp_nm_al += (nm-1);
            tmp_nm_s1 += (nm-1);
            tmp_nm_s2 += (nm-1)*(nm-1);
        }
        else
        {
          tmp_nm_s1 += nm;
          tmp_nm_s2 += (nm * nm);
        }
    }
    else
    {
        if (q>20)
        {
          tmp_oth_exp_q20 += (est->lt->pl2prob(q) * 2. / 3.);
          ++tmp_oth_obs_q20;
          ++tmp_dp_q20;
        }
    }
}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 */
void IndelGenotypingRecord::process_read(AugmentedBAMRecord& as, int32_t sampleIndex, double contam)
{
    if (v_alleles.size()==2)
    {
        if (as.beg1 <= beg1 && end1 <= as.end1)
        {
            bam1_t *s = as.s;

            char strand = bam_is_rev(s) ? 'R' : 'F';
            int32_t allele = 0;
            //uint32_t bpos1 = bam_get_pos1(s);
            uint8_t* seq = bam_get_seq(s);
            uint8_t* qual = bam_get_qual(s);
            uint32_t rlen = bam_get_l_qseq(s);
            uint8_t mapq = bam_get_mapq(s);

            uint32_t q = len*30;
            uint32_t cycle = 10;

            std::vector<uint32_t>& aug_cigar = as.aug_cigar;
            std::vector<std::string>& aug_ref = as.aug_ref;
            std::vector<std::string>& aug_alt = as.aug_alt;

            int32_t vpos1 = pos1;

            int32_t cpos1 = bam_get_pos1(s);
            int32_t rpos0 = 0;

            for (uint32_t i=0; i<aug_cigar.size(); ++i)
            {
                uint32_t oplen = bam_cigar_oplen(aug_cigar[i]);
                char opchr = bam_cigar_opchr(aug_cigar[i]);

                if (opchr=='S')
                {
                    rpos0 += oplen;
                }
                else if (opchr=='=')
                {
                    if (vpos1>=cpos1 && vpos1<=(cpos1+oplen-1))
                    {
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                    }
                    cpos1 += oplen;
                    rpos0 += oplen;
                }
                else if (opchr=='X')
                {
                    if (cpos1-1==vpos1)
                    {
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                    }
                    ++cpos1;
                    ++rpos0;
                }
                else if (opchr=='I')
                {
                    if (dlen>0 && cpos1-1==vpos1)
                    {
                        if (indel==aug_alt[i])
                        {
                            q = len*30;
                            allele = 1;
                        }
                        else
                        {
                            q = abs(len-(int32_t)aug_ref[i].size())*30;
                            allele = -1;
                        }

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }
                    else if (dlen<0 && cpos1-1==vpos1)
                    {
                        q = 30;
                        allele = -3;
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }

                    rpos0 += oplen;
                }
                else if (opchr=='D')
                {
                    if (dlen<0 && cpos1-1==vpos1)
                    {
                        if (indel==aug_ref[i])
                        {
                            q = len*30;
                            allele = 1;
                        }
                        else
                        {
                            q = abs(len-(int32_t)aug_ref[i].size())*30;
                            allele = -1;
                        }

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }
                    else if (dlen>0 && cpos1-1==vpos1)
                    {
                        q = 30;
                        allele = -2;

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }

                }
                else
                {
                    std::cerr << "unrecognized cigar state " << opchr << "\n";
                }
            }

            add_allele( contam, allele, mapq, strand == 'F', q, rand() % 75, as.no_mismatches );
        }
        else
        {
            bam1_t *s = as.s;
            uint8_t mapq = bam_get_mapq(s);

            add_allele( contam, -1, mapq, bam_is_rev(as.s) ? false : true, 20, rand() % 75, as.no_mismatches );
        }
    }
    else //multiallelic
    {
        if (as.beg1 <= beg1 && end1 <= as.end1)
        {
            bam1_t *s = as.s;

            char strand = bam_is_rev(s) ? 'R' : 'F';
            int32_t allele = 0;
            //uint32_t bpos1 = bam_get_pos1(s);
            uint8_t* seq = bam_get_seq(s);
            uint8_t* qual = bam_get_qual(s);
            uint32_t rlen = bam_get_l_qseq(s);
            uint8_t mapq = bam_get_mapq(s);

            uint32_t q = len*30;
            uint32_t cycle = 10;

            std::vector<uint32_t>& aug_cigar = as.aug_cigar;
            std::vector<std::string>& aug_ref = as.aug_ref;
            std::vector<std::string>& aug_alt = as.aug_alt;

            int32_t vpos1 = pos1;

            int32_t cpos1 = bam_get_pos1(s);
            int32_t rpos0 = 0;

            for (uint32_t i=0; i<aug_cigar.size(); ++i)
            {
                uint32_t oplen = bam_cigar_oplen(aug_cigar[i]);
                char opchr = bam_cigar_opchr(aug_cigar[i]);

                if (opchr=='S')
                {
                    rpos0 += oplen;
                }
                else if (opchr=='=')
                {
                    if (vpos1>=cpos1 && vpos1<=(cpos1+oplen-1))
                    {
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                    }
                    cpos1 += oplen;
                    rpos0 += oplen;
                }
                else if (opchr=='X')
                {
                    if (cpos1-1==vpos1)
                    {
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                    }
                    ++cpos1;
                    ++rpos0;
                }
                else if (opchr=='I')
                {
                    if (dlen>0 && cpos1-1==vpos1)
                    {
                        if (indel==aug_alt[i])
                        {
                            q = len*30;
                            allele = 1;
                        }
                        else
                        {
                            q = abs(len-(int32_t)aug_ref[i].size())*30;
                            allele = -1;
                        }

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }
                    else if (dlen<0 && cpos1-1==vpos1)
                    {
                        q = 30;
                        allele = -3;
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }

                    rpos0 += oplen;
                }
                else if (opchr=='D')
                {
                    if (dlen<0 && cpos1-1==vpos1)
                    {
                        if (indel==aug_ref[i])
                        {
                            q = len*30;
                            allele = 1;
                        }
                        else
                        {
                            q = abs(len-(int32_t)aug_ref[i].size())*30;
                            allele = -1;
                        }

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }
                    else if (dlen>0 && cpos1-1==vpos1)
                    {
                        q = 30;
                        allele = -2;

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }

                }
                else
                {
                    std::cerr << "unrecognized cigar state " << opchr << "\n";
                }
            }

            add_allele( contam, allele, mapq, strand == 'F', q, rand() % 75, as.no_mismatches );
        }
        else
        {
            bam1_t *s = as.s;
            uint8_t mapq = bam_get_mapq(s);

            add_allele( contam, -1, mapq, bam_is_rev(as.s) ? false : true, 20, rand() % 75, as.no_mismatches );
        }
    }
}
