/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#include "estimator.h"

/**
 * Computes allele frequencies using hard calls.
 *
 * @gts        - genotypes
 * @no_samples - number of samples
 * @ploidy     - ploidy
 * @no_alleles - number of alleles
 * @AC         - alternate allele counts
 * @AN         - total number of allele counts
 * @AF         - alternate allele frequency
 * @GC         - genotype counts
 * @GN         - total number of genotype counts
 * @GF         - genotype frequency
 * @NS         - number of samples with data
 */
void Estimator::compute_af(int32_t *gts, int32_t no_samples, int32_t ploidy,
                int32_t no_alleles, int32_t *AC, int32_t& AN, float *AF,
                int32_t *GC,  int32_t& GN, float *GF, int32_t& NS)
{
    int32_t iter = 0;

    if (no_alleles==2 && ploidy==2)
    {
        NS = 0;
        int32_t no_genotypes = 3;
        AC[0] = 0;
        AC[1] = 0;

        GC[0] = 0;
        GC[1] = 0;
        GC[2] = 0;

        GN=0;

        for (size_t k=0; k<no_samples; ++k)
        {
            size_t offset = k*ploidy;

            int32_t g1 = bcf_gt_allele(gts[offset]);
            if (g1>=0)
            {
                ++AC[g1];
                ++AN;
            }
            int32_t g2 = bcf_gt_allele(gts[offset+1]);
            if (g2>=0)
            {
                ++AC[g2];
                ++AN;
            }

            if (g1>=0 && g2>=0)
            {
                ++GC[g1+g2];
                ++GN;
            }

            if (g1>=0||g2>=0) ++NS;
        }

        if (!NS)
        {
            return;
        }

        AF[0] = (float)AC[0]/AN;
        AF[1] = (float)AC[1]/AN;

        GF[0] = (float)GC[0]/GN;
        GF[1] = (float)GC[1]/GN;
        GF[2] = (float)GC[2]/GN;
    }
    else
    {
        NS = 0;
        AN = 0;
        GN = 0;
        int32_t no_genotypes = bcf_an2gn(no_alleles);
        for (size_t i=0; i<no_alleles; ++i)
        {
            AC[i] = 0;
        }

        for (size_t i=0; i<no_genotypes; ++i)
        {
            GC[i] = 0;
        }

        int32_t gt_indiv[ploidy];

        for (size_t k=0; k<no_samples; ++k)
        {
            size_t offset = k*ploidy;

            int32_t last_AN = AN;
            for (size_t i=0; i<ploidy; ++i)
            {
                gt_indiv[i] = bcf_gt_allele(gts[offset+i]);

                if (gt_indiv[i]>=0)
                {
                    ++AC[gt_indiv[i]];
                    ++AN;
                }
            }

            if (last_AN<AN) ++NS;

            if (ploidy==2 && gt_indiv[0]>=0 && gt_indiv[1]>=0)
            {
                if (gt_indiv[1]<gt_indiv[0])
                {
                    gt_indiv[0] += gt_indiv[1];
                    gt_indiv[1] = gt_indiv[0] - gt_indiv[1];
                    gt_indiv[0] -= gt_indiv[1];
                }

                ++GC[bcf_alleles2gt(gt_indiv[0],gt_indiv[1])];
                ++GN;
            }
        }

        for (size_t i=0; i<no_alleles; ++i)
        {
            AF[i] = (float)AC[i]/AN;
        }

        for (size_t i=0; i<no_genotypes; ++i)
        {
            GF[i] = (float)GC[i]/GN;
        }
    }
}

/**
 * Computes allele frequencies using EM algorithm from genotype likelihoods
 * under assumption of Hardy-Weinberg Equilibrium.
 *
 * @pls        - PHRED genotype likelihoods
 * @no_samples - number of samples
 * @ploidy     - ploidy
 * @no_alleles - number of alleles
 * @MLE_HWE_AF - estimated AF
 * @MLE_HWE_GF - estimated GF
 * @n          - effective sample size
 * @e          - error
 */
void Estimator::compute_gl_af_hwe(int32_t *pls, int32_t no_samples, int32_t ploidy,
                int32_t no_alleles, float *MLE_HWE_AF, float *MLE_HWE_GF, int32_t& n,
                double e)
{
    int32_t iter = 0;

    if (ploidy!=2)
    {
        return;
    }

    if (no_alleles==2)
    {
        n = 0;
        int32_t imap[no_samples];

        for (size_t i=0; i<no_samples; ++i)
        {
            if (pls[i*3]!=bcf_int32_missing)
            {
                imap[n] = i;
                ++n;
            }
        }

        if (!n)
        {
            return;
        }

        float af[2] = {0.5, 0.5};
        float gf[3];
        float gf_indiv[3];

        float mse = e+1;
        float diff = 0;
        while (mse>e && iter<50)
        {
            gf[0] = af[0]*af[0];
            gf[1] = 2*af[0]*af[1];
            gf[2] = af[1]*af[1];

            MLE_HWE_AF[0] = 0;
            MLE_HWE_AF[1] = 0;

            for (size_t i=0; i<n; ++i)
            {
                size_t offset = imap[i]*3;

                float prob_data = (gf_indiv[0] = gf[0]*lt->pl2prob(pls[offset]));
                prob_data += (gf_indiv[1] = gf[1]*lt->pl2prob(pls[offset+1]));
                prob_data += (gf_indiv[2] = gf[2]*lt->pl2prob(pls[offset+2]));

                gf_indiv[0] /= prob_data;
                gf_indiv[1] /= prob_data;
                gf_indiv[2] /= prob_data;

                MLE_HWE_AF[0] += gf_indiv[0] + 0.5*gf_indiv[1];
                MLE_HWE_AF[1] += gf_indiv[2] + 0.5*gf_indiv[1];
            }

            MLE_HWE_AF[0] /= n;
            MLE_HWE_AF[1] /= n;

            diff = (af[0]-MLE_HWE_AF[0]);
            mse = diff*diff;
            diff = (af[1]-MLE_HWE_AF[1]);
            mse += diff*diff;

            af[0] = MLE_HWE_AF[0];
            af[1] = MLE_HWE_AF[1];

            ++iter;
        }

        MLE_HWE_GF[0] = MLE_HWE_AF[0]*MLE_HWE_AF[0];
        MLE_HWE_GF[1] = 2*MLE_HWE_AF[0]*MLE_HWE_AF[1];
        MLE_HWE_GF[2] = MLE_HWE_AF[1]*MLE_HWE_AF[1];
    }
    else
    {
        n = 0;
        int32_t imap[no_samples];
        int32_t no_genotypes = bcf_an2gn(no_alleles);
        for (size_t k=0; k<no_samples; ++k)
        {
            if (pls[k*no_genotypes]!=bcf_int32_missing)
            {
                imap[n] = k;
                ++n;
            }
        }

        if (!n) return;

        float af[no_alleles];
        float p = 1.0/no_alleles;
        for (size_t i=0; i<no_alleles; ++i)
        {
            af[i] = p;
        }
        float gf[no_genotypes];
        float gf_indiv[no_genotypes];

        bool debug = false;

        float mse = e+1;
        while (mse>e && iter<50)
        {
            for (size_t i=0; i<no_alleles; ++i)
            {
                MLE_HWE_AF[i] = 0;
                for (size_t j=0; j<=i; ++j)
                {
                    gf[bcf_alleles2gt(i,j)] = (i!=j?2:1)*af[i]*af[j];
                }
            }

            //iterate through individuals
            for (size_t k=0; k<n; ++k)
            {
                size_t offset = imap[k]*no_genotypes;

                float prob_data = 0;
                for (size_t i=0; i<no_genotypes; ++i)
                {
                    prob_data += (gf_indiv[i] = gf[i]*lt->pl2prob(pls[offset+i]));
                }

                for (size_t i=0; i<no_genotypes; ++i)
                {
                    gf_indiv[i] /= prob_data;
                }

                for (size_t i=0; i<no_alleles; ++i)
                {
                    for (size_t j=0; j<=i; ++j)
                    {
                        int32_t gf_index = bcf_alleles2gt(i,j);
                        MLE_HWE_AF[i] += 0.5*gf_indiv[gf_index];
                        MLE_HWE_AF[j] += 0.5*gf_indiv[gf_index];
                    }
                }
            }

            //normalize to frequency
            mse = 0;
            float diff;
            for (size_t i=0; i<no_alleles; ++i)
            {
                MLE_HWE_AF[i] /= n;
                diff = af[i]-MLE_HWE_AF[i];
                mse += (diff*diff);
                af[i] = MLE_HWE_AF[i];
            }

            ++iter;
        }

        for (size_t i=0; i<no_alleles; ++i)
        {
            for (size_t j=0; j<=i; ++j)
            {
                MLE_HWE_GF[bcf_alleles2gt(i,j)] = (i!=j?2:1)*MLE_HWE_AF[i]*MLE_HWE_AF[j];
            }
        }
    }
}

/**
 * Computes allele frequencies using EM algorithm from genotype likelihoods.
 *
 * @pls        - PHRED genotype likelihoods
 * @no_samples - number of samples
 * @ploidy     - ploidy
 * @n_alleles  - number of alleles
 * @MLE_AF     - estimated AF
 * @MLE_GF     - estimated GF
 * @n          - effective sample size
 * @e          - error
 */
void Estimator::compute_gl_af(int32_t *pls, int32_t nsamples, int32_t ploidy,
                int32_t n_allele, float *MLE_AF, float *MLE_GF, int32_t& n,
                double e)
{
    int32_t iter = 0;

    if (n_allele==2 && ploidy==2)
    {
        n = 0;
        int32_t imap[nsamples];

        for (size_t i=0; i<nsamples; ++i)
        {
            if (pls[i*3]!=bcf_int32_missing)
            {
                imap[n] = i;
                ++n;
            }
        }

        if (!n)
        {
            return;
        }

        float gf[3];
        float gf_indiv[3];

        float mse = e+1;
        float diff = 0;

        gf[0] = 1.0/3;
        gf[1] = gf[0];
        gf[2] = gf[1];

        while (mse>e && iter<50)
        {
            MLE_GF[0] = 0;
            MLE_GF[1] = 0;
            MLE_GF[2] = 0;

            for (size_t i=0; i<n; ++i)
            {
                size_t offset = imap[i]*3;

                float prob_data = (gf_indiv[0] = gf[0]*lt->pl2prob(pls[offset]));
                prob_data += (gf_indiv[1] = gf[1]*lt->pl2prob(pls[offset+1]));
                prob_data += (gf_indiv[2] = gf[2]*lt->pl2prob(pls[offset+2]));

                MLE_GF[0] += (gf_indiv[0] /= prob_data);
                MLE_GF[1] += (gf_indiv[1] /= prob_data);
                MLE_GF[2] += (gf_indiv[2] /= prob_data);
            }

            MLE_GF[0] /= n;
            MLE_GF[1] /= n;
            MLE_GF[2] /= n;

            diff = (gf[0]-MLE_GF[0]);
            mse = diff*diff;
            diff = (gf[1]-MLE_GF[1]);
            mse += diff*diff;
            diff = (gf[2]-MLE_GF[2]);
            mse += diff*diff;

            gf[0] = MLE_GF[0];
            gf[1] = MLE_GF[1];
            gf[2] = MLE_GF[2];


            ++iter;
        }

        MLE_AF[0] = MLE_GF[0]+0.5*MLE_GF[1];
        MLE_AF[1] = MLE_GF[2]+0.5*MLE_GF[1];
    }
    else
    {
        n = 0;
        int32_t imap[nsamples];
        int32_t no_genotypes = bcf_an2gn(n_allele);
        for (size_t i=0; i<nsamples; ++i)
        {
            if (pls[i*no_genotypes]!=bcf_int32_missing)
            {
                imap[n] = i;
                ++n;
            }
        }

        if (!n) return;

        float gf[no_genotypes];
        float gf_indiv[no_genotypes];

        //initialization
        gf[0] = 1.0/no_genotypes;
        for (size_t i=1; i<no_genotypes; ++i)
        {
            gf[i] = gf[0];
        }

        float mse = e+1;
        while (mse>e && iter<50)
        {
            //initialization
            for (size_t i=0; i<no_genotypes; ++i)
            {
                MLE_GF[i] = 0;
            }

            //iterate through individuals
            for (size_t i=0; i<n; ++i)
            {
                size_t offset = imap[i]*no_genotypes;

                float prob_data = 0;
                for (size_t j=0; j<no_genotypes; ++j)
                {
                    prob_data += (gf_indiv[j] = gf[j]*lt->pl2prob(pls[offset+j]));
                }

                for (size_t j=0; j<no_genotypes; ++j)
                {
                    MLE_GF[j] += (gf_indiv[j] /= prob_data);
                }
            }

            mse = 0;
            float diff;
            for (size_t i=0; i<no_genotypes; ++i)
            {
                MLE_GF[i] /= n;
                diff = gf[i]-MLE_GF[i];
                mse += (diff *= diff);
                gf[i] = MLE_GF[i];
            }

            ++iter;
        }

        for (size_t i=0; i<n_allele; ++i)
        {
            for (size_t j=0; j<=i; ++j)
            {
                int32_t index = bcf_alleles2gt(i,j);
                MLE_AF[i] += 0.5*MLE_GF[index];
                MLE_AF[i] += 0.5*MLE_GF[index];
            }
        }
    }
}

/**
 * Computes the Hardy-Weinberg Likelihood Ratio Test Statistic
 *
 * @pls        - PHRED genotype likelihoods
 * @no_samples - number of samples
 * @ploidy     - ploidy
 * @no_alleles - number of alleles
 * @MLE_HWE_AF - estimated AF assuming HWE
 * @MLE_GF     - estimated GF
 * @n          - effective sample size
 * @lr         - log10 likelihood ratio
 * @logp       - likelihood ratio test log p-value
 * @df         - degrees of freedom
 *
 */
void Estimator::compute_hwe_lrt(int32_t *pls, int32_t no_samples, int32_t ploidy,
            int32_t no_alleles, float *MLE_HWE_GF, float *MLE_GF, int32_t& n,
            float& lr, float& logp, int32_t& df)
{
    n = 0;
    if (ploidy==2)
    {
        if (no_alleles==2)
        {
            int32_t no_genotypes = 3;

            float l0=0, la=0;
            for (size_t k=0; k<no_samples; ++k)
            {
                size_t offset = k*3;

                if (pls[offset]==bcf_int32_missing) continue;

                ++n;

                float p = lt->pl2prob(pls[offset]);
                float l0i = MLE_HWE_GF[0] * p;
                float lai = MLE_GF[0] * p;
                p = lt->pl2prob(pls[offset+1]);
                l0i += MLE_HWE_GF[1] * p;
                lai += MLE_GF[1] * p;
                p = lt->pl2prob(pls[offset+2]);
                l0i += MLE_HWE_GF[2] * p;
                lai += MLE_GF[2] * p;

                l0 += log(l0i);
                la += log(lai);
            }

            if (!n) return;

            lr = l0-la;
            float lrts = lr>0 ? 0 : -2*lr;
            df = 1;
            logp = pchisq(lrts, 1, 0, 1);
        }
        else
        {
            int32_t no_genotypes = bcf_an2gn(no_alleles);

            float l0=0, la=0;
            for (size_t k=0; k<no_samples; ++k)
            {
                size_t offset = k*no_genotypes;

                if (pls[offset]==bcf_int32_missing) continue;

                ++n;

                float l0i=0, lai=0;
                for (size_t j=0; j<no_genotypes; ++j)
                {
                    float p = lt->pl2prob(pls[offset+j]);
                    l0i += MLE_HWE_GF[j]*p;
                    lai += MLE_GF[j]*p;
                }

                l0 += log(l0i);
                la += log(lai);
            }

            if (!n) return;

            lr = l0-la;
            float lrts = lr>0 ? 0 : -2*lr;
            df = no_genotypes-no_alleles;
            logp = pchisq(lrts, df, 0, 1);
        }
    }
};

/**
 * Computes the Inbreeding Coefficient Statistic from Genotype likelihoods.
 *
 * @pls        - PHRED genotype likelihoods
 * @no_samples - number of samples
 * @ploidy     - ploidy
 * @GF         - GF
 * @HWE_AF     - AF under HWE assumption
 * @no_alleles - number of alleles
 * @F          - estimated inbreeding coefficient
 * @n          - effective sample size
 */
void Estimator::compute_gl_fic(int32_t * pls, int32_t no_samples, int32_t ploidy,
                               float* HWE_AF, int32_t no_alleles, float* GF,
                               float& F, int32_t& n)
{
    n = 0;

    if (ploidy!=2)
    {
        return;
    }

    float HWE_GF[3];
    HWE_GF[0] = HWE_AF[0]*HWE_AF[0];
    HWE_GF[1] = 2*HWE_AF[0]*HWE_AF[1];
    HWE_GF[2] = HWE_AF[1]*HWE_AF[1];

    if (no_alleles==2)
    {
        int32_t no_genotypes = 3;

        float num = 0, denum=0;
        for (size_t k=0; k<no_samples; ++k)
        {
            size_t offset = k*no_genotypes;

            if (pls[offset]==bcf_int32_missing)
            {
                continue;
            }

            ++n;

            float o_het_sum = lt->pl2prob(pls[offset+1])*GF[1];
            float o_sum = lt->pl2prob(pls[offset])*GF[0];
            o_sum += o_het_sum;
            o_sum += lt->pl2prob(pls[offset+2])*GF[2];

            float e_het_sum = lt->pl2prob(pls[offset+1])*HWE_GF[1];
            float e_sum = lt->pl2prob(pls[offset])*HWE_GF[0];
            e_sum += e_het_sum;
            e_sum += lt->pl2prob(pls[offset+2])*HWE_GF[2];

            num += o_het_sum/o_sum;
            denum += e_het_sum/e_sum;
        }

        F = 1-num/denum;
    }
    else
    {
        int32_t no_genotypes = bcf_an2gn(no_alleles);

        float HWE_GF[no_genotypes];

        for (size_t i=0; i<no_alleles; ++i)
        {
            for (size_t j=0; j<=i; ++j)
            {
                HWE_GF[bcf_alleles2gt(i,j)] = (i!=j?2:1)*HWE_AF[i]*HWE_AF[j];
            }
        }

        float num=0, denum=0;
        float o_het_sum;
        float o_sum;
        float e_het_sum;
        float e_sum;
        for (size_t k=0; k<no_samples; ++k)
        {
            size_t offset = k*no_genotypes;
            if (pls[offset]==bcf_int32_missing)
            {
                continue;
            }

            ++n;

            o_het_sum = 0;
            o_sum = 0;
            e_het_sum = 0;
            e_sum = 0;
            int32_t gt_index = 0;
            for (size_t i=0; i<no_alleles; ++i)
            {
                for (size_t j=0; j<i; ++j)
                {
                    float p = lt->pl2prob(pls[offset+gt_index]);
                    o_het_sum += p * GF[gt_index];
                    o_sum += p * GF[gt_index];

                    e_het_sum += p * HWE_GF[gt_index];
                    e_sum += p * HWE_GF[gt_index];

                    ++gt_index;
                }

                //for homozygote
                o_sum += lt->pl2prob(pls[offset+gt_index]) * GF[gt_index];
                e_sum += lt->pl2prob(pls[offset+gt_index]) * HWE_GF[gt_index];
                ++gt_index;
            }

            num += o_het_sum/o_sum;
            denum += e_het_sum/e_sum;
        }

        F = 1-num/denum;
    }
};

/**
 * Computes Allele Balance from genotype likelihoods.
 *
 * @pls        - PHRED genotype likelihoods
 * @no_samples - number of samples
 * @ploidy     - ploidy
 * @dps        - depths
 * @GF         - estimated GF
 * @no_alleles - number of alleles
 * @ab         - estimate of allele balance
 * @n          - effective sample size
 */
void Estimator::compute_gl_ab(int32_t *pls, int32_t no_samples, int32_t ploidy,
                              int32_t *dps,
                              float* GF, int32_t no_alleles,
                              float& ab, int32_t& n)
{
    n = 0;

    if (ploidy!=2) return;

    if (no_alleles==2)
    {
        float num = 0, denum = 0;
        for (size_t k=0; k<no_samples; ++k)
        {
            size_t offset = k*3;
            if(pls[offset]!=bcf_int32_missing && dps[k]!=0)
            {
                float nrefnum = pls[offset+2]-pls[offset+0];
                float nrefdenum = pls[offset]+pls[offset+2]-2*pls[offset+1] +6*dps[k];
                float nref = 0.5*dps[k]*(1+(nrefdenum?nrefnum/nrefdenum:0));
                float phet = lt->pl2prob(pls[offset+1])*GF[1] /
                             ( lt->pl2prob(pls[offset])*GF[0]
                              +lt->pl2prob(pls[offset+1])*GF[1]
                              +lt->pl2prob(pls[offset+2])*GF[2]);

                num += phet*nref;
                denum += phet*dps[k];
                ++n;
            }
        }

        ab = (0.05+num)/(0.10+denum);
    }
    else
    {
        int32_t no_genotypes = bcf_an2gn(no_alleles);
        float num = 0, denum = 0;
        for (size_t k=0; k<no_samples; ++k)
        {
            size_t offset = k*no_genotypes;
            if(pls[offset]!=bcf_int32_missing)
            {
                float prob_data, p_ref;
                int32_t gt_index = 0;

                for (size_t j=1; j<no_alleles; ++j)
                {
                    size_t het_index = bcf_alleles2gt(0,j);
                    size_t homalt_index = bcf_alleles2gt(j,j);
                    float nrefnum = pls[offset+homalt_index]-pls[offset];
                    float nrefdenum = pls[offset]+pls[offset+homalt_index]-2*pls[offset+het_index] +6*dps[k];
                    float nref = 0.5*dps[k]*(1+(nrefdenum?nrefnum/nrefdenum:0));

                    float n = lt->pl2prob(pls[offset+het_index])*GF[het_index] ;
                    float d = (lt->pl2prob(pls[offset])*GF[0]
                               +n
                               +lt->pl2prob(pls[offset+homalt_index])*GF[homalt_index]);
                    float phet = d?n/d:0.333;
                    num += phet*nref;
                    denum += phet*dps[k];
                }

                ++n;
            }
        }

        ab = (0.05+num)/(0.10+denum);
    }
};

/**
 * Computes the phred scaled QUAL for a variant.
 *
 * @pls        - PHRED genotype likelihoods
 * @no_samples - number of samples
 * @ploidy     - ploidy
 * @no_alleles - number of alleles
 * @n          - effective sample size
 * @qual       - PHRED scaled QUAL
 */
void Estimator::compute_qual(int32_t *pls, int32_t no_samples, int32_t ploidy,
            int32_t no_alleles, float &qual, int32_t &n)
{
    n = 0;
    if (ploidy==2)
    {
        if (no_alleles==2)
        {
            int32_t no_genotypes = 3;

            qual = 0;
            for (size_t k=0; k<no_samples; ++k)
            {
                size_t offset = k*3;

                if (pls[offset]==bcf_int32_missing) continue;

                ++n;

                qual += lt->log10((1-lt->pl2prob(pls[offset])/(lt->pl2prob(pls[offset])+lt->pl2prob(pls[offset+1])+lt->pl2prob(pls[offset+2]))));

            }

            if (!n) return;

            qual = lt->round(-qual*10);
        }
        else
        {
            //works only for ploidy of 2

            int32_t no_genotypes = bcf_an2gn(no_alleles);

            float gq = 0;
            for (size_t k=0; k<no_samples; ++k)
            {
                size_t offset = k*no_genotypes;

                if (pls[offset]==bcf_int32_missing) continue;

                ++n;

                float denom = 0;
                for (size_t j=0; j<no_genotypes; ++j)
                {
                    denom += lt->pl2prob(pls[offset]);
                }

                qual += lt->log10((1-lt->pl2prob(pls[offset])/denom));
            }

            if (!n) return;

            qual = lt->round(-qual*10);
        }
    }
};