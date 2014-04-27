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
 * Computes allele frequencies using EM algorithm from genotype likelihoods 
 * under assumption of Hardy-Weinberg Equilibrium.
 *
 * @pls        - PHRED genotype likelihoods
 * @nsamples   - number of samples
 * @ploidy     - ploidy
 * @n_alleles  - number of alleles
 * @MLE_HWE_AF - estimated AF
 * @MLE_HWE_GF - estimated GF
 * @n          - effective sample size
 * @e          - error
 */
void Estimator::compute_gl_af_hwe(int32_t *pls, int32_t nsamples, int32_t ploidy,
                int32_t n_allele, float *MLE_HWE_AF, float *MLE_HWE_GF, int32_t& n, 
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

        float af[2] = {0.5, 0.5};
        float gf[3];
        float gf_indiv[3];

        float mse = e+1;
        float diff = 0;
        while (mse>e)
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
        int32_t imap[nsamples];
        int32_t n_genotype = bcf_an2gn(n_allele);
        for (size_t i=0; i<nsamples; ++i)
        {
            if (pls[i*n_genotype]!=bcf_int32_missing)
            {
                imap[n] = i;
                ++n;
            }
        }

        if (!n) return;

        float af[n_allele];
        float p = 1.0/n_allele;
        for (size_t i=0; i<n_allele; ++i)
        {
            af[i] = p;
        }
        float gf[n_genotype];
        float gf_indiv[n_genotype];

        bool debug = false;

        float mse = e+1;
        while (mse>e && iter<50)
        {
            //initialization
            for (size_t i=0; i<n_genotype; ++i)
            {
                gf[i] = 0;
            }

            for (size_t i=0; i<n_allele; ++i)
            {
                MLE_HWE_AF[i] = 0;
                for (size_t j=0; j<=i; ++j)
                {
                    gf[bcf_alleles2gt(i,j)] += (i!=j?2:1)*af[i]*af[j];
                }
            }

            //iterate through individuals
            for (size_t i=0; i<n; ++i)
            {
                size_t offset = imap[i]*n_genotype;

                float prob_data = 0;
                for (size_t j=0; j<n_genotype; ++j)
                {
                    prob_data += (gf_indiv[j] = gf[j]*lt->pl2prob(pls[offset+j]));
                }

                for (size_t j=0; j<n_genotype; ++j)
                {
                    gf_indiv[j] /= prob_data;
                }

                for (size_t j=0; j<n_allele; ++j)
                {
                    for (size_t k=0; k<j; ++k)
                    {
                        int32_t gf_index = bcf_alleles2gt(j,k);
                        MLE_HWE_AF[j] += 0.5*gf_indiv[gf_index];
                        MLE_HWE_AF[k] += 0.5*gf_indiv[gf_index];
                    }
                    MLE_HWE_AF[j] += gf_indiv[bcf_alleles2gt(j,j)];
                }
            }

            //normalize to frequency
            mse = 0;
            float diff;
            for (size_t i=0; i<n_allele; ++i)
            {
                MLE_HWE_AF[i] /= n;
                diff = af[i]-MLE_HWE_AF[i];
                mse += (diff *= diff);
                af[i] = MLE_HWE_AF[i];
            }

            ++iter;
        }

        for (size_t i=0; i<n_allele; ++i)
        {
            MLE_HWE_GF[i] = 0;
        }

        for (size_t i=0; i<n_allele; ++i)
        {
            for (size_t j=0; j<=i; ++j)
            {
                MLE_HWE_GF[bcf_alleles2gt(i,j)] += (i!=j?2:1)*MLE_HWE_AF[i]*MLE_HWE_AF[j];
            }
        }
    }
}

/**
 * Computes allele frequencies using EM algorithm from genotype likelihoods.
 *
 * @pls        - PHRED genotype likelihoods
 * @nsamples   - number of samples
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

                MLE_GF[0] = (gf_indiv[0] /= prob_data);
                MLE_GF[1] = (gf_indiv[1] /= prob_data);
                MLE_GF[2] = (gf_indiv[2] /= prob_data);
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
        int32_t n_genotype = bcf_an2gn(n_allele);
        for (size_t i=0; i<nsamples; ++i)
        {
            if (pls[i*n_genotype]!=bcf_int32_missing)
            {
                imap[n] = i;
                ++n;
            }
        }

        if (!n) return;

        float gf[n_genotype];
        float gf_indiv[n_genotype];

        float mse = e+1;
        while (mse>e && iter<50)
        {
            //initialization
            for (size_t i=0; i<n_genotype; ++i)
            {
                MLE_GF[i] = 0;
            }

            //iterate through individuals
            for (size_t i=0; i<n; ++i)
            {
                size_t offset = imap[i]*n_genotype;

                float prob_data = 0;
                for (size_t j=0; j<n_genotype; ++j)
                {
                    prob_data += (gf_indiv[j] = gf[j]*lt->pl2prob(pls[offset+j]));
                }

                for (size_t j=0; j<n_genotype; ++j)
                {
                    MLE_GF[j] = (gf_indiv[j] /= prob_data);
                }
            }

            mse = 0;
            float diff;
            for (size_t i=0; i<n_allele; ++i)
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
            MLE_AF[i] = 0;
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
 * @nsamples   - number of samples
 * @ploidy     - ploidy
 * @n_allele   - number of alleles
 * @MLE_HWE_AF - estimated AF assuming HWE
 * @MLE_GF     - estimated GF
 * @n          - effective sample size
 * @lrts       - log10 likelihood ratio test statistic p value
 * @dof        - degrees of freedom
 *
 */
void Estimator::compute_hwe_lrt(int32_t *pls, int32_t nsamples, int32_t ploidy,
            int32_t n_allele, float *MLE_HWE_GF, float *MLE_GF, int32_t& n,
            float& lrts, float& logp, int32_t& df)                
{
    if (n_allele==2 && ploidy==2)
    {    
        int32_t n_genotype = 3;
        df = 1;
    
        //compute LRT statistic
        float l0=0, l0i=0, la=0, lai=0;
        int32_t n = 0;
        for (size_t i=0; i<nsamples; ++i)
        {
            size_t offset = i*3;
    
            if (pls[offset]==bcf_int32_missing) continue;
            ++n;
                
            l0i=0;
            lai=0;
            for (size_t j=0; j<n_genotype; ++j)
            {
                int32_t pl = lt->pl2prob(pls[offset+j]);
                l0i += MLE_HWE_GF[j] * pl;
                lai += MLE_GF[j] * pl;
            }
    
            l0 += log10(l0i);
            la += log10(lai);
        }
    
        if (!n) return;
        
        lrts = -2*(l0-la)<0 ? 0 : -2*(l0-la);
        logp = pchisq(lrts, df, 0, 1);
    }
};

/**
 * Computes the Inbreeding Coefficient Statistic from Genotype likelihoods.
 *
 * @pls        - PHRED genotype likelihoods
 * @nsamples   - number of samples
 * @ploidy     - ploidy
 * @n_allele   - number of alleles
 * @MLE_HWE_AF - estimated AF
 * @MLE_HWE_GF - estimated GF
 * @n          - effective sample size
 * @lrts       - log10 likelihood ratio test statistic p value
 * @dof        - degrees of freedom
 *
 */
bool Estimator::compute_gl_fic(int32_t * pls, float* GF, float* HWE_AF, float& F, int32_t n_allele)
{
    uint32_t noGenotypes = noAlleles*(noAlleles+1)/2;

    if (pG.size()!=noGenotypes)
    {
        return false;
    }

    //==========================================
    //sets up map from genotype code to genotype
    //==========================================
    std::vector<std::vector<uint32_t> > index2Genotype(noGenotypes);
    int32_t genotypeCode = 0;
    for (uint32_t i=0; i<noAlleles; ++i)
    {
        for (uint32_t j=0; j<=i; ++j)
        {
            index2Genotype[genotypeCode].resize(2);
            index2Genotype[genotypeCode][0] = j;
            index2Genotype[genotypeCode][1] = i;

            ++genotypeCode;
        }
    }

    //compute genotype frequencies under HWE
    std::vector<double> pG_HWE(noGenotypes);
    for (uint32_t j=0; j<noGenotypes; ++j)
    {
        pG_HWE[j] = pA_HWE[index2Genotype[j][0]]*pA_HWE[index2Genotype[j][1]];

        if(index2Genotype[j][0]!=index2Genotype[j][1])
        {
            pG_HWE[j] *= 2;
        }
    }

    //compute Fst
    double FPNum=0, FPDenum=0;
    double pG_Reads=0;
    double pHet_Reads_sum=0;
    double pG_Reads_sum=0;
    double pHet_sum=0;

    double pG_ReadsHWE=0;
    double pHet_ReadsHWE_sum=0;
    double pG_ReadsHWE_sum=0;
    double pHet_HWE_sum=0;

    for (uint32_t i=0; i<GLs.size(); ++i)
    {
        //missing data
        if (GLs[i].size()==0)
        {
            continue;
        }

        pHet_Reads_sum=0;
        pG_Reads_sum=0;
        pHet_sum=0;

        pHet_ReadsHWE_sum=0;
        pG_ReadsHWE_sum=0;
        pHet_HWE_sum=0;

        for (uint32_t j=0; j<noGenotypes; ++j)
        {
            pG_Reads = GLs[i][j] * pG[j];
            pG_ReadsHWE = GLs[i][j] * pG_HWE[j];

            //hets
            if(index2Genotype[j][0]!=index2Genotype[j][1])
            {
                pHet_Reads_sum += pG_Reads;
                pHet_sum += pG[j];

                pHet_ReadsHWE_sum += pG_ReadsHWE;
                pHet_HWE_sum += pG_HWE[j];
            }

            pG_Reads_sum += pG_Reads;
            pG_ReadsHWE_sum += pG_ReadsHWE;
        }

        FPNum += pHet_ReadsHWE_sum/pG_ReadsHWE_sum;
        FPDenum += pHet_HWE_sum;
    }

    F = 1-FPNum/FPDenum;
    return true;
};

/**
 * Computes the Allele Balance from genotype likelihoods.
 *
 * @pls        - PHRED genotype likelihoods
 * @nsamples   - number of samples
 * @ploidy     - ploidy
 * @n_allele   - number of alleles
 * @MLE_HWE_AF - estimated AF
 * @MLE_HWE_GF - estimated GF
 * @n          - effective sample size
 * @lrts       - log10 likelihood ratio test statistic p value
 * @dof        - degrees of freedom
 *
 */
void estimateAlleleBalance(std::vector<std::vector<uint32_t> >& PLs, std::vector<std::vector<double> >& GLs, std::vector<uint32_t>& DPs, std::vector<double>& genotypeFreq, double& ab)
{
    double num = 0, denum = 0;
    uint32_t n=0;
    for (uint32_t i=0; i<GLs.size(); ++i)
    {
        if(GLs[i].size()!=0 && DPs[i]!=0)
        {
            double nrefnum = (double)PLs[i][2]-(double)PLs[i][0];
            double nrefdenum = (double)PLs[i][0]+(double)PLs[i][2]-2*(double)PLs[i][1] +6*DPs[i];
            double nref = 0.5*DPs[i]*(1+nrefnum/nrefdenum);
            double phet = GLs[i][1]*genotypeFreq[1] / (GLs[i][0]*genotypeFreq[0]+GLs[i][1]*genotypeFreq[1]+GLs[i][2]*genotypeFreq[2]);

//          std::cerr << "nrefnum " << nrefnum << "\n";
//          std::cerr << "nrefdenum " << nrefdenum << "\n";
//          std::cerr << "phet " << phet << "\n";
            //std::cerr << "GLs " << GLs[i][0] << "," << GLs[i][1] << "," << GLs[i][2] << "\n";

            //double phet = genotypeFreq[1];
            num += phet*nref;
            denum += phet*DPs[i];
            ++n;
        }
    }

    //std::cerr << "num/denum " << num  << ", "  << denum << "\n";
    ab = (0.05+num)/(0.10+denum);
};
