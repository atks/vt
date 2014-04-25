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
 * Computes HWE allele frequencies using EM algorithm
 *
 * @gts        - genotypes
 * @pls        - PHRED genotype likelihoods
 * @nsamples   - number of samples
 * @ploidy     - ploidy
 * @n_alleles  - number of alleles
 * @MLE_HWE_AF - estimated AF
 * @MLE_HWE_GF - estimated GF
 * @n          - effective sample size
 * @eps        - error
 */
void Estimator::compute_hwe_af(int32_t *gts, int32_t *pls, int32_t nsamples, int32_t ploidy,
                int32_t n_allele, float *MLE_HWE_AF, float *MLE_HWE_GF, int32_t& n, double e)
{
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

        int32_t iter = 0;
        float mse = e+1;
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

            mse = (af[0]-MLE_HWE_AF[0]);
            mse *= mse;

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

        if (!n)
        {
            return;
        }

//        std::cerr << "n: " << n  << "\n";

        float af[n_allele];
        float p = 1.0/n_allele;
        for (size_t i=0; i<n_allele; ++i)
        {
            af[i] = p;
        }
        float gf[n_genotype];
        float gf_indiv[n_genotype];

        int32_t iter = 0;
        float mse = e+1;
        while (mse>e)
        {
            //initialization
            for (size_t i=0; i<n_allele; ++i)
            {
                MLE_HWE_AF[i] = 0;
                for (size_t j=0; j<=i; ++j)
                {
                    gf[bcf_alleles2gt(i,j)] = (i==j?2:1)*af[i]*af[j];
                }
            }

            //iterate through individuals
            for (size_t i=0; i<n; ++i)
            {
                size_t offset = imap[i]*n_genotype;

                std::cerr << "Priot GFs : "; 
                for (size_t i=0; i<n_genotype; ++i)
                {
                    std::cerr << " " << gf[i] ;
                }
                std::cerr << "\n";
                    

                for (size_t i=0; i<n_allele; ++i)
                {
                    std::cerr << " " << MLE_HWE_AF[i] ;
                }
                std::cerr << "\n";

                float prob_data = 0;
                for (size_t i=0; i<n_genotype; ++i)
                {
                    prob_data += (gf_indiv[i] = gf[i]*lt->pl2prob(pls[offset+i]));
                }
                
                 
//                std::cerr << imap[i] << ") PSEUDO COUNTS " << gf_indiv[0] << " " << gf_indiv[1] << " " << gf_indiv[2]
//                    <<  " prob_data : " << prob_data
//                    << pls[offset] << " " << pls[offset+1] << " " << pls[offset+2] <<"\n";

                for (size_t i=0; i<n_genotype; ++i)
                {
                    gf_indiv[i] /= prob_data;
                }
//                  
                std::cerr << "Individual GLs : "; 
                for (size_t i=0; i<n_genotype; ++i)
                {
                    std::cerr << " " << gf_indiv[i] ;
                }
                std::cerr << "\n";

//                for (size_t i=0; i<n_allele; ++i)
//                {
//                    MLE_HWE_AF[i] = 0;
//                }
                
                for (size_t i=0; i<n_allele; ++i)
                {
                    for (size_t j=0; j<i; ++j)
                    {
                        int32_t gf_index = bcf_alleles2gt(i,j);
                        MLE_HWE_AF[i] += 0.5*gf[gf_index];
                        MLE_HWE_AF[j] += 0.5*gf[gf_index];
                    }
                    MLE_HWE_AF[i] += gf[bcf_alleles2gt(i,i)];
                }

//                std::cerr << "\tMLE COUNTS: " << MLE_HWE_AF[0] << ":" << MLE_HWE_AF[1] << " " << n << "\n";
//
                if (i==10)exit(1);
            }

//            std::cerr << "MLE_PSEUDO_COUNTS : " << MLE_HWE_AF[0] << ":" << MLE_HWE_AF[0] << "\n";

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
            
//            std::cerr << af[0] << ":" << af[1] << "\n";

            ++iter;
        }

        for (size_t i=0; i<n_allele; ++i)
        {
            for (size_t j=0; j<=i; ++j)
            {
                MLE_HWE_GF[bcf_alleles2gt(i,j)] = (i==j?2:1)*MLE_HWE_AF[i]*MLE_HWE_AF[j];
            }
        }
    }
}


/**
Computes HWE allele frequencies using EM algorithm
Input:
1)Genotype Likelihoods

Output:
1) estimated HWE allele frequencies
2) sample size
*/
bool estimateHWEAlleleFrequencies(std::vector<std::vector<double> >& GLs, double eps, std::vector<double>& mleHWEAlleleFreq, uint32_t& N, uint32_t noAlleles)
{
    if (noAlleles==2)
    {
        N = 0;
        for (uint32_t i=0; i<GLs.size(); ++i)
        {
            //count non missing data
            if (GLs[i].size()!=0)
            {
                ++N;
            }
        }

        if (N==0)
        {
            return false;
        }

        std::vector<double> genotypeFreqPrior(3);
        std::vector<double> alleleFreqPrior(noAlleles, 0.5);
        std::vector<double> individualGenotypePosterior(3);
        mleHWEAlleleFreq.resize(2);
        double probData = 0;

        int iter = 0;
        double mse = eps+1;
        while (mse>eps)
        {
            genotypeFreqPrior[0] = alleleFreqPrior[0]*alleleFreqPrior[0];
            genotypeFreqPrior[1] = 2*alleleFreqPrior[0]*alleleFreqPrior[1];
            genotypeFreqPrior[2] = alleleFreqPrior[1]*alleleFreqPrior[1];

            mleHWEAlleleFreq[0] = 0;
            mleHWEAlleleFreq[1] = 0;

            for (uint32_t i=0; i<GLs.size(); ++i)
            {
                if (GLs[i].size()==0)
                {
                    continue;
                }

                //std::cerr << "GLS " << GLs[i][0] << " " << GLs[i][0] << " "  << GLs[i][0] << "\n" ;

                individualGenotypePosterior[0] = genotypeFreqPrior[0]*GLs[i][0];
                individualGenotypePosterior[1] = genotypeFreqPrior[1]*GLs[i][1];
                individualGenotypePosterior[2] = genotypeFreqPrior[2]*GLs[i][2];
                probData = individualGenotypePosterior[0];
                probData += individualGenotypePosterior[1];
                probData += individualGenotypePosterior[2];

                //std::cerr << "tot " << probData <<  "\n" ;

                individualGenotypePosterior[0] /= probData;
                individualGenotypePosterior[1] /= probData;
                individualGenotypePosterior[2] /= probData;

                mleHWEAlleleFreq[0] += individualGenotypePosterior[0] + 0.5*individualGenotypePosterior[1];
                mleHWEAlleleFreq[1] += individualGenotypePosterior[2] + 0.5*individualGenotypePosterior[1];

                //std::cerr << "AF " << mleHWEAlleleFreq[0] <<  " " << mleHWEAlleleFreq[1] <<"\n" ;
            }

            mleHWEAlleleFreq[0] /= N;
            mleHWEAlleleFreq[1] /= N;

            mse = (alleleFreqPrior[0]-mleHWEAlleleFreq[0])*(alleleFreqPrior[0]-mleHWEAlleleFreq[0]);
            mse += (alleleFreqPrior[1]-mleHWEAlleleFreq[1])*(alleleFreqPrior[1]-mleHWEAlleleFreq[1]);

            alleleFreqPrior[0] = mleHWEAlleleFreq[0];
            alleleFreqPrior[1] = mleHWEAlleleFreq[1];

            ++iter;
        }

        return true;
    }
    else
    {
        //effective size
        N = 0;
        uint32_t noGenotypes = (noAlleles*(noAlleles+1))/2;

        for (uint32_t i=0; i<GLs.size(); ++i)
        {
            //count non missing data
            if (GLs[i].size()!=0)
            {
                ++N;
            }
        }

        if (N==0)
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

        std::vector<double> genotypeFreqPrior(noGenotypes);
        std::vector<double> alleleFreqPrior(noAlleles, 1.0/noAlleles);
        std::vector<double> individualGenotypePosterior(noGenotypes);
        mleHWEAlleleFreq.resize(noAlleles);
        double probData;

        //track convergence
        int iter = 0;
        //double oldRelDiff = 0;
        //double relDiff = 0;
        //mse convergence used as it is guaranteed to converge thanks to MLE consistency property
        double mse = DBL_MAX;
        //check for convergence on prior
        //if(relDiff<eps || check_tol(oldRelDiff, relDiff, eps))
        while (mse>eps)
        {
            //initialize genotype prior under HWE
            for (uint32_t j=0; j<noGenotypes; ++j)
            {
                genotypeFreqPrior[j] = alleleFreqPrior[index2Genotype[j][0]]*alleleFreqPrior[index2Genotype[j][1]];

                if(index2Genotype[j][0]!=index2Genotype[j][1])
                {
                    genotypeFreqPrior[j] *= 2;
                }
            }

            //initialize MLE
            for (uint32_t j=0; j<noAlleles; ++j)
            {
                mleHWEAlleleFreq[j] = 0;
            }

            for (uint32_t i=0; i<GLs.size(); ++i)
            {
                //missing data
                if (GLs[i].size()==0)
                {
                    continue;
                }

                //======
                //E Step
                //======
                //compute genotype posterior probabilities
                probData = 0;
                for (uint32_t j=0; j<noGenotypes; ++j)
                {
                    individualGenotypePosterior[j] = genotypeFreqPrior[j]*GLs[i][j];
                    probData += individualGenotypePosterior[j];
                }

                for (uint32_t j=0; j<noGenotypes; ++j)
                {
                    individualGenotypePosterior[j] /= probData;

                    //==============
                    //M Step Summing
                    //==============
                    mleHWEAlleleFreq[index2Genotype[j][0]] += 0.5*individualGenotypePosterior[j];
                    mleHWEAlleleFreq[index2Genotype[j][1]] += 0.5*individualGenotypePosterior[j];
                }
            }

            mse = 0;
            for (uint32_t j=0; j<noAlleles; ++j)
            {
                //===============
                //M Step Dividing
                //===============
                mleHWEAlleleFreq[j] /= N;
                mse += (alleleFreqPrior[j]-mleHWEAlleleFreq[j])*(alleleFreqPrior[j]-mleHWEAlleleFreq[j]);
                alleleFreqPrior[j] = mleHWEAlleleFreq[j];
                //oldRelDiff = relDiff;
                //relDiff += update(alleleFreqPrior[j], mleHWEAlleleFreq[j]);
            }

            /*
            print("ALLELE ");
            print(iter);
            print(" ");
            print(oldRelDiff);
            print(" ");
            print(relDiff);
            print(" ");
            println(mleHWEAlleleFreq);
            */

            ++iter;
        }

        return true;
    }
};


/**
Computes genotype frequencies using EM algorithm
Input:
1)Genotype Likelihoods (qscores)

Output:
1) estimated genotype frequencies
2) maximum likelihood (Q score)

*/
bool estimateGenotypeFrequencies(std::vector<std::vector<double> >& GLs,
                                 double eps, std::vector<double>& mleGenotypeFreq,
                                 uint32_t& N, uint32_t noAlleles)
{
    //std::cerr << "BEFORE In GL# " << GLs.size() << "\n";
    //std::cerr << "BEFORE In GF estimate allele# " << noAlleles << "\n";
    //std::cerr << "BEFORE In GF estimate genotype# " << mleGenotypeFreq.size() << "\n";

    if (noAlleles==2)
    {
        N = 0;
        for (uint32_t i=0; i<GLs.size(); ++i)
        {
            if (GLs[i].size()!=0)
            {
                ++N;
            }
        }

        if(N==0)
        {
            return false;
        }

        std::vector<double> genotypeFreqPrior(3, 1.0/3.0);
        std::vector<double> individualGenotypeFreqPosterior(3);
        mleGenotypeFreq.resize(3);
        double probData = 0;

        int iter = 0;
        double mse=eps+1;
        while (mse>eps && iter<1000)
        {
            //initialize MLE
            mleGenotypeFreq[0] = 0;
            mleGenotypeFreq[1] = 0;
            mleGenotypeFreq[2] = 0;

            for (uint32_t i=0; i<GLs.size(); ++i)
            {
                //missing data
                if (GLs[i].size()==0)
                    continue;

                individualGenotypeFreqPosterior[0] = genotypeFreqPrior[0] * GLs[i][0];
                individualGenotypeFreqPosterior[1] = genotypeFreqPrior[1] * GLs[i][1];
                individualGenotypeFreqPosterior[2] = genotypeFreqPrior[2] * GLs[i][2];

                probData = individualGenotypeFreqPosterior[0];
                probData += individualGenotypeFreqPosterior[1];
                probData += individualGenotypeFreqPosterior[2];

                mleGenotypeFreq[0] += individualGenotypeFreqPosterior[0]/probData;
                mleGenotypeFreq[1] += individualGenotypeFreqPosterior[1]/probData;
                mleGenotypeFreq[2] += individualGenotypeFreqPosterior[2]/probData;
            }

            mleGenotypeFreq[0] /= N;
            mleGenotypeFreq[1] /= N;
            mleGenotypeFreq[2] /= N;
            mse = (genotypeFreqPrior[0]-mleGenotypeFreq[0])*(genotypeFreqPrior[0]-mleGenotypeFreq[0]);
            mse += (genotypeFreqPrior[1]-mleGenotypeFreq[1])*(genotypeFreqPrior[1]-mleGenotypeFreq[1]);
            mse += (genotypeFreqPrior[2]-mleGenotypeFreq[2])*(genotypeFreqPrior[2]-mleGenotypeFreq[2]);

            genotypeFreqPrior[0] = mleGenotypeFreq[0];
            genotypeFreqPrior[1] = mleGenotypeFreq[1];
            genotypeFreqPrior[2] = mleGenotypeFreq[2];

            ++iter;
        }

        return true;
    }
    else
    {
        //effective size
        uint32_t N = 0;
        uint32_t noGenotypes = (noAlleles*(noAlleles+1))/2;

        for (uint32_t i=0; i<GLs.size(); ++i)
        {
            //count non missing data
            if (GLs[i].size()!=0)
            {
                ++N;

                //ensures same number of genotypes
                if(noGenotypes!=GLs[i].size())
                {
                    std::cerr << "Number of genotypes not consistent for individuals\n";
                    exit(0);
                }
            }
        }

        if (N==0)
        {
            return false;
            //throw(Exception("No observations to estimate frequencies"));
        }

        //==========================================
        //sets up map from genotype code to genotype
        //this is a critical part for multi allellic
        //loci
        //==========================================
//      uint32_t index2genotype[noGenotypes][2];
//      uint32_t genotypeCode = 0;
//      for (uint32_t i=0; i<noAlleles; ++i)
//      {
//          for (uint32_t j=0; j<=i; ++j)
//          {
//              index2genotype[genotypeCode][0] = j;
//              index2genotype[genotypeCode][1] = i;
//              ++genotypeCode;
//          }
//      }
        //by default start with a uniform prior
        std::vector<double> genotypeFreqPrior(noGenotypes, 1.0/noGenotypes);

        std::vector<double> individualGenotypeFreqPosterior(noGenotypes);
        mleGenotypeFreq.resize(noGenotypes);
        double probData;

        //track convergence
        //2 ways to terminate, MSE or by reldiff
        //mse is used right now as it terminates more reliably
        int iter = 0;
        //double oldRelDiff = 0;
        //double relDiff = 0;
        double mse=DBL_MAX;
        //check for convergence on prior
        //if(relDiff<eps || check_tol(oldRelDiff, relDiff, eps))
        while (mse>eps && iter<1000)
        {
            //initialize MLE
            for (uint32_t j=0; j<noGenotypes; ++j)
            {
                mleGenotypeFreq[j] = 0;
            }

            //across samples
            for (uint32_t i=0; i<GLs.size(); ++i)
            {
                //missing data
                if (GLs[i].size()==0)
                    continue;

                //======
                //E Step
                //======
                //compute genotype posterior probabilities
                probData = 0;
                for (uint32_t j=0; j<noGenotypes; ++j)
                {
                    individualGenotypeFreqPosterior[j] = genotypeFreqPrior[j] * GLs[i][j];
                    probData += individualGenotypeFreqPosterior[j];
                }

                for (uint32_t j=0; j<noGenotypes; ++j)
                {
                    individualGenotypeFreqPosterior[j] /= probData;

                    //======
                    //M Step
                    //======
                    mleGenotypeFreq[j] += individualGenotypeFreqPosterior[j];
                }
            }

            //std::cerr << "In GF estimate " << noGenotypes << "\n";

            mse=0;
            for (uint32_t j=0; j<noGenotypes; ++j)
            {
                //======
                //M Step
                //======
                mleGenotypeFreq[j] /= N;
                //oldRelDiff = relDiff;
                //relDiff += update(genotypeFreqPrior[j], mleGenotypeFreq[j]);
                mse += (genotypeFreqPrior[j]-mleGenotypeFreq[j])*(genotypeFreqPrior[j]-mleGenotypeFreq[j]);
                genotypeFreqPrior[j] = mleGenotypeFreq[j];
            }

            ++iter;
        }

        return true;
        //std::cerr << "In GF estimate " << noGenotypes << "\n";
    }

    //std::cerr << "In GF estimate allele#" << noAlleles << "\n";
    //std::cerr << "In GF estimate genotype#" << mleGenotypeFreq.size() << "\n";
};


/**
Performs the HWE Likelihood Ratio Test.
Input:
1)Diploid
2)Multi-allelic
3)Genotype Likelihoods (qscores)

Output:
1)Degrees of freedom
2)P-values
*/
bool hweLRT(std::vector<std::vector<double> >& GLs,
            std::vector<double>& mleGenotypeFreq,
            std::vector<double>& mleHWEAlleleFreq,
            double& lrts, double& pValue, uint32_t& dof, uint32_t noAlleles)
{
    uint32_t noGenotypes = noAlleles*(noAlleles+1)/2;
    dof = noAlleles*(noAlleles-1)/2;

    if (noAlleles!=mleHWEAlleleFreq.size() || noGenotypes!=mleGenotypeFreq.size())
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
    std::vector<double> mleHWEGenotypeFreq(noGenotypes);
    for (uint32_t j=0; j<noGenotypes; ++j)
    {
        mleHWEGenotypeFreq[j] = mleHWEAlleleFreq[index2Genotype[j][0]]*mleHWEAlleleFreq[index2Genotype[j][1]];

        if(index2Genotype[j][0]!=index2Genotype[j][1])
        {
            mleHWEGenotypeFreq[j] *= 2;
        }
    }


    //compute LRT statistic
    double l0=0, l0i=0, la=0, lai=0;
    uint32_t N = 0;
    for (uint32_t i=0; i<GLs.size(); ++i)
    {
        //missing data
        if (GLs[i].size()==0)
        {
            continue;
        }

        ++N;

        //std::cerr << "GL: ";
        l0i=0;
        lai=0;
        for (uint32_t j=0; j<noGenotypes; ++j)
        {
            //std::cerr << j << ")" << GLs[i][j] << "\t";
            //std::cerr << i  << ": " <<  GLs[i][j] << "\n";
            l0i += mleHWEGenotypeFreq[j] * GLs[i][j];
            lai += mleGenotypeFreq[j] * GLs[i][j];
            //std::cerr << i  << ": " <<  lai << "\n";
        }


        //std::cerr << "\n";
        l0 += log(l0i);
        la += log(lai);
    }
    lrts = -2*(l0-la);

    if (N==0)
    {
        return false;
    }


    if (lrts<0)
    {
        //std::cerr << "lrts less than 0 " << lrts << "\t";
        lrts = 0;
    }

 //   boost::math::chi_squared chisqDist(dof);
  //  pValue = boost::math::cdf(complement(chisqDist, lrts));

    return true;
};

/**
 Estimate FIC.
Input:
1)Diploid
2)Multi-allelic
3)Genotype Likelihoods (qscores)

Output:
1)Fis
*/
bool estimateFIC(std::vector<std::vector<double> >& GLs, //genotype likelihoods
                 std::vector<double>& pG, //prior genotype frequencies
                 std::vector<double>& pA_HWE,
                 double& F, uint32_t noAlleles)
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
    //std::cerr << "no alleles in FIC estimation " << noAlleles << " " << noGenotypes << "\n";
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
Allele Balance Statistic developed by Tom Blackwell.
Works only for biallelic variants
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

/**
Ratio of observe variance against expected variance - RSQ.
Works only for biallelic variants
*/
bool computeRSQ(std::vector<std::vector<double> >& GLs, std::vector<double>& mleHWEAlleleFreq, double& RSQ, uint32_t noAlleles)
{
    if (noAlleles!=mleHWEAlleleFreq.size())
    {
        return false;
    }

    if (noAlleles==2)
    {
        double pRR = mleHWEAlleleFreq[0]*mleHWEAlleleFreq[0];
        double pRA = 2*mleHWEAlleleFreq[0]*mleHWEAlleleFreq[1];
        double pAA = mleHWEAlleleFreq[1]*mleHWEAlleleFreq[1];

        double postRR = 0;
        double postRA = 0;
        double postAA = 0;
        double total = 0;

        double sumD2 = 0;
        double sumD = 0;
        double n = 0;
        double D_i = 0;

        for (uint32_t i=0; i<GLs.size(); ++i)
        {
            if(GLs[i].size()==0)
            {
                continue;
            }

            postRR = GLs[i][0] * pRR;
            postRA = GLs[i][1] * pRA;
            postAA = GLs[i][2] * pAA;
            total = postRR + postRA + postAA;
            postRR /= total;
            postRA /= total;
            postAA /= total;

            D_i = postRA + 2*postAA;
            sumD2 += D_i * D_i;
            sumD += D_i;

            ++n;
        }

        if (n==0)
        {
            return false;
        }

        double meanD = sumD/n;

        RSQ = ((sumD2-n*meanD*meanD)/(n-1)) / (2*mleHWEAlleleFreq[0]*mleHWEAlleleFreq[1]);
    }
    else
    {
        RSQ = 1;
    }

    return true;
};


/**
convert PLs to probabilities.
*/
double logFact(uint32_t n, std::vector<double>& LOGFACTS)
{
    if(LOGFACTS.size()==0)
    {
        LOGFACTS.push_back(0);
    }

    if (n>=LOGFACTS.size())
    {
        for (uint32_t i=LOGFACTS.size(); i<=n; ++i)
        {
            LOGFACTS.push_back(LOGFACTS[i-1] + log(i));
        }
    }

    return LOGFACTS[n];
};

double logHypergeometricProb(std::vector<double>& logFacs, uint32_t a, uint32_t b, uint32_t c, uint32_t d, std::vector<double>& LOGFACTS)
{
    return logFact(a+b, LOGFACTS) + logFact(c+d, LOGFACTS) + logFact(a+c, LOGFACTS) + logFact(b+d, LOGFACTS) - logFact(a, LOGFACTS) - logFact(b, LOGFACTS) - logFact(c, LOGFACTS) - logFact(d, LOGFACTS) - logFact(a+b+c+d, LOGFACTS);
};

double fisher1(uint32_t a, uint32_t b, uint32_t c, uint32_t d, std::vector<double>& LOGFACTS)
{
    uint32_t n = a + b + c + d;

    std::vector<double> logFacs(n+1);

    double logpCutoff = logHypergeometricProb(logFacs, a, b, c, d, LOGFACTS);
    double p2SidedFraction = 0;

    for(uint32_t x=0; x<=n; ++x)
    {
        if ( a+b >= x && a+c >= x && d+x >=a )
        {
            ///std::cerr << "ABCDX\t" << a  << "\t" << b  << "\t" << c << "\t" << d << "\t" << x << "\n";
            double l = logHypergeometricProb(logFacs, x, a+b-x, a+c-x, d-a+x, LOGFACTS);
            //std::cerr << "\t" << l  << "\t" << p2SidedFraction  << "\t" << logpCutoff << "\n";
            if (l<=logpCutoff)
            {
                p2SidedFraction += exp(l-logpCutoff);
            }
        }
    }

    return exp(log(p2SidedFraction)+logpCutoff);
};


//pearson correlation applied to a 2x2 table ()
double cor(uint32_t a, uint32_t b, uint32_t c, uint32_t d)
{

    double N = a+b+c+d;
    double sumXY = a;
    double sumX = a+b;
    double sumY = a+c;

    //std::cerr << N << "\t" << sumXY << "\t" << sumX << "\t" << sumY << "\n";

    double num = sumXY-sumX*sumY/N;
    double denum = sqrt((sumX-sumX*sumX/N)*(sumY-sumY*sumY/N));
    return denum!=0 ? num/denum : 0;
};

//pearson correlation applied to 2 vectors
double cor(std::vector<double> a, std::vector<double> b)
{
    double N = a.size();
    double sumXY = 0;
    double sumXX = 0;
    double sumX = 0;
    double sumYY = 0;
    double sumY = 0;

    for (uint32_t i=0; i<a.size(); ++i)
    {
        //std::cerr << a[i] << "\t" << b[i] << "\n";
        sumXY += a[i] * b[i];
        sumXX += a[i] * a[i];
        sumX += a[i];
        sumYY += b[i] * b[i];
        sumY += b[i];
    }

    double num = sumXY-sumX*sumY/N;
    double denum = sqrt((sumXX-sumX*sumX/N)*(sumYY-sumY*sumY/N));
    return denum!=0 ? num/denum : 0;
};


bool computeQualAndBF(std::vector<std::vector<double> >& GLs,
                 std::vector<double>& pG,
                 double& qual, double& bf, uint32_t noAlleles)
{
    if (noAlleles==2)
    {
        uint32_t noGenotypes = noAlleles*(noAlleles+1)/2;
        if (noGenotypes!=pG.size())
        {
            return false;
        }


        qual = 0;
        bf = 0;
        uint32_t N = 0;
        double paa=0, pab=0, pbb=0, total=0;
        double homRefPrior=pG[0], hetPrior=pG[1], homAltPrior=pG[2];
        for (uint32_t i=0; i<GLs.size(); ++i)
        {
            //missing data
            if (GLs[i].size()==0)
            {
                continue;
            }

            total = (paa = GLs[i][0] * homRefPrior);
            total += (pab = GLs[i][1] * hetPrior);
            total += (pbb = GLs[i][2] * homAltPrior);

            qual +=  log10((pab + pbb)/total);
            bf +=  log10(total) - log10(GLs[i][0]);
            ++N;
        }

        return (N==0? false : true);
    }
    else
    {
        qual = 0;
        bf = 0;
        uint32_t N = 0;
        uint32_t noGenotypes = (noAlleles * (noAlleles+1))/2;
        double total=0;

        std::vector<double> p(noGenotypes, 0);
        for (uint32_t i=0; i<GLs.size(); ++i)
        {
            //missing data
            if (GLs[i].size()==0)
            {
                continue;
            }

            total = 0;

            for (uint32_t j=0; j<noGenotypes; ++j)
            {
                total += (p[j] =  GLs[i][j] * pG[j]);
            }

            for (uint32_t j=0; j<noGenotypes; ++j)
            {
                p[j] /= total;
            }
             ++N;
            //evidence for alternate alleles
            qual +=  log10((total-p[0])/total);
            bf +=  log10(total) - log10(GLs[i][0]);
        }

        return (N==0? false : true);
    }
};