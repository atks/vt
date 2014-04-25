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

#ifndef ESTIMATORS_H
#define ESTIMATORS_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <limits.h>
#include <float.h>
#include "utils.h"
#include "hts_utils.h"
#include "log_tool.h"

/**
 * Class containing estimating methods of common features.
 */
class Estimator
{
    public:
    //log tool for log space computations.
    LogTool *lt;

    /**
     * Constructor.
     */
    Estimator()
    {
        lt = new LogTool();
    };

    /**
     * Constructor.
     */
    Estimator(LogTool *lt)
    {
        this->lt = lt;
    };

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
    void compute_hwe_af(int32_t *gts, int32_t *pls, int32_t nsamples, int32_t ploidy,
                    int32_t n_allele, float *MLE_HWE_AF, float *MLE_HWE_GF, int32_t& n, double e);

    private:
};



/**
Computes HWE allele frequencies using EM algorithm
Input:
1)Genotype Likelihoods

Output:
1) estimated HWE allele frequencies
2) sample size
*/
bool estimateHWEAlleleFrequencies(std::vector<std::vector<double> >& GLs, double eps, std::vector<double>& mleHWEAlleleFreq, uint32_t& N, uint32_t noAlleles);


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
                                 uint32_t& N, uint32_t noAlleles);


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
            double& lrts, double& pValue, uint32_t& dof, uint32_t noAlleles);

/**
 Estimate Fis.
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
                 double& F, uint32_t noAlleles);

/**
Allele Balance Statistic developed by Tom Blackwell.
Works only for biallelic variants
*/
void estimateAlleleBalance(std::vector<std::vector<uint32_t> >& PLs, std::vector<std::vector<double> >& GLs, std::vector<uint32_t>& DPs, std::vector<double>& genotypeFreq, double& ab);

/**
Ratio of observe variance against expected variance - RSQ.
Works only for biallelic variants
*/
bool computeRSQ(std::vector<std::vector<double> >& GLs, std::vector<double>& mleHWEAlleleFreq, double& RSQ, uint32_t noAlleles);

/**
convert PLs to probabilities.
*/
double logFact(uint32_t n, std::vector<double>& LOGFACTS);

double logHypergeometricProb(std::vector<double>& logFacs, uint32_t a, uint32_t b, uint32_t c, uint32_t d, std::vector<double>& LOGFACTS);

double fisher1(uint32_t a, uint32_t b, uint32_t c, uint32_t d, std::vector<double>& LOGFACTS);

//pearson correlation applied to a 2x2 table ()
double cor(uint32_t a, uint32_t b, uint32_t c, uint32_t d);

//pearson correlation applied to 2 vectors
double cor(std::vector<double> a, std::vector<double> b);


bool computeQualAndBF(std::vector<std::vector<double> >& GLs,
                 std::vector<double>& pG,
                 double& qual, double& bf, uint32_t noAlleles);



#endif