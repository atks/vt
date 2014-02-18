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

#ifndef LHMM_H
#define LHMM_H

#include "utils.h"
#include "log_tool.h"
#include <regex.h>

//#define X 0
//#define Y 0
//#define M 0
//#define I 0
//#define D 0
//#define W 0
//#define Z 0
#define MAXLEN 2000

class LHMM
{
    public:
    const char* x;
    const char* y;
    const char* qual;

	

    double scoreX[MAXLEN][MAXLEN];
    double scoreY[MAXLEN][MAXLEN];
    double scoreM[MAXLEN][MAXLEN];
    double scoreI[MAXLEN][MAXLEN];  
	double scoreD[MAXLEN][MAXLEN];
	double scoreW[MAXLEN][MAXLEN];
    double scoreZ[MAXLEN][MAXLEN];

	double pathX[MAXLEN][MAXLEN];
    double pathY[MAXLEN][MAXLEN];
    double pathM[MAXLEN][MAXLEN];
    double pathI[MAXLEN][MAXLEN];  
	double pathD[MAXLEN][MAXLEN];
	double pathW[MAXLEN][MAXLEN];
    double pathZ[MAXLEN][MAXLEN];

    std::vector<double> PLs;

    uint32_t maxLength;
    uint32_t xlen;
    uint32_t ylen;
    std::string path;
    double maxLogOdds;

    double delta;
    double epsilon;
    double tau;
    double eta;
    double logOneSixteenth;

//    double a;
//    double e;
//    double d;
//    double c;
    double logEta;
    double logTau;

    double tsx;
    double txx;
    double tsy;
    double txy;
    double tyy;
    double txm;
    double tym;
    double txi;
    double txd;
    double tsm;
    double tmm;
    double tim;
    double tdm;
    double tmi;
    double tmd;
    double tii;
    double tdd;
    double tyi;
    double tyd;
    double tmw;
    double tiw;
    double tdw;
    double tww;
    double tmz;
    double tiz;
    double tdz;
    double twz;
    double tzz;

    double tMM;
    double tMI;
    double tMD;
    double tII;
    double tDD;
    double tDM;
    double tIM;

    int32_t matchStartX;
    int32_t matchEndX;
    int32_t matchStartY;
    int32_t matchEndY;
    int32_t matchedBases;
    int32_t mismatchedBases;

    std::vector<uint32_t> indelStartsInX;
    std::vector<uint32_t> indelEndsInX;
    std::vector<uint32_t> indelStartsInY;
    std::vector<uint32_t> indelEndsInY;
    std::vector<uint32_t> indelStartsInPath;
    std::vector<uint32_t> indelEndsInPath;
    std::vector<char> indelStatusInPath;

    uint32_t noBasesAligned;

    LogTool lt;

    /*Constructor*/
    LHMM();

    bool containsIndel();

    /**
     * Convert PLs to probabilities.
     */
    double pl2prob(uint32_t PL);

    /**
     * Updates matchStart, matchEnd, globalMaxPath and path
     * Updates locations of insertions and deletions
     */
    void tracePath();
    void tracePath(std::stringstream& ss, char state, uint32_t i, uint32_t j);

    /**
     * Align and compute genotype likelihood.
     */
    void align(double& llk, const char* _x, const char* _y, const char* qual, bool debug=false);

    //computes log likelihood based on path given in arguments
    void computeLogLikelihood(double& llk, double& perfectllk, std::string& _path, const char* qual);

    //computes log likelihood based on path given in arguments
    void computeLogLikelihood(double& llk, std::string& _path, const char* qual);

    //computes log likelihood based on path saved
    void computeLogLikelihood(double& llk, const char* qual);

    //get path
    std::string& getPath();

    /**
     * Left align indels in an alignment
     */
    void left_align();

    double logEmissionOdds(char readBase, char probeBase, double e);

    //compute log emission based on equal error probability distribution
    double logEmission(char readBase, char probeBase, double e);

    double emission(char readBase, char probeBase, double e);

    std::string reverse(std::string s);

    void printVector(double v[][MAXLEN], uint32_t xLen, uint32_t yLen);

    void printVector(char v[][MAXLEN], uint32_t xLen, uint32_t yLen);

    void printVector(double v[][MAXLEN]);

    void printAlignment();

    void printAlignment(std::string& pad);

    void printAlignment(std::string& pad, std::stringstream& log);

    double score(char a, char b);
    
    /**
     * Checks if deletion exists in alignment.
     */
    bool deletion_start_exists(uint32_t pos, uint32_t& rpos);

    /**
     * Checks if insertion exists in alignment.
     */
    bool insertion_start_exists(uint32_t pos, uint32_t& rpos);
        
};

#endif