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

#ifndef CHMM_H
#define CHMM_H

#include "utils.h"
#include "log_tool.h"
#include "lhmm.h"

//#define MAXLEN 250

//#define S  0
//#define X  1
//#define Y  2
//#define M  3
//#define I  4
//#define D  5
//#define W  6
//#define Z  7
//#define E  8
#define LM 9
#define LI 10
#define LD 11

#define CHMM_NSTATES 12

class CHMM
{
    public:
    const char* x;
    const char* y;
    const char* qual;

    int32_t xlen;
    int32_t ylen;
    std::string path;
    double maxLogOdds;

    double delta;
    double epsilon;
    double tau;
    double eta;
    double logOneSixteenth;

    double logEta;
    double logTau;

    double transition[CHMM_NSTATES][CHMM_NSTATES];

    //scoring matrix
    double *scoreX;
    double *scoreY;
    double *scoreM;
    double *scoreI;
	double *scoreD;
	double *scoreW;
    double *scoreZ;

    char *pathX;
    char *pathY;
    char *pathM;
    char *pathI;
    char *pathD;
    char *pathW;
    char *pathZ;

    int32_t matchStartX;
    int32_t matchEndX;
    int32_t matchStartY;
    int32_t matchEndY;
    int32_t matchedBases;
    int32_t mismatchedBases;

    uint32_t noBasesAligned;

    LogTool lt;

    CHMM();  
    
    ~CHMM()
    {   
        delete scoreX;
        delete scoreY;
        delete scoreM;
        delete scoreI;
        delete scoreD;
        delete scoreW;
        delete scoreZ;
        
        delete pathX;
        delete pathY;
        delete pathM;
        delete pathI;
        delete pathD;
        delete pathW;
        delete pathZ;
    };

    /**
     * Align and compute genotype likelihood.
     */
    void align(double& llk, const char* _x, const char* _y, const char* qual, bool debug=false);

    bool containsIndel();

    /**
     * Updates matchStart, matchEnd, globalMaxPath and path
     * Updates locations of insertions and deletions
     */
    void tracePath();
    void tracePath(std::stringstream& ss, char state, uint32_t i, uint32_t j);

    //get path
    std::string& getPath();

    double logEmissionOdds(char readBase, char probeBase, double e);

    //compute log emission based on equal error probability distribution
    double logEmission(char readBase, char probeBase, double e);

    double emission(char readBase, char probeBase, double e);

    std::string reverse(std::string s);

    void printVector(double (*v)[MAXLEN], uint32_t xLen, uint32_t yLen);

    void printVector(char v[][MAXLEN], uint32_t xLen, uint32_t yLen);

    void printVector(double v[][MAXLEN]);

    void printAlignment();

    void printAlignment(std::string& pad);

    void printAlignment(std::string& pad, std::stringstream& log);

    double score(char a, char b);


};

#endif