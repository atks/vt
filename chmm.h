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

#include <sstream>
#include "log_tool.h"

#define NSTATES 15

class CHMM
{
    public:
    const char* read;
    const char* qual;
    const char* lflank;
    const char* ru;
    const char* rflank;
    
    int32_t rlen, lflen, rulen, rflen;

    std::string path;
    double maxLogOdds;

    double delta;
    double epsilon;
    double tau;
    double eta;

    double logEta;
    double logTau;
    double logOneSixteenth;

    double transition[NSTATES][NSTATES];

    double *scoreX;
    double *scoreY;
    double *scoreLM;
    double *scoreLI;
    double *scoreLD;
    double *scoreM;
    double *scoreI;
    double *scoreD;
    double *scoreRM;
    double *scoreRI;
    double *scoreRD;
    double *scoreW;
    double *scoreZ;

    int32_t *positionLM;
    int32_t *positionLI;
    int32_t *positionLD;
    int32_t *positionM;
    int32_t *positionI;
    int32_t *positionD;
    int32_t *positionRM;
    int32_t *positionRI;
    int32_t *positionRD;
    
    char *pathX;
    char *pathY;
    char *pathLM;
    char *pathLI;
    char *pathLD;
    char *pathM;
    char *pathI;
    char *pathD;
    char *pathRM;
    char *pathRI;
    char *pathRD;
    char *pathW;
    char *pathZ;

    LogTool *lt;

    /**
     * Constructor.
     */
    CHMM();

    /**
     * Constructor.
     */
    CHMM(LogTool *lt);

    /**
     * Destructor.
     */
    ~CHMM();

    /**
     * Initializes object, helper function for constructor.
     */
    void initialize(const char* lflank, const char* ru, const char* rflank);

    /**
     * Align and compute genotype likelihood.
     */
    void align(const char* y, const char* qual=NULL, bool debug=false);

    /**
     * Trace path after alignment.
     */
    void trace_path();

    /**
     * Compute log10 emission odds based on equal error probability distribution contrasted against log10(1/16).
     */
    double log10_emission_odds(char readBase, char probeBase, double e);

    /**
     * Prints an alignment.
     */
    void print_alignment();

    /**
     * Prints an alignment with padding.
     */
    void print_alignment(std::string& pad);

    /**
     * Prints a double matrix.
     */
    void print(double *v, uint32_t rlen);

    /**
     * Prints a char matrix.
     */
    void print(char *v, uint32_t rlen);
};

#undef NSTATES

#endif