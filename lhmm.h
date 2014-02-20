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

#include <sstream>
#include "log_tool.h"

#define NSTATES 9

class LHMM
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

    double logEta;
    double logTau;
    double logOneSixteenth;

    double transition[NSTATES][NSTATES];

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

    //tracking of features
    int32_t matchStartX;
    int32_t matchEndX;
    int32_t matchStartY;
    int32_t matchEndY;
    int32_t matchedBases;
    int32_t mismatchedBases;

    //for left alignment
    std::vector<uint32_t> indelStartsInX;
    std::vector<uint32_t> indelEndsInX;
    std::vector<uint32_t> indelStartsInY;
    std::vector<uint32_t> indelEndsInY;
    std::vector<uint32_t> indelStartsInPath;
    std::vector<uint32_t> indelEndsInPath;
    std::vector<char> indelStatusInPath;

    std::stringstream ss;

    uint32_t noBasesAligned;

    LogTool *lt;

    /**
     * Constructor.
     */
    LHMM();

    /**
     * Constructor.
     */
    LHMM(LogTool *lt);

    /**
     * Destructor.
     */
    ~LHMM();

    /**
     * Initializes object, helper function for constructor.
     */
    void initialize();

    /**
     * Align and compute genotype likelihood.
     */
    void align(double& llk, const char* _x, const char* _y, const char* qual, bool debug=false);

    /**
     * Updates matchStart, matchEnd, globalMaxPath and path
     */
    void trace_path();

    /**
     * Recursive call for trace_path
     */
    void trace_path(char state, uint32_t i, uint32_t j);

    /**
     * Left align indels in an alignment
     */
    void left_align();

    /**
     * Compute log10 emission odds based on equal error probability distribution.
     * Substracting log10(1/16).
     */
    double log10_emission_odds(char readBase, char probeBase, double e);

    /**
     * Reverses a string.
     */
    std::string reverse(std::string s);

    /**
     * Checks if deletion exists in alignment.
     */
    bool deletion_start_exists(uint32_t pos, uint32_t& rpos);

    /**
     * Checks if insertion exists in alignment.
     */
    bool insertion_start_exists(uint32_t pos, uint32_t& rpos);

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
    void print(double *v, uint32_t xlen, uint32_t ylen);

    /**
     * Prints a char matrix.
     */
    void print(char *v, uint32_t xlen, uint32_t ylen);
};

#undef NSTATES

#endif