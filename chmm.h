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

#define S  0
#define X  1
#define Y  2
#define ML 3
#define DL 4
#define IL 5
#define M  6
#define D  7
#define I  8
#define MR 9
#define DR 10
#define IR 11
#define W  12
#define Z  13
#define E  14
#define N  15
#define NSTATES 15

/*functions for getting trace back information*/
#define track_get_u(t) (((t)&0xFF000000)>>24)
#define track_get_d(t) (((t)&0x00FF0000)>>16)
#define track_get_c(t) (((t)&0x0000FF00)>>8)
#define track_get_p(t) (((t)&0x000000FF))
#define track_set_u(t,u) (((t)&0x00FFFFFF)|((u)<<24))
#define track_set_d(t,d) (((t)&0xFF00FFFF)|((d)<<16))
#define track_set_c(t,c) (((t)&0xFFFF00FF)|((c)<<8))
#define track_set_p(t,p) (((t)&0xFFFFFF00)|(p))
#define make_track(t,u,d,c,p) {\
        track_set_u(t,u); \
        track_set_d(t,d); \
        track_set_c(t,c); \
        track_set_p(t,p); }
#define NULL_TRACK 0x0F000000

class CHMM
{
    public:
    const char* read;
    const char* qual;
    const char* lflank;
    const char* ru;
    const char* rflank;

    int32_t rlen, plen, lflen, rulen, rflen;

    std::string path;
    double maxLogOdds;

    double delta;
    double epsilon;
    double tau;
    double eta;

    double logEta;
    double logTau;
    double logOneSixteenth;

    double T[NSTATES][NSTATES];

    double *V_X;
    double *V_Y;
    double *V_ML;
    double *V_IL;
    double *V_DL;
    double *V_M;
    double *V_I;
    double *V_D;
    double *V_MR;
    double *V_IR;
    double *V_DR;
    double *V_W;
    double *V_Z;

    int32_t *U_X;
    int32_t *U_Y;
    int32_t *U_ML;
    int32_t *U_IL;
    int32_t *U_DL;
    int32_t *U_M;
    int32_t *U_I;
    int32_t *U_D;
    int32_t *U_MR;
    int32_t *U_IR;
    int32_t *U_DR;
    int32_t *U_W;
    int32_t *U_Z;

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
    double log10_emission_odds(char readBase, char probeBase, uint32_t pl);

    /**
     * Advance position in model.
     */
    inline int32_t advance_X(int32_t state, int32_t t)
    {
        if (track_get_d(t)==N && track_get_p(t)<lflen)
        {
            return track_set_d(t, track_get_p(t)+1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    /**
     * Prints an alignment.
     */
    std::string state2string(int32_t state)
    {
        if (state==S)
        {
            return "S";
        }
        else if (state==X)
        {
            return "X";
        }
        else if (state==Y)
        {
            return "Y";
        }
        else if (state==ML)
        {
            return "L";
        }
        else if (state==DL)
        {
            return "D";
        }
        else if (state==IL)
        {
            return "I";
        }
        else if (state==M)
        {
            return "M";
        }
        else if (state==D)
        {
            return "D";
        }
        else if (state==I)
        {
            return "I";
        }
        else if (state==MR)
        {
            return "R";
        }
        else if (state==DR)
        {
            return "D";
        }
        else if (state==IR)
        {
            return "I";
        }
        else if (state==W)
        {
            return "W";
        }
        else if (state==Z)
        {
            return "Z";
        }
        else if (state==E)
        {
            return "E";
        }
        else if (state==N)
        {
            return "N";
        }
        else
        {
            return "!";
        }
    }

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
    void print(double *v, size_t plen, size_t rlen);

    /**
     * Prints a int32_t matrix.
     */
    void print(int32_t *v, size_t plen, size_t rlen);
    
    /**
     * Prints U.
     */
    void print_U(int32_t *U, size_t plen, size_t rlen);
    
};

#undef MAXLEN
#undef S
#undef X
#undef Y
#undef ML
#undef IL
#undef DL
#undef M
#undef I
#undef D
#undef MR
#undef IR
#undef DR
#undef W
#undef Z
#undef E
#undef N
#undef NSTATES

#endif