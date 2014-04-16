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
#include <iomanip>
#include "log_tool.h"

#define MAXLEN 256
#define MAXLEN_NBITS 8

#define S   0
#define X   1
#define Y   2
#define ML  3
#define DL  4
#define IL  5
#define M   6
#define D   7
#define I   8
#define MR  9
#define DR  10
#define IR  11
#define W   12
#define Z   13
#define E   14
#define N   15
#define TBD 16
#define NSTATES 16

#define LFLANK    0
#define MOTIF     1
#define RFLANK    2
#define READ      3
#define UNMODELED 4
#define UNCERTAIN 5

/*for indexing single array*/
#define index(i,j) (((i)<<MAXLEN_NBITS)+(j))

/*functions for getting trace back information*/
#define track_get_u(t)    (((t)&0xFF000000)>>24)
#define track_get_d(t)    (((t)&0x00FF0000)>>16)
#define track_get_c(t)    (((t)&0x0000FF00)>>8)
#define track_get_p(t)    (((t)&0x000000FF))
#define track_get_base(t) (model[track_get_d(t)][track_get_p(t)-1])
#define track_set_u(t,u)  (((t)&0x00FFFFFF)|((u)<<24))
#define track_set_d(t,d)  (((t)&0xFF00FFFF)|((d)<<16))
#define track_set_c(t,c)  (((t)&0xFFFF00FF)|((c)<<8))
#define track_set_p(t,p)  (((t)&0xFFFFFF00)|(p))
#define make_track(u,d,c,p) (((u)<<24)|((d)<<16)|((c)<<8)|(p))

#define NULL_TRACK  0x0F040000
#define START_TRACK 0x0F000000

class CHMM
{
    public:
    const char* read;
    const char* qual;
    //array indexed by LFLANK, MOTIF, RFLANK
    char **model;
    //length of read, probe and components in the model
    int32_t rlen, plen, lflen, mlen, rflen;

    std::string path;
    double maxLogOdds;

    //for track intermediate scores during Viterbi algorithm
    double max_score;
    double max_track;
    int32_t *max_path;

    double delta;
    double epsilon;
    double tau;
    double eta;

    double logEta;
    double logTau;
    double logOneSixteenth;

    double T[NSTATES][NSTATES];

    double **V;
    int32_t **U;

    LogTool *lt;

    typedef int32_t (CHMM::*move) (int32_t t, int32_t j);
    move **moves;

    std::stringstream ss;

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
     *
     *
     * @A - start state
     * @B - end state
     * @i - index
     */
    void proc_comp(int32_t A, int32_t B, int32_t i, int32_t j, int32_t match_type);

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
    double log10_emission_odds(char read_base, char probe_base, uint32_t pl);

    /**
     * Advance position in model.
     */
    inline int32_t track_next(int32_t start, int32_t end, int32_t t)
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
     * Converts state to string representation.
     */
    std::string state2string(int32_t state);

    /**
     * Converts model component to string representation.
     */
    std::string component2string(int32_t component);

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
    
    /**
     * Prints U and V.
     */
    void print_trace(int32_t state, size_t plen, size_t rlen);

    /**
     * Returns a string representation of track.
     */
    std::string track2string(int32_t t);

    /**
     * Prints track.
     */
    void print_track(int32_t t);

    /**
     * move_A_B
     *
     * A - source state
     * B - destination state
     * t - source track (model position)
     * j - source read position
     */
    int32_t move_S_X(int32_t t, int32_t j)
    {
        return make_track(S,LFLANK,0,1);
    }

    int32_t move_X_X(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen)
        {
            return make_track(X,LFLANK,0,p+1);
        }
    
        return NULL_TRACK;
    }

    int32_t move_S_Y(int32_t t, int32_t j)
    {
        return make_track(S,LFLANK,0,1);
    }

    int32_t move_X_Y(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return t;
        }
    
        return NULL_TRACK;
    }

    int32_t move_Y_Y(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return t;
        }
        
        return NULL_TRACK;
    }

    //////////////
    //Left flank//
    //////////////
    int32_t move_S_ML(int32_t t, int32_t j)
    {
        if (t==START_TRACK)
        {
            return make_track(S,LFLANK,0,1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_X_ML(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen && j<rlen)
        {
            return make_track(X,LFLANK,0,p+1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_Y_ML(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen && j<rlen)
        {
            return make_track(Y,LFLANK,0,p+1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_ML_ML(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen && j<rlen)
        {
            return make_track(ML,LFLANK,0,p+1);
        }
    
        return NULL_TRACK;
    }

    int32_t move_DL_ML(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen && j<rlen)
        {
            return make_track(DL,LFLANK,0,p+1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_IL_ML(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen && j<rlen)
        {
            return make_track(IL,LFLANK,0,p+1);
        }
    
        return NULL_TRACK;
    }

    int32_t move_ML_DL(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen)
        {
            return make_track(ML,LFLANK,0,p+1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_DL_DL(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen)
        {
            return make_track(DL,LFLANK,0,p+1);
        }
    
        return NULL_TRACK;
    }

    int32_t move_ML_IL(int32_t t, int32_t j)
    {
        int32_t p;
        if (j<rlen)
        {
            return track_set_u(t, IL);
        }
        
        return NULL_TRACK;
    }

    int32_t move_IL_IL(int32_t t, int32_t j)
    {
        int32_t p;
        if (j<rlen)
        {
            return track_set_u(t, IL);
        }
    
        return NULL_TRACK;
    }

    ///////////////////
    //Repeating Motif//
    ///////////////////

    int32_t move_S_M(int32_t t, int32_t j)
    {
        if (t==START_TRACK)
        {
            return make_track(S,MOTIF,1,1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_X_M(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))==lflen && j<rlen)
        {
            return make_track(X,MOTIF,1,1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_Y_M(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))==lflen && j<rlen)
        {
            return make_track(Y,MOTIF,1,1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_ML_M(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))==lflen && j<rlen)
        {
            return make_track(ML,MOTIF,1,1);
        }
        
        return NULL_TRACK;
    }

    int32_t move_M_M(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<=mlen && j<rlen)
        {    
            if (p==mlen)
            {    
                return make_track(M,MOTIF,track_get_c(t)+1,1);
            }
            else
            {
                return make_track(M,MOTIF,track_get_c(t),p+1);
            }
        }
        else
        {    
            return NULL_TRACK;
        }
    }

    int32_t move_D_M(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<=mlen && j<rlen)
        {    
            if (p==mlen)
            {    
                return make_track(D,MOTIF,track_get_c(t)+1,1);
            }
            else
            {
                return make_track(D,MOTIF,track_get_c(t),p+1);
            }
        }
        else
        {    
            return NULL_TRACK;
        }
    }

    int32_t move_I_M(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<=mlen)
        {    
            if (p==mlen)
            {    
                return make_track(I,MOTIF,track_get_c(t)+1,1);
            }
            else
            {
                return make_track(I,MOTIF,track_get_c(t),p+1);
            }
        }
        else
        {    
            return NULL_TRACK;
        }
    }

    int32_t move_M_D(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<=mlen)
        {    
            if (p==mlen)
            {    
                return make_track(M,MOTIF,track_get_c(t)+1,1);
            }
            else
            {
                return make_track(M,MOTIF,track_get_c(t),p+1);
            }
        }
        else
        {    
            return NULL_TRACK;
        }
    }

    int32_t move_D_D(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<=mlen)
        {    
            if (p==mlen)
            {    
                return make_track(D,MOTIF,track_get_c(t)+1,1);
            }
            else
            {
                return make_track(D,MOTIF,track_get_c(t),p+1);
            }
        }
        else
        {    
            return NULL_TRACK;
        }
    }

    int32_t move_M_I(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return track_set_u(t,M);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_I_I(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return track_set_u(t,I);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    ///////////////
    //Right flank//
    ///////////////

    int32_t move_S_MR(int32_t t, int32_t j)
    {
        return make_track(S,MOTIF,0,1);
    }

    int32_t move_X_MR(int32_t t, int32_t j)
    {
        if (track_get_d(t)==LFLANK && track_get_p(t)==lflen && j<rlen)
        {
            return make_track(Y,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }        
    }

    int32_t move_Y_MR(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            if ((track_get_d(t)==LFLANK && track_get_p(t)==lflen) ||
                (track_get_d(t)==MOTIF && track_get_p(t)==mlen))
            {
                return make_track(Y,RFLANK,0,1);
            }
            else
            {
                return NULL_TRACK;
            }
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_ML_MR(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))==lflen && j<rlen)
        {
            return make_track(ML,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_M_MR(int32_t t, int32_t j)
    {
        if (track_get_p(t)==mlen && j<rlen)
        {
            return make_track(M,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_MR_MR(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<rflen && j<rlen)
        {
            return make_track(MR,RFLANK,0,p+1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_DR_MR(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<rflen && j<rlen)
        {
            return make_track(IR,RFLANK,0,p+1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_IR_MR(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<rflen && j<rlen)
        {
            return make_track(IR,RFLANK,0,p+1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_MR_DR(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<rflen)
        {
            return make_track(MR,RFLANK,0,p+1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_DR_DR(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<rflen)
        {
            return make_track(DR,RFLANK,0,p+1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_MR_IR(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return t;
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_IR_IR(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return t;
        }
        else
        {
            return NULL_TRACK;
        }
    }

    ////////////////////////
    //Unmapped right flank//
    ////////////////////////

    int32_t move_X_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)==LFLANK && track_get_p(t)==lflen && j<rlen)
        {
            return make_track(X,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_Y_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)==LFLANK && track_get_p(t)==lflen)
        {
            return make_track(Y,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_ML_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)==LFLANK && track_get_p(t)==lflen)
        {
            return make_track(ML,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_M_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)==MOTIF && track_get_p(t)==mlen)
        {
            return make_track(M,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_D_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)==MOTIF && track_get_p(t)==mlen)
        {
            return make_track(D,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_I_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)==MOTIF && track_get_p(t)==mlen)
        {
            return make_track(I,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_MR_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)==RFLANK && track_get_p(t)<rflen)
        {
            return make_track(MR,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_W_W(int32_t t, int32_t j)
    {
        if (track_get_p(t)<rflen)
        {
            return make_track(W,RFLANK,0,1);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    /////
    //Z//
    /////
    int32_t move_MR_Z(int32_t t, int32_t j)
    {
        if (track_get_d(t)==RFLANK && track_get_p(t)==rflen && j<rlen)
        {
            return t;
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_W_Z(int32_t t, int32_t j)
    {
        int32_t p;
        if (track_get_d(t)==RFLANK && track_get_p(t)==rflen && j<rlen)
        {
            return t;
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_Z_Z(int32_t t, int32_t j)
    {
        int32_t p;
        if (j<rlen)
        {
            return t;
        }
        else
        {
            return NULL_TRACK;
        }
    }

};

#undef MAXLEN
#undef MAXLEN_NBITS
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

#undef LFLANK
#undef MOTIF
#undef RFLANK
#undef READ
#undef UNMODELED
#undef UNCERTAIN

#endif