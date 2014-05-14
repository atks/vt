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
#ifndef RFHMM_H
#define RFHMM_H

#include <iomanip>
#include "htslib/kstring.h"
#include "log_tool.h"

#define MAXLEN 1024
#define MAXLEN_NBITS 10

#define S   0
#define Y   1
#define M   2
#define D   3
#define I   4
#define E   5
#define N   6
#define TBD 7
#define NSTATES 6

#define MOTIF     0
#define READ      1
#define UNMODELED 2
#define UNCERTAIN 3

//match type
#define MATCH      0
#define READ_ONLY  1
#define PROBE_ONLY 2

/*for indexing single array*/
#define index(i,j) (((i)<<MAXLEN_NBITS)+(j))

/*functions for getting trace back information*/
#define track_get_u(t)    (((t)&0xF8000000)>>27)
#define track_get_d(t)    (((t)&0x07000000)>>24)
#define track_get_c(t)    (((t)&0x00FFF000)>>12)
#define track_get_p(t)    (((t)&0x00000FFF))
#define track_get_base(t) (model[track_get_d(t)][track_get_p(t)-1])
#define track_valid(t) ((track_get_d(t)==RFLANK||track_get_d(t)==MOTIF||track_get_d(t)==LFLANK)&&track_get_p(t)!=0)
#define track_set_u(t,u)  (((t)&0x07FFFFFF)|((u)<<27))
#define track_set_d(t,d)  (((t)&0xF8FFFFFF)|((d)<<24))
#define track_set_c(t,c)  (((t)&0xFF000FFF)|((c)<<12))
#define track_set_p(t,p)  (((t)&0xFFFFF000)|(p))
#define make_track(u,d,c,p) (((u)<<27)|((d)<<24)|((c)<<12)|(p))

//[N|!|0|0]
#define NULL_TRACK  0x7C000000
//[N|l|0|0]
#define START_TRACK 0x78000000

class RFHMM
{
    public:
    const char* read;
    const char* qual;
    //array indexed by LFLANK, MOTIF, RFLANK
    char **model;
    //length of read, probe and components in the model
    int32_t rlen, plen, mlen;

    double maxLogOdds;

    //for track intermediate scores during Viterbi algorithm
    double max_score;
    int32_t max_track;
    
    //for optimal alignment
    bool optimal_path_traced;
    double optimal_score;
    int32_t optimal_state;
    int32_t optimal_track;
    int32_t optimal_probe_len;
    int32_t *optimal_path; // for storage
    int32_t *optimal_path_ptr; //just a pointer
    int32_t optimal_path_len;
        
    double delta;
    double epsilon;
    double tau;
    double eta;

    double logEta;
    double logTau;
    double log10OneSixteenth;

    double T[NSTATES][NSTATES];

    double **V;
    int32_t **U;

    LogTool *lt;

    typedef int32_t (RFHMM::*move) (int32_t t, int32_t j);
    move **moves;

    /**
     * Constructor.
     */
    RFHMM();

    /**
     * Constructor.
     */
    RFHMM(LogTool *lt);

    /**
     * Destructor.
     */
    ~RFHMM();

    /**
     * Initializes object, helper function for constructor.
     */
    void initialize(const char* ru);

    /**
     * Computes the score associated with the move from A to B
     * Updates the max_score and associated max_track.
     *
     * @A      - start state
     * @B      - end state
     * @index1 - flattened index of the one dimensional array of start state
     * @j      - 1 based position of read of start state
     * @m      - base match required (MATCH, PROBE_ONLY, READ_ONLY)
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
     * Converts state to string representation.
     */
    std::string state2string(int32_t state);

    /**
     * Converts state to cigar string representation.
     */
    std::string state2cigarstring(int32_t state);

    /**
     * Converts state to cigar string representation.
     */
    std::string track2cigarstring1(int32_t state, int32_t j);
        
    /**
     * Converts state to cigar string representation.
     */
    std::string track2cigarstring2(int32_t state);
            
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
     * Prints the transition matrix.
     */
    void print_T();
    
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
     * move_A_B(int32_t t, int32_t j)
     *
     * A - source state
     * B - destination state
     * t - source track (model position)
     * j - source read position
     */

    /////////////////
    //Unmapped Read//
    /////////////////
    int32_t move_S_Y(int32_t t, int32_t j)
    {
        if (track_get_p(t)==mlen && j<rlen)
        {
            return track_set_u(t,S);
        }
        
        return NULL_TRACK;
    }
    
    int32_t move_Y_Y(int32_t t, int32_t j)
    {
        int32_t p;
        if (j<rlen)
        {
            return track_set_u(t,Y);
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

    int32_t move_Y_M(int32_t t, int32_t j)
    {
        if (track_get_d(t)==MOTIF && track_get_p(t)==0)
        {
            return make_track(Y,MOTIF,1,1);
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
        
        return NULL_TRACK;
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
         
        return NULL_TRACK;
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

        return NULL_TRACK;
    }

    int32_t move_S_D(int32_t t, int32_t j)
    {
        if (t==START_TRACK)
        {    
            return make_track(S,MOTIF,1,1);
        }
            
        return NULL_TRACK;
    }
    
    int32_t move_Y_D(int32_t t, int32_t j)
    {
        if (t==START_TRACK)
        {    
            return make_track(S,MOTIF,1,1);
        }
            
        return NULL_TRACK;
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
   
        return NULL_TRACK;
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
  
        return NULL_TRACK;
    }

    int32_t move_S_I(int32_t t, int32_t j)
    {
        if (j<rlen)
        {    
            return make_track(S,MOTIF,1,0);
        }
    
        return NULL_TRACK;
    }
    
    int32_t move_Y_I(int32_t t, int32_t j)
    {
        if (j<rlen)
        {    
            return make_track(S,MOTIF,1,0);
        }
    
        return NULL_TRACK;
    }
    
    int32_t move_M_I(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return track_set_u(t,M);
        }

        return NULL_TRACK;
    }

    int32_t move_I_I(int32_t t, int32_t j)
    {
        if (j<rlen)
        {
            return track_set_u(t,I);
        }

        return NULL_TRACK;
    }
};

#undef MAXLEN
#undef MAXLEN_NBITS

#undef S
#undef Y
#undef M
#undef I
#undef D
#undef E
#undef N
#undef TBD
#undef NSTATES

#undef index

#undef MOTIF
#undef READ
#undef UNMODELED
#undef UNCERTAIN

#undef MATCH
#undef READ_ONLY
#undef PROBE_ONLY

#undef index
#undef track_get_u
#undef track_get_d
#undef track_get_d
#undef track_get_c
#undef track_get_p
#undef track_get_base
#undef track_valid
#undef track_set_u
#undef track_set_d
#undef track_set_c
#undef track_set_p
#undef make_track

#undef NULL_TRACK
#undef START_TRACK

#endif