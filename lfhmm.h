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
#ifndef LFHMM_H
#define LFHMM_H

#include <iomanip>
#include "htslib/kstring.h"
#include "log_tool.h"

#define MAXLEN 1024
#define MAXLEN_NBITS 10

#define S   0
#define ML  1
#define M   2
#define D   3
#define I   4
#define Z   5
#define E   6
#define N   7
#define TBD 8
#define NSTATES 8

#define LFLANK    0
#define MOTIF     1
#define UNMODELED 2
#define UNCERTAIN 3

//match type
#define PROBE 0
#define READ  1
#define MATCH 2


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

class LFHMM
{
    public:
        
    const char* read;
    const char* qual;
    
    //model variables
    //array indexed by LFLANK, MOTIF, RFLANK
    char **model;
    //length of read, probe and components in the model
    int32_t rlen, plen, lflen, mlen;

    /*result variables*/    
    int32_t lflank_start[2], lflank_end[2], motif_start[2], motif_end[2], rflank_start[2], rflank_end[2];
    int32_t motif_count, exact_motif_count, motif_m, motif_xid;
    int32_t* motif_discordance;
    float motif_concordance, maxLogOdds;

    //for track intermediate scores during Viterbi algorithm
    float max_score;
    int32_t max_track;
    
    //for optimal alignment
    bool optimal_path_traced;
    float optimal_score;
    int32_t optimal_state;
    int32_t optimal_track;
    int32_t optimal_probe_len;
    int32_t *optimal_path;     // for storage
    int32_t *optimal_path_ptr; //just a pointer
    int32_t optimal_path_len;
    
    float delta;
    float epsilon;
    float tau;
    float eta;

    float T[NSTATES][NSTATES];

    float **V;
    int32_t **U;

    LogTool *lt;

    typedef int32_t (LFHMM::*move) (int32_t t, int32_t j);
    move **moves;
           
    /**
     * Constructor.
     */
    LFHMM();

    /**
     * Constructor.
     */
    LFHMM(LogTool *lt);

    /**
     * Destructor.
     */
    ~LFHMM();

    /**
     * Initializes object, helper function for constructor.
     */
    void initialize(const char* lflank, const char* ru);

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
    float log10_emission_odds(char read_base, char probe_base, uint32_t pl);

    /**
     * Converts state to string representation.
     */
    std::string state2string(int32_t state);

    /**
     * Converts state to cigar string representation.
     */
    std::string state2cigarstring(int32_t state);

    /**
     * Converts track to cigar string representation.
     */
    std::string track2cigarstring1(int32_t t, int32_t j);
        
    /**
     * Converts state to cigar string representation.
     */
    std::string track2cigarstring2(int32_t t);
            
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
     * Prints a float matrix.
     */
    void print(float *v, size_t plen, size_t rlen);

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
     * Collect alignment summary statistics.
     */
    void collect_statistics(int32_t t1, int32_t t2, int32_t j);
   
    /**
     * Clear alignment statistics.
     */
    void clear_statistics();
    
    /**
     * Update alignment statistics after collection.
     */
    void update_statistics();
    
    /**
     * Returns true if flanks are mapped.
     */
    bool flanks_are_mapped();
        
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

    int32_t move_ML_ML(int32_t t, int32_t j)
    {
        int32_t p;
        if ((p=track_get_p(t))<lflen && j<rlen)
        {
            return make_track(ML,LFLANK,0,p+1);
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

    int32_t move_ML_D(int32_t t, int32_t j)
    {
        if (track_get_p(t)==lflen)
        {    
            return make_track(ML,MOTIF,1,1);
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

    int32_t move_ML_I(int32_t t, int32_t j)
    {
        if (j<rlen)
        {    
            return make_track(ML,MOTIF,1,0);
        }
        else
        {    
            return NULL_TRACK;
        }
    }
    
    int32_t move_M_I(int32_t t, int32_t j)
    {
        if (t!=NULL_TRACK && j<rlen)
        {
            return make_track(D,MOTIF,track_get_c(t)+1,1); track_set_u(t,M);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_I_I(int32_t t, int32_t j)
    {
        if (t!=NULL_TRACK && j<rlen)
        {
            return track_set_u(t,I);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    /////
    //Z//
    /////
    int32_t move_M_Z(int32_t t, int32_t j)
    {
        if (track_get_d(t)==MOTIF && track_get_p(t)==mlen && j<rlen)
        {
            return track_set_u(t,M);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_D_Z(int32_t t, int32_t j)
    {
        int32_t p;
        if (track_get_d(t)==MOTIF && track_get_p(t)==mlen && j<rlen)
        {
            return track_set_u(t,D);
        }
        else
        {
            return NULL_TRACK;
        }
    }

    int32_t move_I_Z(int32_t t, int32_t j)
    {
        int32_t p;
        if (track_get_d(t)==MOTIF && track_get_p(t)==mlen && j<rlen)
        {
            return track_set_u(t,D);
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
            return track_set_u(t,Z);
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
#undef M
#undef I
#undef D
#undef Z
#undef E
#undef N
#undef TBD
#undef NSTATES

#undef LFLANK
#undef MOTIF
#undef RFLANK
#undef UNMODELED
#undef UNCERTAIN

#undef READ
#undef PROBE
#undef MATCH

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