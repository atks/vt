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
#ifndef AHMM_H
#define AHMM_H

#include <iomanip>
#include "htslib/kstring.h"
#include "log_tool.h"

#define MAXLEN 1024
#define MAXLEN_NBITS 10

#define S       0
#define M       1
#define D       2
#define I       3
#define E       4
#define N       5
#define TBD     6
#define NSTATES 6

#define LFLANK    0
#define MOTIF     1
#define RFLANK    2
#define UNMODELED 3
#define UNCERTAIN 4

//match type
#define MODEL 0
#define READ  1
#define MATCH 2

/*for indexing single array*/
#define index(i,j) (((i)<<MAXLEN_NBITS)+(j))

/*functions for getting trace back information*/
#define track_get_u(t)      (((t)&0xF8000000)>>27)
#define track_get_d(t)      (((t)&0x07000000)>>24)
#define track_get_c(t)      (((t)&0x00FFF000)>>12)
#define track_get_p(t)      (((t)&0x00000FFF))
#define track_get_base(t)   (model[track_get_d(t)][track_get_p(t)-1])
#define track_valid(t)      ((track_get_d(t)==RFLANK||track_get_d(t)==MOTIF||track_get_d(t)==LFLANK)&&track_get_p(t)!=0)
#define track_set_u(t,u)    (((t)&0x07FFFFFF)|((u)<<27))
#define track_set_d(t,d)    (((t)&0xF8FFFFFF)|((d)<<24))
#define track_set_c(t,c)    (((t)&0xFF000FFF)|((c)<<12))
#define track_set_p(t,p)    (((t)&0xFFFFF000)|(p))
#define make_track(u,d,c,p) (((u)<<27)|((d)<<24)|((c)<<12)|(p))

//[N|!|0|0]
#define NULL_TRACK  0x7C000000
//[N|l|0|0]
#define START_TRACK 0x78000000

class AHMMParameters
{
    public: 
        
    float delta;
    float epsilon;
    float tau;
    float eta;
    float mismatch_penalty;
        
    AHMMParameters()
    {
        delta = 0.0001;
        epsilon = 0.001;
        tau = 0.01;
        eta = 0.01;
        mismatch_penalty = 1;
    };
};

class AHMM
{
    public:

    const char* read;
    const char* qual;

    //model variables
    //array indexed by LFLANK, MOTIF
    char **model;
    //length of read, probe and components in the model
    int32_t rlen, plen, lflen, mlen;

    /*result variables*/
    int32_t lflank_start[2], lflank_end[2], motif_start[2], motif_end[2], rflank_start[2], rflank_end[2];
    int32_t motif_count, exact_motif_count, motif_m, motif_xid;
    int32_t last_motif_pos;
    float repeat_tract_len;
    int32_t* motif_discordance;
    float motif_concordance, maxLogOdds;
    
    AHMMParameters par;

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

    float T[NSTATES][NSTATES];

    float **V;
    int32_t **U;

    LogTool *lt;

    typedef int32_t (AHMM::*move) (int32_t t, int32_t j);
    move **moves;

    bool debug;

    /**
     * Constructor.
     */
    AHMM(bool debug=false);

    /**
     * Constructor.
     */
    AHMM(LogTool *lt, bool debug=false);

    /**
     * Destructor.
     */
    ~AHMM();

    /**
     * Initializes object, helper function for constructor.
     */
    void initialize();

    /**
     * Initializes objects for constructor.
     */
    void initialize_structures();
        
    /**
     * Initialize transition matrix based on parameters.
     */
    void initialize_T();
    
    /**
     * Initializes U and V.
     */
    void initialize_UV();

    /**
     * Sets a model.
     */
    void set_model(const char* motif);

    /**
     * Sets delta.
     */
    void set_delta(float delta);

    /**
     * Sets epsilon.
     */
    void set_epsilon(float epsilon);
    
    /**
     * Sets tau.
     */
    void set_tau(float tau);
    
    /**
     * Sets eta.
     */
    void set_eta(float eta);

    /**
     * Sets mismatch penalty.
     */
    void set_mismatch_penalty(float mismatch_penalty);

    /**
     * Sets debug.
     */
    void set_debug(bool debug);
    
    /**
     * Get left flank start position for model.
     */
    int32_t get_lflank_model_spos1();

    /**
     * Get left flank end position for model.
     */
    int32_t get_lflank_model_epos1();
    
    /**
     * Get motif start position for model.
     */
    int32_t get_motif_model_spos1();

    /**
     * Get motif end position for model.
     */
    int32_t get_motif_model_epos1();

    /**
     * Get right flank start position for model.
     */
    int32_t get_rflank_model_spos1();

    /**
     * Get right flank end position for model.
     */
    int32_t get_rflank_model_epos1();

    /**
     * Get left flank start position for read.
     */
    int32_t get_lflank_read_spos1();

    /**
     * Get left flank end position for read.
     */
    int32_t get_lflank_read_epos1();

    /**
     * Get motif start position for read.
     */
    int32_t get_motif_read_spos1();

    /**
     * Get motif end position for read.
     */
    int32_t get_motif_read_epos1();
    
    /**
     * Get right flank start position for read.
     */
    int32_t get_rflank_read_spos1();

    /**
     * Get right flank end position for read.
     */
    int32_t get_rflank_read_epos1();
    
    /**
     * Get motif concordance.
     */
    float get_motif_concordance();
    
    /**
     * Get exact motif count.
     */
    uint32_t get_exact_motif_count();
    
    /**
     * Get motif count.
     */
    uint32_t get_motif_count();
    
    /**
     * Computes the score associated with the move from A to B
     * Updates the max_score and associated max_track.
     *
     * @A      - start state
     * @B      - end state
     * @index1 - flattened index of the one dimensional array of start state
     * @j      - 1 based position of read of start state
     * @m      - base match required (MATCH, MODEL_ONLY, READ_ONLY)
     */
    void proc_comp(int32_t A, int32_t B, int32_t i, int32_t j, int32_t match_type);

    /**
     * Align and compute genotype likelihood.
     */
    void align(const char* y, const char* qual=NULL);

    /**
     * Trace path after alignment.
     */
    void trace_path();

    /**
     * Compute log10 emission odds based on equal error probability distribution contrasted against log10(1/16).
     */
    float log10_emission_odds(char read_base, char probe_base, uint32_t pl, float mismatch_penalty);

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
        if (t==START_TRACK)
        {
            return make_track(S,MOTIF,1,1);
        }

        return NULL_TRACK;
    }
    
    int32_t move_M_I(int32_t t, int32_t j)
    {
        if (t!=NULL_TRACK && j<rlen)
        {
            return track_set_u(t,M);
        }

        return NULL_TRACK;
    }

    int32_t move_I_I(int32_t t, int32_t j)
    {
        if (t!=NULL_TRACK && j<rlen)
        {
            return track_set_u(t,I);
        }

        return NULL_TRACK;
    }
};

#undef MAXLEN
#undef MAXLEN_NBITS
#undef S
#undef M
#undef I
#undef D
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
#undef MODEL
#undef MATCH

#undef index
#undef track_get_u
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