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

#include "lfhmm.h"

#define MAXLEN 1024
#define MAXLEN_NBITS 10

#define S       0
#define ML      1
#define M       2
#define D       3
#define I       4
#define Z       5
#define E       6
#define N       7
#define TBD     8
#define NSTATES 7

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

/**
 * Constructor.
 */
LFHMM::LFHMM(bool debug)
{
    lt = new LogTool();
    this->debug = debug;
    initialize();
};

/**
 * Constructor.
 */
LFHMM::LFHMM(LogTool *lt, bool debug)
{
    this->lt = lt;
    this->debug = debug;
    initialize();
};

/**
 * Destructor.
 */
LFHMM::~LFHMM()
{
    delete optimal_path;

    for (size_t state=S; state<=E; ++state)
    {
        delete V[state];
        delete U[state];
    }

    delete V;
    delete U;
};

/**
 * Initializes objects; helper function for constructor.
 */
void LFHMM::initialize()
{
    initialize_structures();
    initialize_T();
    initialize_UV();
};

/**
 * Initializes objects for constructor.
 */
void LFHMM::initialize_structures()
{
    model = new char*[3];
    model[LFLANK] = NULL;
    model[MOTIF] = NULL;
    model[RFLANK] = NULL;

    lflen = 0;
    mlen = 0;

    motif_discordance = new int32_t[MAXLEN];
    optimal_path = new int32_t[MAXLEN<<2];
    optimal_path_traced = false;

    typedef int32_t (LFHMM::*move) (int32_t t, int32_t j);
    V = new float*[NSTATES];
    U = new int32_t*[NSTATES];
    moves = new move*[NSTATES];
    for (size_t state=S; state<=E; ++state)
    {
        V[state] = new float[MAXLEN*MAXLEN];
        U[state] = new int32_t[MAXLEN*MAXLEN];
        moves[state] = new move[NSTATES];
    }

    for (size_t state=S; state<=E; ++state)
    {
        moves[state] = new move[NSTATES];
    }

    moves[S][ML]  = &LFHMM::move_S_ML;
    moves[ML][ML] = &LFHMM::move_ML_ML;

    moves[ML][M] = &LFHMM::move_ML_M;
    moves[M][M]  = &LFHMM::move_M_M;
    moves[D][M]  = &LFHMM::move_D_M;
    moves[I][M]  = &LFHMM::move_I_M;
    moves[ML][D] = &LFHMM::move_ML_D;
    moves[M][D]  = &LFHMM::move_M_D;
    moves[D][D]  = &LFHMM::move_D_D;
    moves[ML][I] = &LFHMM::move_ML_I;
    moves[M][I]  = &LFHMM::move_M_I;
    moves[I][I]  = &LFHMM::move_I_I;

    moves[M][Z] = &LFHMM::move_M_Z;
    moves[D][Z] = &LFHMM::move_D_Z;
    moves[I][Z] = &LFHMM::move_I_Z;
    moves[Z][Z] = &LFHMM::move_Z_Z;
};

/**
 * Initialize transition matrix based on parameters.
 */
void LFHMM::initialize_T()
{
    float delta = par.delta;
    float epsilon = par.epsilon;
    float tau = par.tau;
    float eta = par.eta;

    for (size_t i=S; i<=E; ++i)
    {
        for (size_t j=S; j<=E; ++j)
        {
            T[i][j] = -INFINITY;
        }
    }

    T[S][ML] = 0;
    T[ML][ML] = 0;

    T[ML][M] = log10(((1-2*delta-tau))/(eta*(1-eta)*(1-eta)));
    T[M][M] = log10(((1-2*delta-tau))/((1-eta)*(1-eta)));
    T[D][M] = log10(((1-epsilon-tau))/((1-eta)*(1-eta)));
    T[I][M] = T[D][M];

    T[ML][D] = log10((tau*delta)/(eta*(1-eta)));
    T[M][D] = log10(delta/(1-eta));
    T[D][D] = log10(epsilon/(1-eta));;

    T[ML][I] = T[ML][D];
    T[M][I] = T[M][D];
    T[I][I] = T[D][D];

    T[ML][Z] = log10((tau*tau*(1-eta))/(eta*eta*eta*(1-eta)));
    T[M][Z] = log10((tau*(1-eta))/(eta*(1-eta)));
    T[D][Z] = log10((tau*(1-eta))/(eta*(1-eta)));
    T[I][Z] = log10((tau*(1-eta))/(eta*(1-eta)));
    T[Z][Z] = log10((1-eta)/(1-eta));
};

/**
 * Initializes U and V.
 */
void LFHMM::initialize_UV()
{
    for (size_t i=0; i<MAXLEN; ++i)
    {
        for (size_t j=0; j<MAXLEN; ++j)
        {
            size_t c = index(i,j);

            V[S][c] = -INFINITY;
            U[S][c] = NULL_TRACK;

            V[ML][c] = -INFINITY;
            if (!i || !j)
            {
                U[ML][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[ML][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            V[M][c] = -INFINITY;
            if (!i || !j)
            {
                U[M][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[M][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            V[D][c] = -INFINITY;
            if (!i || !j)
            {
                U[D][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[D][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            V[I][c] = -INFINITY;
            if (!i || !j)
            {
                U[I][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[I][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            V[Z][c] = -INFINITY;
            if (!i || !j)
            {
                U[Z][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[Z][c] = make_track(TBD,UNCERTAIN,0,0);
            }
        }
    }

    V[S][index(0,0)] = 0;
    U[S][index(0,0)] = START_TRACK;

    V[ML][index(0,0)] = -INFINITY;
    U[ML][index(0,0)] = START_TRACK;

    V[M][index(0,0)] = -INFINITY;
    V[Z][index(0,0)] = -INFINITY;
};

/**
 * Sets a model.
 */
void LFHMM::set_model(const char* lflank, const char* motif)
{
    if (model[LFLANK]) free(model[LFLANK]);
    if (model[MOTIF]) free(model[MOTIF]);
    if (model[RFLANK]) free(model[RFLANK]);

    model[LFLANK] = strdup(lflank);
    model[MOTIF] = strdup(motif);

    lflen = strlen(model[LFLANK]);
    mlen = strlen(model[MOTIF]);
}

/**
 * Sets delta.
 */
void LFHMM::set_delta(float delta)
{
    par.delta = delta;
}

/**
 * Sets epsilon.
 */
void LFHMM::set_epsilon(float epsilon)
{
    par.epsilon = epsilon;
}

/**
 * Sets tau.
 */
void LFHMM::set_tau(float tau)
{
    par.tau = tau;
}

/**
 * Sets eta.
 */
void LFHMM::set_eta(float eta)
{
    par.eta = eta;
}

/**
 * Sets mismatch penalty.
 */
void LFHMM::set_mismatch_penalty(float mismatch_penalty)
{
    par.mismatch_penalty = mismatch_penalty;
}

/**
 * Sets debug.
 */
void LFHMM::set_debug(bool debug)
{
    this->debug = debug;
}

/**
 * Get left flank start position for model.
 */
int32_t LFHMM::get_lflank_model_spos1()
{
    return lflank_start[MODEL];
};

/**
 * Get left flank end position for model.
 */
int32_t LFHMM::get_lflank_model_epos1()
{
    return lflank_end[MODEL];
};

/**
 * Get motif start position for model.
 */
int32_t LFHMM::get_motif_model_spos1()
{
    return motif_start[MODEL];
};

/**
 * Get motif end position for model.
 */
int32_t LFHMM::get_motif_model_epos1()
{
    return motif_end[MODEL];
};

/**
 * Get right flank start position for model.
 */
int32_t LFHMM::get_rflank_model_spos1()
{
    return rflank_start[MODEL];
};

/**
 * Get right flank end position for model.
 */
int32_t LFHMM::get_rflank_model_epos1()
{
    return rflank_end[MODEL];
};

/**
 * Get left flank start position for read.
 */
int32_t LFHMM::get_lflank_read_spos1()
{
    return lflank_start[READ];
};

/**
 * Get left flank end position for read.
 */
int32_t LFHMM::get_lflank_read_epos1()
{
    return lflank_end[READ];
};

/**
 * Get motif start position for read.
 */
int32_t LFHMM::get_motif_read_spos1()
{
    return motif_start[READ];
};

/**
 * Get motif end position for read.
 */
int32_t LFHMM::get_motif_read_epos1()
{
    return motif_end[READ];
};

/**
 * Get right flank start position for read.
 */
int32_t LFHMM::get_rflank_read_spos1()
{
    return rflank_start[READ];
};

/**
 * Get right flank end position for read.
 */
int32_t LFHMM::get_rflank_read_epos1()
{
    return rflank_end[READ];
};

/**
 * Computes the score associated with the move from A to B
 * Updates the max_score and associated max_track.
 *
 * @A      - start state
 * @B      - end state
 * @index1 - flattened index of the one dimensional array of start state
 * @j      - 1 based position of read of start state
 * @m      - base match required (MATCH, MODEL, READ)
 */
void LFHMM::proc_comp(int32_t A, int32_t B, int32_t index1, int32_t j, int32_t match_type)
{
    float emission = 0, score = 0, valid = 0;

    //t is the new track
    int32_t  t =  (this->*(moves[A][B]))(U[A][index1],j);

    if (t==NULL_TRACK)
    {
        valid = -INFINITY;
    }
    else if (match_type==MATCH)
    {
        emission = log10_emission_odds(track_get_base(t), read[j], qual[j]-33);
    }

    score = V[A][index1] + T[A][B] + emission + valid;

    if (score>max_score)
    {
        max_score = score;
        max_track = t;
    }

    if (debug)
    {
        std::cerr << "\t" << state2string(A) << "=>" << state2string(B);
        std::cerr << " (" << ((index1-j)>>MAXLEN_NBITS) << "," << j << ") ";
        std::cerr << track2string(U[A][index1]) << "=>";
        std::cerr << track2string(t) << " ";
        std::cerr << emission << " (e: " << (track_get_d(t)<=MOTIF?track_get_base(t):'N') << " vs " << (j!=rlen?read[j]:'N')  << ") + ";
        std::cerr << T[A][B] << " (t) + ";
        std::cerr << V[A][index1] << " (p) + ";
        std::cerr << valid << " (v) = ";
        std::cerr << score << "\n";
    }
}

/**
 * Align read against model.
 */
void LFHMM::align(const char* read, const char* qual)
{
    clear_statistics();
    optimal_path_traced = false;
    this->read = read;
    this->qual = qual;
    rlen = strlen(read);
    plen = lflen + rlen;

    if (rlen>MAXLEN)
    {
        fprintf(stderr, "[%s:%d %s] Sequence to be aligned is greater than %d currently supported: %d\n", __FILE__, __LINE__, __FUNCTION__, MAXLEN, rlen);
        exit(1);
    }

    float max = 0;
    char maxPath = 'X';

    size_t c,d,u,l;

    //alignment
    //take into consideration
    for (size_t i=1; i<=plen; ++i)
    {
        //break;

        for (size_t j=1; j<=rlen; ++j)
        {
            c = index(i,j);
            d = index(i-1,j-1);
            u = index(i-1,j);
            l = index(i,j-1);

            //////
            //ML//
            //////
            if (debug) std::cerr << "(" << i << "," << j << ")\n";
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i<=lflen)
            {
                proc_comp(S, ML, d, j-1, MATCH);
                proc_comp(ML, ML, d, j-1, MATCH);
            }
            V[ML][c] = max_score;
            U[ML][c] = max_track;
            if (debug) std::cerr << "\tset ML " << max_score << " - " << track2string(max_track) << "\n";

            /////
            //M//
            /////
            //only need to update this i>rflen
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(ML, M, d, j-1, MATCH);
                proc_comp(M, M, d, j-1, MATCH);
                proc_comp(D, M, d, j-1, MATCH);
                proc_comp(I, M, d, j-1, MATCH);
            }
            V[M][c] = max_score;
            U[M][c] = max_track;
            if (debug) std::cerr << "\tset M " << max_score << " - " << track2string(max_track) << "\n";

            /////
            //D//
            /////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(ML, D, u, j, MODEL);
                proc_comp(M, D, u, j, MODEL);
                proc_comp(D, D, u, j, MODEL);
            }
            V[D][c] = max_score;
            U[D][c] = max_track;
            if (debug) std::cerr << "\tset D " << max_score << " - " << track2string(max_track) << "\n";

            /////
            //I//
            /////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(ML, I, l, j-1, READ);
                proc_comp(M, I, l, j-1, READ);
                proc_comp(I, I, l, j-1, READ);
            }
            V[I][c] = max_score;
            U[I][c] = max_track;
            if (debug) std::cerr << "\tset I " << max_score << " - " << track2string(max_track) << "\n";

            /////
            //Z//
            /////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(M, Z, l, j-1, READ);
                proc_comp(D, Z, l, j-1, READ);
                proc_comp(Z, Z, l, j-1, READ);
            }
            V[Z][c] = max_score;
            U[Z][c] = max_track;
            if (debug) std::cerr << "\tset Z " << max_score << " - " << track2string(max_track) << "\n";
        }
    }

    if (debug)
    {
        std::cerr << "\n   =V[S]=\n";
        print(V[S], plen+1, rlen+1);
        std::cerr << "\n   =U[S]=\n";
        print_U(U[S], plen+1, rlen+1);

        std::cerr << "\n   =V[ML]=\n";
        print(V[ML], plen+1, rlen+1);
        std::cerr << "\n   =U[ML]=\n";
        print_U(U[ML], plen+1, rlen+1);

        std::cerr << "\n   =V[M]=\n";
        print(V[M], plen+1, rlen+1);
        std::cerr << "\n   =U[M]=\n";
        print_U(U[M], plen+1, rlen+1);
        std::cerr << "\n   =V[D]=\n";
        print(V[D], plen+1, rlen+1);
        std::cerr << "\n   =U[D]=\n";
        print_U(U[D], plen+1, rlen+1);
        std::cerr << "\n   =V[I]=\n";
        print(V[I], plen+1, rlen+1);
        std::cerr << "\n   =U[I]=\n";
        print_U(U[I], plen+1, rlen+1);

        std::cerr << "\n   =V[Z]=\n";
        print(V[Z], plen+1, rlen+1);
        std::cerr << "\n   =U[Z]=\n";
        print_U(U[Z], plen+1, rlen+1);

        std::cerr << "\n";
        std::cerr << "\n";

        print_T();
    }

    trace_path();

    exact_motif_count = motif_count;
    motif_concordance = 0;
    for (int32_t k=1; k<=motif_count; ++k)
    {
        if (motif_discordance[k])
        {
            --exact_motif_count;
        }

        if (mlen>=motif_discordance[k])
        {
            motif_concordance += (float)(mlen-motif_discordance[k]) / mlen;
        }
    }
    motif_concordance *= 1.0/motif_count;
};

/**
 * Trace path after alignment.
 */
void LFHMM::trace_path()
{
    //search for a complete path in MR or W or Z
    size_t c;
    optimal_score = -INFINITY;
    optimal_track = NULL_TRACK;
    optimal_state = TBD;
    optimal_probe_len = 0;
    for (size_t i=lflen; i<=plen; ++i)
    {
        c = index(i,rlen);

        if (V[Z][c]>=optimal_score)
        {
            optimal_score = V[Z][c];
            optimal_track = U[Z][c];
            optimal_state = Z;
            optimal_probe_len = i;
        }

        if (V[M][c]>=optimal_score)
        {
            optimal_score = V[M][c];
            optimal_track = U[M][c];
            optimal_state = M;
            optimal_probe_len = i;
        }

        if (V[D][c]>=optimal_score)
        {
            optimal_score = V[D][c];
            optimal_track = U[D][c];
            optimal_state = D;
            optimal_probe_len = i;
        }

        if (V[I][c]>=optimal_score)
        {
            optimal_score = V[I][c];
            optimal_track = U[I][c];
            optimal_state = I;
            optimal_probe_len = i;
        }
    }

    //trace path
    optimal_path_ptr = optimal_path+(MAXLEN<<2)-1;
    int32_t i = optimal_probe_len, j = rlen;
    int32_t last_t = make_track(optimal_state, UNMODELED, 0, 0); //dummy end track for E
    optimal_path_len = 0;
    int32_t u;
    int32_t des_t, src_t = make_track(E, UNMODELED, 0, 0);

    do
    {
        u = track_get_u(last_t);
        last_t = U[u][index(i,j)];
        *optimal_path_ptr = track_set_u(last_t, u);

        des_t = *optimal_path_ptr;
        collect_statistics(src_t, des_t, j);
        if (debug) std::cerr << track2string(src_t) << " (" << i << "," << j << ") => " << track2string(des_t) << " :  " << track2string(last_t) << "\n";
        src_t = des_t;

        if (u==ML || u==M)
        {
            --i; --j;
        }
        else if (u==D)
        {
            --i;
        }
        else if (u==I || u==Z)
        {
            --j;
        }

        --optimal_path_ptr;
        ++optimal_path_len;

    } while (track_get_u(last_t)!=S);

    collect_statistics(src_t, last_t, j);

    ++optimal_path_ptr;
    optimal_path_traced = true;
};

/**
 * Collect alignment summary statistics.
 */
void LFHMM::collect_statistics(int32_t src_t, int32_t des_t, int32_t j)
{
    //std::cerr << "\t " << track2string(src_t) << " (" << j << ") => " << track2string(des_t) << "\n";

    int32_t src_u = track_get_u(src_t);
    int32_t des_u = track_get_u(des_t);

    if (src_u==E)
    {
        if (des_u==Z)
        {
            rflank_end[MODEL] = NAN;
            rflank_end[READ] = j+1;
        }
        else if (des_u==M || des_u==D || des_u==I)
        {
            rflank_start[MODEL] = NAN;
            rflank_start[READ] = NAN;
            rflank_end[MODEL] = INT32_MAX;
            rflank_end[READ] = INT32_MAX;

//            std::cerr << "SET TO INFINITY " << rflank_end[READ] << "\n";

//            std::cerr << std::setprecision(1) << std::fixed;
//            std::cerr << std::setw(8) << std::setprecision(2) << std::fixed  << ((float)rflank_start[MODEL]) << " " <<  ((float)rflank_end[MODEL]) << "\n";

            motif_end[MODEL] = track_get_c(des_t);
            motif_count = track_get_c(des_t);
            motif_end[READ] = j;

            //initialize array for tracking inexact repeats
            for (int32_t k=1; k<=motif_count; ++k)
            {
                motif_discordance[k] = 0;
            }

            if (des_u==D || track_get_base(des_t)!=read[j-1])
            {
                ++motif_discordance[motif_count];
            }
        }
    }
    else if (src_u==Z)
    {
        if (des_u==M || des_u==D)
        {
            rflank_start[MODEL] = track_get_p(src_t);
            rflank_start[READ] = j+1;

            motif_end[MODEL] = track_get_c(des_t);
            motif_count = track_get_c(des_t);
            motif_end[READ] = j;

            //initialize array for tracking inexact repeats
            for (int32_t k=1; k<=motif_count; ++k)
            {
                motif_discordance[k] = 0;
            }

            if (des_u==D || track_get_base(des_t)!=read[j-1])
            {
                ++motif_discordance[motif_count];
            }
        }
    }
    else if (src_u==M)
    {
        if (des_u==ML)
        {
            motif_start[MODEL] = track_get_c(src_t);
            motif_start[READ] = j+1;
            lflank_end[MODEL] = track_get_p(des_t);
            lflank_end[READ] = j;
        }
    }
    else if (src_u==ML)
    {
        if (des_u==S)
        {
            lflank_start[MODEL] = track_get_p(src_t);
            lflank_start[READ] = j+1;
        }
    }

    if (des_u==M)
    {
        if (track_get_base(des_t)!=read[j-1])
        {
            ++motif_discordance[track_get_c(des_t)];
        }
    }

    if (des_u==D || des_u==I)
    {
        ++motif_discordance[track_get_c(des_t)];
    }
};

/**
 * Clear alignment statistics.
 */
void LFHMM::clear_statistics()
{
    lflank_start[MODEL] = NAN;
    lflank_start[READ] = NAN;
    lflank_end[MODEL] = NAN;
    lflank_end[READ] = NAN;
    motif_start[MODEL] = NAN;
    motif_start[READ] = NAN;
    motif_end[MODEL] = NAN;
    motif_end[READ] = NAN;
    motif_count = NAN;
    exact_motif_count = NAN;
    motif_m = NAN;
    motif_xid = NAN;
    motif_concordance = NAN;
    rflank_start[MODEL] = NAN;
    rflank_start[READ] = NAN;
    rflank_end[MODEL] = NAN;
    rflank_end[READ] = NAN;
}

/**
 * Update alignment statistics after collection.
 */
void LFHMM::update_statistics()
{
    motif_concordance = (float)motif_m/(motif_m+motif_xid);
}

/**
 * Returns true if flanks are mapped.
 */
bool LFHMM::flanks_are_mapped()
{
    return lflank_end[MODEL]==lflen;
}

/**
 * Compute log10 emission odds based on equal error probability distribution.
 */
float LFHMM::log10_emission_odds(char probe_base, char read_base, uint32_t pl, float mismatch_penalty)
{
    if (read_base=='N' || probe_base=='N')
    {
        return -INFINITY;  //is this appropriate in this case?
    }

    if (read_base!=probe_base)
    {
        return lt->pl2log10_varp(pl);
    }
    else //match
    {
        return -(lt->pl2log10_varp(pl)-mismatch_penalty);
    }
};

/**
 * Compute log10 emission odds based on equal error probability distribution.
 */
float LFHMM::log10_emission_odds(char probe_base, char read_base, uint32_t pl)
{
//    if (read_base=='N' || probe_base=='N')
//    {
//        return -INFINITY;  //is this appropriate in this case?
//    }

    if (read_base!=probe_base)
    {
        return lt->pl2log10_varp(pl)-par.mismatch_penalty;
    }
    else //match
    {
        return -(lt->pl2log10_varp(pl));
    }
};

/**
 * Converts state to string representation.
 */
std::string LFHMM::state2string(int32_t state)
{
    if (state==S)
    {
        return "S";
    }
    else if (state==ML)
    {
        return "ML";
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
    else if (state==TBD)
    {
        return "*";
    }
    else
    {
        return "!";
    }
}

/**
 * Converts state to cigar string representation.
 */
std::string LFHMM::state2cigarstring(int32_t state)
{
    if (state==S)
    {
        return "S";
    }
    else if (state==ML)
    {
        return "L";
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
    else if (state==TBD)
    {
        return "*";
    }
    else
    {
        return "!";
    }
}

/**
 * Converts state to cigar string representation.
 */
std::string LFHMM::track2cigarstring1(int32_t t, int32_t j)
{
    int32_t state = track_get_u(t);

    if (state==S)
    {
        return "S";
    }
    else if (state==ML || state==M)
    {
        if (track_get_base(t)==read[j-1])
        {
            return "M";
        }
        else
        {
            return "*";
        }
    }
    else if (state==D)
    {
        return "D";
    }
    else if (state==I)
    {
        return "I";
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
    else if (state==TBD)
    {
        return "*";
    }
    else
    {
        return "!";
    }
}

/**
 * Converts track to cigar string representation.
 */
std::string LFHMM::track2cigarstring2(int32_t t)
{
    int32_t state = track_get_u(t);

    if (state==ML)
    {
        return "L";
    }
    else if (state==M)
    {
        return (track_get_c(t)%2==0?"+":"o");
    }
    else if (state==D)
    {
        return (track_get_c(t)%2==0?"+":"o");
    }
    else if (state==I)
    {
        return (track_get_c(t)%2==0?"+":"o");
    }
    else
    {
        return " ";
    }
}

/**
 * Converts model component to string representation.
 */
std::string LFHMM::component2string(int32_t component)
{
    if (component==LFLANK)
    {
        return "l";
    }
    else if (component==MOTIF)
    {
        return "m";
    }
    else if (component==RFLANK)
    {
        return "r";
    }
    else if (component==UNMODELED)
    {
        return "!";
    }
    else if (component==READ)
    {
        return "s";
    }
    else if (component==UNCERTAIN)
    {
        return "?";
    }
    else
    {
        return "!";
    }
}

/**
 * Prints an alignment.
 */
void LFHMM::print_alignment()
{
    std::string pad = "\t";
    print_alignment(pad);
};

/**
 * Prints an alignment with padding.
 */
void LFHMM::print_alignment(std::string& pad)
{
    std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    if (!optimal_path_traced)
    {
        std::cerr << "path not traced\n";
    }

    std::cerr << "\n";
//    print_T();
//    std::cerr << "\n";
    std::cerr << "QUAL+33\tMATCH\tMISMATCH\tPENALTY\n";
    int32_t qual = 'K' - 33;
    std::cerr << qual << "\t"<< log10_emission_odds('A', 'A', qual) << "\t"
                             << log10_emission_odds('A', 'C', qual) << "\t"
                             << log10_emission_odds('A', 'C', qual, 2) << "\n";
    std::cerr << "\n";

    std::cerr << "lflank       : " << model[LFLANK] << "\n";
    std::cerr << "repeat motif : " << model[MOTIF] << "\n";
    std::cerr << "lflen        : " << lflen << "\n";
    std::cerr << "mlen         : " << mlen << "\n";
    std::cerr << "plen         : " << plen << "\n";
    std::cerr << "\n";
    std::cerr << "read         : " << read << "\n";
    std::cerr << "rlen         : " << rlen << "\n";
    std::cerr << "\n";
    std::cerr << "optimal score: " << optimal_score << "\n";
    std::cerr << "optimal state: " << state2string(optimal_state) << "\n";
    std::cerr << "optimal track: " << track2string(optimal_track) << "\n";
    std::cerr << "optimal probe len: " << optimal_probe_len << "\n";
    std::cerr << "optimal path length : " << optimal_path_len << "\n";
    std::cerr << "max j: " << rlen << "\n";
    std::cerr << "mismatch penalty: " << par.mismatch_penalty << "\n";
    std::cerr << "\n";

    std::cerr << "model: " << "(" << lflank_start[MODEL] << "~" << lflank_end[MODEL] << ") "
                          << "[" << motif_start[MODEL] << "~" << motif_end[MODEL] << "]\n";
    std::cerr << "read : " << "(" << lflank_start[READ] << "~" << lflank_end[READ] << ") "
                          << "[" << motif_start[READ] << "~" << motif_end[READ] << "]"
                          << "[" << rflank_start[READ] << "~" << rflank_end[READ] << "]\n";
    std::cerr << "\n";
    std::cerr << "motif #           : " << motif_count << " [" << motif_start[READ] << "," << motif_end[READ] << "]\n";

    std::cerr << "motif concordance : " << motif_concordance << "% (" << exact_motif_count << "/" << motif_count << ")\n";
    std::cerr << "motif discordance : ";
    for (int32_t k=1; k<=motif_count; ++k)
    {
        std::cerr << motif_discordance[k] << (k==motif_count?"\n":"|");
    }
    std::cerr << "\n";

    //print path
    int32_t* path;
    path = optimal_path_ptr;
    std::cerr << "Model:  ";
    int32_t t = NULL_TRACK;
    int32_t j = 0;
    while (path<optimal_path+(MAXLEN<<2))
    {
        int32_t u = track_get_u(*path);
        if (u==ML || u==M || u==D)
        {
            std::cerr << track_get_base(*path);
        }
        else
        {
            std::cerr << '-';
        }
        ++path;
    }
    std::cerr << " \n";

    std::cerr << "       S";
    path = optimal_path_ptr;
    j=1;
    while (path<optimal_path+(MAXLEN<<2))
    {
        std::cerr << track2cigarstring1(*path,j);
        int32_t u = track_get_u(*path);
        if (u==ML || u==M || u==I || u==Z)
        {
            ++j;
        }
        ++path;
    }
    std::cerr << "E\n";

    path = optimal_path_ptr;
    std::cerr << "        ";
    while (path<optimal_path+(MAXLEN<<2))
    {
        std::cerr << track2cigarstring2(*path);
        ++path;
    }
    std::cerr << " \n";

    path = optimal_path_ptr;
    j=1;
    std::cerr << "Read:   ";
    while (path<optimal_path+(MAXLEN<<2))
    {
        int32_t u = track_get_u(*path);
        if (u==ML || u==M || u==I || u==Z)
        {
            std::cerr << read[j-1];
            ++j;
        }
        else
        {
            std::cerr << '-';
        }
        ++path;
    }
    std::cerr << " \n";
    std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
};

/**
 * Prints a float matrix.
 */
void LFHMM::print(float *v, size_t plen, size_t rlen)
{
    float val;
    std::cerr << std::setprecision(1) << std::fixed;
    for (size_t i=0; i<plen; ++i)
    {
        for (size_t j=0; j<rlen; ++j)
        {
            val =  v[index(i,j)];
            std::cerr << (val<0?"  ":"   ") << val;
        }

        std::cerr << "\n";
    }
};

/**
 * Prints a char matrix.
 */
void LFHMM::print(int32_t *v, size_t plen, size_t rlen)
{
    float val;
    std::cerr << std::setprecision(1) << std::fixed << std::setw(6);
    for (size_t i=0; i<plen; ++i)
    {
        for (size_t j=0; j<rlen; ++j)
        {
          val =  v[index(i,j)];
          std::cerr << (val<0?"  ":"   ") << val;
        }

        std::cerr << "\n";
    }
};

/**
 * Prints the transition matrix.
 */
void LFHMM::print_T()
{
    for (size_t j=S; j<=Z; ++j)
    {
        std::cerr << std::setw(8) << std::setprecision(2) << std::fixed << state2string(j);
    }
    std::cerr << "\n";

    for (size_t i=S; i<=Z; ++i)
    {
        for (size_t j=S; j<=Z; ++j)
        {
            if (j)
            {
                std::cerr << std::setw(8) << std::setprecision(2) << std::fixed << T[i][j];
            }
            else
            {
                std::cerr << state2string(i) << std::setw(8) << std::setprecision(2) << std::fixed << T[i][j];
            }
        }
        std::cerr << "\n";
    }
};

/**
 * Prints U.
 */
void LFHMM::print_U(int32_t *U, size_t plen, size_t rlen)
{
    std::cerr << std::setprecision(1) << std::fixed;
    std::string state;
    for (size_t i=0; i<plen; ++i)
    {
        for (size_t j=0; j<rlen; ++j)
        {
            int32_t t = U[index(i,j)];
            state = state2string(track_get_u(t));
            std::cerr << (state.size()==1 ? "   " : "  ")
                      << state << "|"
                      << component2string(track_get_d(t)) << "|"
                      << track_get_c(t) << "|"
                      << track_get_p(t) << (j==rlen-1?"\n":"   ");
        }
    }
};

/**
 * Prints U and V.
 */
void LFHMM::print_trace(int32_t state, size_t plen, size_t rlen)
{
    std::cerr << std::setprecision(1) << std::fixed;
    int32_t *u = U[state];
    float *v = V[state];
    std::string s;
    for (size_t i=0; i<plen; ++i)
    {
        for (size_t j=0; j<rlen; ++j)
        {
            int32_t t = u[index(i,j)];
            s = state2string(track_get_u(t));
            std::cerr << (s.size()==1 ? "   " : "  ")
                      << s << "|"
                      << component2string(track_get_d(t)) << "|"
                      << track_get_c(t) << "|"
                      << track_get_p(t) << "|"
                      << v[index(i,j)];
        }

        std::cerr << "\n";
    }
};

/**
 * Returns a string representation of track.
 */
std::string LFHMM::track2string(int32_t t)
{
    kstring_t str = {0,0,0};

    kputs(state2string(track_get_u(t)).c_str(), &str);
    kputc('|', &str);
    kputs(component2string(track_get_d(t)).c_str(), &str);
    kputc('|', &str);
    kputw(track_get_c(t), &str);
    kputc('|', &str);
    kputw(track_get_p(t), &str);

    std::string s(str.s);
    if (str.m) free(str.s);

    return s;
}

/**
 * Prints track.
 */
void LFHMM::print_track(int32_t t)
{
    std::cerr << track2string(t) << "\n";
}

#undef MAXLEN
#undef MAXLEN_NBITS
#undef S
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
#undef MODEL
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