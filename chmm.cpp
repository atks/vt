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

#include "chmm.h"

#define MAXLEN 1024
#define MAXLEN_NBITS 10

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
#define NSTATES 15

#define LFLANK    0
#define MOTIF     1
#define RFLANK    2
#define UNMODELED 3
#define UNCERTAIN 4

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
/**
 * Constructor.
 */
CHMM::CHMM()
{
    lt = new LogTool();
};

/**
 * Constructor.
 */
CHMM::CHMM(LogTool *lt)
{
    this->lt = lt;
};

/**
 * Destructor.
 */
CHMM::~CHMM()
{
    delete optimal_path;

    //the best alignment V_ for subsequence (i,j)
    for (size_t state=X; state<=Z; ++state)
    {
        delete V[state];
        delete U[state];
    }

    delete V;
    delete U;
};

/**
 * Initializes object, helper function for constructor.
 */
void CHMM::initialize(const char* lflank, const char* motif, const char* rflank)
{
    debug = true;
    
    model = new char*[3];
    model[LFLANK] = strdup(lflank);
    model[MOTIF] = strdup(motif);
    model[RFLANK] = strdup(rflank);

    motif_discordance = new int32_t[MAXLEN];

    lflen = strlen(model[LFLANK]);
    mlen = strlen(model[MOTIF]);
    rflen = strlen(model[RFLANK]);

    optimal_path = new int32_t[MAXLEN<<2];
    optimal_path_traced = false;

    delta = 0.01;
    epsilon = 0.05;
    tau = 0.01;
    eta = 0.01;

    for (size_t i=S; i<=Z; ++i)
    {
        for (size_t j=S; j<=Z; ++j)
        {
            T[i][j] = -INFINITY;
        }
    }

    T[S][X] = 0;
    T[X][X] = 0;

    T[S][Y] = 0;
    T[X][Y] = 0;
    T[Y][Y] = 0;

    T[S][ML] = log10((1-tau)/(eta*(1-eta)*(1-eta)));
    T[X][ML] = T[S][ML];
    T[Y][ML] = T[S][ML];
    T[ML][ML] = log10(((1-2*delta-tau))/((1-eta)*(1-eta)));
    T[DL][ML] = log10(((1-epsilon))/((1-eta)*(1-eta)));
    T[IL][ML] = T[DL][ML];

    T[ML][DL] = log10((delta)/((1-eta)));
    T[DL][DL] = log10((epsilon)/((1-eta)));

    T[ML][IL] = T[ML][DL];
    T[IL][IL] = T[DL][DL];

    T[S][M] = log10((tau*(1-2*delta-tau))/(eta*eta*eta*(1-eta)*(1-eta)));
    T[X][M] = T[S][M];
    T[Y][M] = T[S][M];
    T[ML][M] = log10((tau*(1-2*delta-tau))/(eta*eta*eta*(1-eta)*(1-eta)));
    T[M][M] = log10(((1-2*delta-tau))/((1-eta)*(1-eta)));
    T[D][M] = log10(((1-epsilon-tau))/((1-eta)*(1-eta)));
    T[I][M] = T[I][M];

    T[S][D] = log10((tau*delta)/(eta*eta*(1-eta)));
    T[X][D] = T[S][D];
    T[Y][D] = T[S][D];
    T[ML][D] = log10((tau*delta)/((1-eta)));
    T[M][D] = log10(delta/(1-eta));
    T[D][D] = T[DL][DL];

    T[S][I] = log10((tau*delta)/(eta*eta*eta*(1-eta)));
    T[X][I] = T[S][I];
    T[Y][I] = T[S][I];
    T[ML][I] = log10((tau*delta)/(eta*(1-eta)));
    T[M][I] = log10(delta/(1-eta));
    T[I][I] = T[IL][IL];

    T[S][MR] = log10((tau*tau*(1-tau))/(eta*eta*eta*eta*eta*(1-eta)*(1-eta)));
    T[X][MR] = T[S][MR];
    T[Y][MR] = T[S][MR];
    T[ML][MR] = log10((tau*tau*(1-tau))/(eta*eta*eta*eta*(1-eta)*(1-eta)));
    T[M][MR] = log10((tau*(1-tau))/(eta*eta*eta*(1-eta)*(1-eta)));
    T[D][MR] = T[M][MR];
    T[I][MR] = log10((tau*(1-tau))/(eta*eta*(1-eta)*(1-eta)));
    T[MR][MR] = log10(((1-2*delta-tau))/((1-eta)*(1-eta)));
    T[DR][MR] = log10(((1-epsilon))/((1-eta)*(1-eta)));
    T[IR][MR] = T[DR][MR];

    T[MR][DR] = log10(delta/(1-eta));
    T[DR][DR] = T[DL][DL];

    T[MR][IR] = T[MR][DR];
    T[IR][IR] = T[IL][IL];

    T[S][W] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta*eta));
    T[X][W] = T[S][W];
    T[Y][W] = T[S][W];
    T[ML][W] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta));
    T[M][W] = log10((tau*tau)/(eta*eta*eta*eta));
    T[D][W] = T[M][W];
    T[I][W] = log10((tau*tau)/(eta*eta*eta));
    T[MR][W] = log10(tau/eta);
    T[W][W] = 0;

    T[MR][Z] = log10(tau/eta);
    T[W][Z] = 0;
    T[Z][Z] = 0;

    typedef int32_t (CHMM::*move) (int32_t t, int32_t j);
    V = new float*[NSTATES];
    U = new int32_t*[NSTATES];
    moves = new move*[NSTATES];
    for (size_t state=S; state<=Z; ++state)
    {
        V[state] = new float[MAXLEN*MAXLEN];
        U[state] = new int32_t[MAXLEN*MAXLEN];
        moves[state] = new move[NSTATES];
    }

    for (size_t state=S; state<E; ++state)
    {
        moves[state] = new move[NSTATES];
    }

    moves[S][X] = &CHMM::move_S_X;
    moves[X][X] = &CHMM::move_X_X;

    moves[S][Y] = &CHMM::move_S_Y;
    moves[X][Y] = &CHMM::move_X_Y;
    moves[Y][Y] = &CHMM::move_Y_Y;

    moves[S][ML] = &CHMM::move_S_ML;
    moves[X][ML] = &CHMM::move_X_ML;
    moves[Y][ML] = &CHMM::move_Y_ML;
    moves[ML][ML] = &CHMM::move_ML_ML;
    moves[DL][ML] = &CHMM::move_DL_ML;
    moves[IL][ML] = &CHMM::move_IL_ML;
    moves[ML][DL] = &CHMM::move_ML_DL;
    moves[DL][DL] = &CHMM::move_DL_DL;
    moves[ML][IL] = &CHMM::move_ML_IL;
    moves[IL][IL] = &CHMM::move_IL_IL;

    moves[X][M] = &CHMM::move_X_M;
    moves[Y][M] = &CHMM::move_Y_M;
    moves[ML][M] = &CHMM::move_ML_M;
    moves[M][M] = &CHMM::move_M_M;
    moves[D][M] = &CHMM::move_D_M;
    moves[I][M] = &CHMM::move_I_M;
    moves[ML][D] = &CHMM::move_ML_D;
    moves[M][D] = &CHMM::move_M_D;
    moves[D][D] = &CHMM::move_D_D;
    moves[ML][I] = &CHMM::move_ML_I;
    moves[M][I] = &CHMM::move_M_I;
    moves[I][I] = &CHMM::move_I_I;

    moves[X][MR] = &CHMM::move_X_MR;
    moves[Y][MR] = &CHMM::move_Y_MR;
    moves[ML][MR] = &CHMM::move_ML_MR;
    moves[M][MR] = &CHMM::move_M_MR;
    moves[D][MR] = &CHMM::move_D_MR;
    moves[I][MR] = &CHMM::move_I_MR;
    moves[MR][MR] = &CHMM::move_MR_MR;
    moves[DR][MR] = &CHMM::move_DR_MR;
    moves[IR][MR] = &CHMM::move_IR_MR;
    moves[MR][DR] = &CHMM::move_MR_DR;
    moves[DR][DR] = &CHMM::move_DR_DR;
    moves[MR][IR] = &CHMM::move_MR_IR;
    moves[IR][IR] = &CHMM::move_IR_IR;

    moves[X][W] = &CHMM::move_X_W;
    moves[Y][W] = &CHMM::move_Y_W;
    moves[ML][W] = &CHMM::move_ML_W;
    moves[M][W] = &CHMM::move_M_W;
    moves[D][W] = &CHMM::move_D_W;
    moves[I][W] = &CHMM::move_I_W;
    moves[MR][W] = &CHMM::move_MR_W;
    moves[W][W] = &CHMM::move_W_W;

    moves[MR][Z] = &CHMM::move_MR_Z;
    moves[W][Z] = &CHMM::move_W_Z;
    moves[Z][Z] = &CHMM::move_Z_Z;

    //used for back tracking, this points to the state prior to the alignment for subsequence (i,j)
    //that ends with the corresponding state

    int32_t t=0;
    for (size_t i=0; i<MAXLEN; ++i)
    {
        for (size_t j=0; j<MAXLEN; ++j)
        {
            size_t c = index(i,j);

            V[S][c] = -INFINITY;
            U[S][c] = NULL_TRACK;

            //X
            if (j) //(i,j)
            {
                V[X][c] = -INFINITY;
                U[X][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                V[X][c] = 0;
                if (i) // (i,0)
                {
                    if (i==1)
                    {
                        t = make_track(S,LFLANK,0,i);
                    }
                    else
                    {
                        t = make_track(X,LFLANK,0,i);
                    }

                    U[X][c] = t;
                }
                else // (0,0)
                {
                    U[X][c] = make_track(N,UNMODELED,0,0);
                }
            }

            //Y
            V[Y][c] = 0;
            if (i)
            {
                if (j) // (i,j)
                {
                    U[Y][c] = j==1? make_track(X,LFLANK,0,i) : make_track(Y,LFLANK,0,i);
                }
                else // (i,0)
                {
                    V[Y][c] = -INFINITY;
                    U[Y][c] = make_track(N,UNMODELED,0,0);
                }
            }
            else
            {
                if (j) // (0,j)
                {
                    U[Y][c] = j==1? make_track(S,LFLANK,0,0) : make_track(Y,LFLANK,0,0);
                }
                else // (0,0)
                {
                    U[Y][c] = make_track(N,UNMODELED,0,0);
                }
            }

            //ML
            V[ML][c] = -INFINITY;
            if (!i || !j)
            {
                U[ML][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[ML][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //DL
            V[DL][c] = -INFINITY;
            if (!i || !j)
            {
                U[DL][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[DL][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //IL
            V[IL][c] = -INFINITY;
            if (!i || !j)
            {
                U[IL][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[IL][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //M
            V[M][c] = -INFINITY;
            if (!i || !j)
            {
                U[M][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[M][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //D
            V[D][c] = -INFINITY;
            if (!i || !j)
            {
                U[D][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[D][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //I
            V[I][c] = -INFINITY;
            if (!i || !j)
            {
                U[I][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[I][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //MR
            V[MR][c] = -INFINITY;
            if (!i || !j)
            {
                U[MR][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[MR][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //DR
            V[DR][c] = -INFINITY;
            if (!i || !j)
            {
                U[DR][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[DR][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //IR
            V[IR][c] = -INFINITY;
            if (!i || !j)
            {
                U[IR][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[IR][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //W
            V[W][c] = -INFINITY;
            if (!i || !j)
            {
                U[W][c] = make_track(N,UNMODELED,0,0);
            }
            else
            {
                U[W][c] = make_track(TBD,UNCERTAIN,0,0);
            }

            //Z
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

    V[X][index(0,0)] = -INFINITY;
    U[X][index(0,0)] = START_TRACK;

    U[Y][index(1,1)] = make_track(S,LFLANK,0,1);
    V[Y][index(0,0)] = -INFINITY;
    U[Y][index(0,0)] = START_TRACK;

    V[ML][index(0,0)] = -INFINITY;
    U[ML][index(0,0)] = START_TRACK;

    V[M][index(0,0)] = -INFINITY;
    V[MR][index(0,0)] = -INFINITY;
    V[W][index(0,0)] = -INFINITY;
    V[Z][index(0,0)] = -INFINITY;
};

/**
 * Computes the score associated with the move from A to B
 * Updates the max_score and associated max_track.
 *
 * @A      - start state
 * @B      - end state
 * @index1 - flattened index of the one dimensional array of start state
 * @j      - 1 based position of read of start state
 * @m      - base match required (MATCH, PROBE, READ)
 */
void CHMM::proc_comp(int32_t A, int32_t B, int32_t index1, int32_t j, int32_t match_type)
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

    if (0)
    {
        std::cerr << "\t" << state2string(A) << "=>" << state2string(B);
        std::cerr << " (" << ((index1-j)>>MAXLEN_NBITS) << "," << j << ") ";
        std::cerr << track2string(U[A][index1]) << "=>";
        std::cerr << track2string(t) << " ";
        std::cerr << emission << " (e: " << (track_get_d(t)<=RFLANK?track_get_base(t):'N') << " vs " << (j!=rlen?read[j]:'N')  << ") + ";
        std::cerr << T[A][B] << " (t) + ";
        std::cerr << V[A][index1] << " (p) + ";
        std::cerr << valid << " (v) = ";
        std::cerr << score << "\n";
    }
}

/**
 * Align y against x.
 */
void CHMM::align(const char* read, const char* qual, bool debug)
{
    optimal_path_traced = false;
    this->read = read;
    this->qual = qual;
    rlen = strlen(read);
    plen = lflen + rlen + rflen;

    if (rlen>MAXLEN)
    {
        fprintf(stderr, "[%s:%d %s] Sequence to be aligned is greater than %d currently supported: %d\n", __FILE__, __LINE__, __FUNCTION__, MAXLEN, rlen);
        exit(1);
    }

    debug = false;

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

            /////
            //X//
            /////
            //invariant

            /////
            //Y//
            /////
            //invariant

            //////
            //ML//
            //////
            if (debug) std::cerr << "(" << i << "," << j << ")\n";
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i<=lflen)
            {
                proc_comp(S, ML, d, j-1, MATCH);
                proc_comp(X, ML, d, j-1, MATCH);
                proc_comp(Y, ML, d, j-1, MATCH);
                proc_comp(ML, ML, d, j-1, MATCH);
                proc_comp(DL, ML, d, j-1, MATCH);
                proc_comp(IL, ML, d, j-1, MATCH);
            }
            V[ML][c] = max_score;
            U[ML][c] = max_track;
            if (debug) std::cerr << "\tset ML " << max_score << " - " << track2string(max_track) << "\n";

            //////
            //DL//
            //////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i<=lflen)
            {
                proc_comp(ML, DL, u, j, PROBE);
                proc_comp(DL, DL, u, j, PROBE);
            }
            V[DL][c] = max_score;
            U[DL][c] = max_track;
            if (debug) std::cerr << "\tset DL " << max_score << " - " << track2string(max_track) << "\n";

            //////
            //IL//
            //////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>=2 && i<=lflen)
            {
                proc_comp(ML, IL, l, j-1, READ);
                proc_comp(IL, IL, l, j-1, READ);
            }
            V[IL][c] = max_score;
            U[IL][c] = max_track;
            if (debug) std::cerr << "\tset IL " << max_score << " - " << track2string(max_track) << "\n";

            /////
            //M//
            /////
            //only need to update this i>rflen
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(X, M, d, j-1, MATCH);
                proc_comp(Y, M, d, j-1, MATCH);
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
                proc_comp(ML, D, u, j, PROBE);
                proc_comp(M, D, u, j, PROBE);
                proc_comp(D, D, u, j, PROBE);
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

            //////
            //MR//
            //////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(X, MR, d, j-1, MATCH);
                proc_comp(Y, MR, d, j-1, MATCH);
                proc_comp(ML, MR, d, j-1, MATCH);
                proc_comp(M, MR, d, j-1, MATCH);
                proc_comp(D, MR, d, j-1, MATCH);
                proc_comp(I, MR, d, j-1, MATCH);
                proc_comp(MR, MR, d, j-1, MATCH);
                proc_comp(DR, MR, d, j-1, MATCH);
                proc_comp(IR, MR, d, j-1, MATCH);
            }
            V[MR][c] = max_score;
            U[MR][c] = max_track;
            if (debug) std::cerr << "\tset MR " << max_score << " - " << track2string(max_track) << "\n";

            //////
            //DR//
            //////
            if (i>lflen)
            {
                max_score = -INFINITY;
                max_track = NULL_TRACK;
            }
            proc_comp(MR, DR, u, j, PROBE);
            proc_comp(DR, DR, u, j, PROBE);
            V[DR][c] = max_score;
            U[DR][c] = max_track;
            if (debug) std::cerr << "\tset DR " << max_score << " - " << track2string(max_track) << "\n";

            //////
            //IR//
            //////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(MR, IR, l, j-1, READ);
                proc_comp(IR, IR, l, j-1, READ);
            }
            V[IR][c] = max_score;
            U[IR][c] = max_track;
            if (debug) std::cerr << "\tset IR " << max_score << " - " << track2string(max_track) << "\n";

            /////
            //W//
            /////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(X, W, u, j, PROBE);
                proc_comp(Y, W, u, j, PROBE);
                proc_comp(ML, W, u, j, PROBE);
                proc_comp(M, W, u, j, PROBE);
                proc_comp(D, W, u, j, PROBE);
                proc_comp(I, W, u, j, PROBE);
                proc_comp(MR, W, u, j, PROBE);
                proc_comp(W, W, u, j, PROBE);
            }
            V[W][c] = max_score;
            U[W][c] = max_track;
            if (debug) std::cerr << "\tset W " << max_score << " - " << track2string(max_track) << "\n";

            /////
            //Z//
            /////
            max_score = -INFINITY;
            max_track = NULL_TRACK;
            if (i>lflen)
            {
                proc_comp(MR, Z, l, j-1, READ);
                proc_comp(W, Z, l, j-1, READ);
                proc_comp(Z, Z, l, j-1, READ);
            }
            V[Z][c] = max_score;
            U[Z][c] = max_track;
            if (debug) std::cerr << "\tset Z " << max_score << " - " << track2string(max_track) << "\n";
        }
    }

    if (0)
    {
        std::cerr << "\n   =V[S]=\n";
        print(V[S], plen+1, rlen+1);
        std::cerr << "\n   =U[S]=\n";
        print_U(U[S], plen+1, rlen+1);
        std::cerr << "\n   =V[X]=\n";
        print(V[X], plen+1, rlen+1);
        std::cerr << "\n   =U[X]=\n";
        print_U(U[X], plen+1, rlen+1);
        std::cerr << "\n   =V[Y]=\n";
        print(V[Y], plen+1, rlen+1);
        std::cerr << "\n   =U[Y]=\n";
        print_U(U[Y], plen+1, rlen+1);

        std::cerr << "\n   =V[ML]=\n";
        print(V[ML], plen+1, rlen+1);
        std::cerr << "\n   =U[ML]=\n";
        print_U(U[ML], plen+1, rlen+1);
        std::cerr << "\n   =V[DL]=\n";
        print(V[DL], plen+1, rlen+1);
        std::cerr << "\n   =U[DL]=\n";
        print_U(U[DL], plen+1, rlen+1);
        std::cerr << "\n   =V[IL]=\n";
        print(V[IL], plen+1, rlen+1);
        std::cerr << "\n   =U[IL]=\n";
        print_U(U[IL], plen+1, rlen+1);

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

        std::cerr << "\n   =V[MR]=\n";
        print(V[MR], plen+1, rlen+1);
        std::cerr << "\n   =U[MR]=\n";
        print_U(U[MR], plen+1, rlen+1);
        std::cerr << "\n   =V[DR]=\n";
        print(V[DR], plen+1, rlen+1);
        std::cerr << "\n   =U[DR]=\n";
        print_U(U[DR], plen+1, rlen+1);
        std::cerr << "\n   =V[IR]=\n";
        print(V[IR], plen+1, rlen+1);
        std::cerr << "\n   =U[IR]=\n";
        print_U(U[IR], plen+1, rlen+1);

        std::cerr << "\n   =V[W]=\n";
        print(V[W], plen+1, rlen+1);
        std::cerr << "\n   =U[W]=\n";
        print_U(U[W], plen+1, rlen+1);
        std::cerr << "\n   =V[Z]=\n";
        print(V[Z], plen+1, rlen+1);
        std::cerr << "\n   =U[Z]=\n";
        print_U(U[Z], plen+1, rlen+1);

        std::cerr << "\n";

        std::cerr << "\n";
        std::cerr << "Emission log odds\n";    
        std::cerr << "K = " << 'K' << " " << log10_emission_odds('A', 'A', 'K'-33) << " " << log10_emission_odds('A', 'G', 'K'-33)  << "\n";
        std::cerr << "\n";
    }

    trace_path();

};

/**
 * Trace path after alignment.
 */
void CHMM::trace_path()
{
    //search for a complete path in MR or W or Z
    size_t c;
    optimal_score = -INFINITY;
    optimal_track = NULL_TRACK;
    optimal_state = TBD;
    optimal_probe_len = 0;
    for (size_t i=(lflen+rflen); i<=plen; ++i)
    {
        c = index(i,rlen);
        if (V[W][c]>=optimal_score)
        {
            optimal_score = V[W][c];
            optimal_track = U[W][c];
            optimal_state = W;
            optimal_probe_len = i;
        }

        if (V[Z][c]>=optimal_score)
        {
            optimal_score = V[Z][c];
            optimal_track = U[Z][c];
            optimal_state = Z;
            optimal_probe_len = i;
        }
        
        if (V[MR][c]>=optimal_score)
        {
            optimal_score = V[MR][c];
            optimal_track = U[MR][c];
            optimal_state = MR;
            optimal_probe_len = i;
        }
        
    }

    //trace path
    optimal_path_ptr = optimal_path+(MAXLEN<<2)-1;
    int32_t i = optimal_probe_len, j = rlen;
    int32_t last_t = make_track(optimal_state, RFLANK, 0, rflen+1); //dummy end track for E
    optimal_path_len = 0;
    int32_t u;
    int32_t des_t, src_t = make_track(E, RFLANK, 0, rflen+1);

    do
    {
        u = track_get_u(last_t);
        last_t = U[u][index(i,j)];
        *optimal_path_ptr = track_set_u(last_t, u);

        des_t = *optimal_path_ptr;
        collect_statistics(src_t, des_t, j);
        //std::cerr << track2string(src_t) << " (" << i << "," << j << ") => " << track2string(des_t) << " :  " << track2string(last_t) << "\n";
        src_t = des_t;

        if (u==ML || u==M || u==MR)
        {
            --i; --j;
        }
        else if (u==X || u==DL || u==D || u==DR || u==W)
        {
            --i;
        }
        else if (u==Y || u==IL || u==I || u==IR || u==Z)
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
void CHMM::collect_statistics(int32_t src_t, int32_t des_t, int32_t j)
{
    //std::cerr << "\t " << track2string(src_t) << " (" << j << ") => " << track2string(des_t) << "\n";

    int32_t src_u = track_get_u(src_t);
    int32_t des_u = track_get_u(des_t);

    if (src_u==E)
    {
        if (des_u==MR || des_u==W)
        {
            rflank_end[PROBE] = track_get_p(des_t);
            rflank_end[READ] = j;
        }
    }
    else if (src_u==Z)
    {
        if (des_u==MR || des_u==W)
        {
            rflank_end[PROBE] = track_get_p(des_t);
            rflank_end[READ] = j;
        }
    }
    else if (src_u==MR)
    {
        if (des_u==M || des_u==D)
        {
            rflank_start[PROBE] = track_get_p(src_t);
            rflank_start[READ] = j+1;

            motif_end[PROBE] = track_get_c(des_t);
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
        else if (des_u==ML || des_u==X)
        {
            rflank_start[PROBE] = track_get_p(src_t);
            rflank_start[READ] = j+1;
            motif_start[PROBE] = 0;
            motif_start[READ] = 0;
            motif_end[PROBE] = 0;
            motif_end[READ] = 0;
            lflank_end[PROBE] = track_get_p(des_t);
            lflank_end[READ] = j;
        }
    }
    else if (src_u==M)
    {
        if (des_u==ML || des_u==X)
        {
            motif_start[PROBE] = track_get_c(src_t);
            motif_start[READ] = j+1;
            lflank_end[PROBE] = track_get_p(des_t);
            lflank_end[READ] = j;
        }
    }
    else if (src_u==ML)
    {
        if (des_u==S)
        {
            lflank_start[PROBE] = track_get_p(src_t);
            lflank_start[READ] = j+1;
        }
        else if (des_u==X || des_u==Y)
        {
            lflank_start[PROBE] = track_get_p(src_t);
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
void CHMM::clear_statistics()
{
    lflank_start[PROBE] = -1;
    lflank_start[READ] = -1;
    lflank_end[PROBE] = -1;
    lflank_end[READ] = -1;
    motif_start[PROBE] = -1;
    motif_start[READ] = -1;
    motif_end[PROBE] = -1;
    motif_end[READ] = -1;
    rflank_start[PROBE] = -1;
    rflank_start[READ] = -1;
    rflank_end[PROBE] = -1;
    rflank_end[READ] = -1;
    motif_count = 0;
    exact_motif_count = 0;
    motif_m = 0;
    motif_xid = 0;
    motif_concordance = 0;
}

/**
 * Update alignment statistics after collection.
 */
void CHMM::update_statistics()
{
    motif_concordance = (float)motif_m/(motif_m+motif_xid);
}

/**
 * Returns true if flanks are mapped.
 */
bool CHMM::flanks_are_mapped()
{
    return lflank_end[PROBE]==lflen && rflank_start[PROBE]==rflen;
}

/**
 * Compute log10 emission odds based on equal error probability distribution.
 */
float CHMM::log10_emission_odds(char probe_base, char read_base, uint32_t pl)
{
    //4 encodes for N
    if (read_base=='N' || probe_base=='N')
    {
        //silent match
        return -INFINITY;
    }

    if (read_base!=probe_base)
    {
        return lt->pl2log10_varp(pl);
    }
    else
    {
        return -lt->pl2log10_varp(pl);
    }
};

/**
 * Converts state to string representation.
 */
std::string CHMM::state2string(int32_t state)
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
        return "ML";
    }
    else if (state==DL)
    {
        return "DL";
    }
    else if (state==IL)
    {
        return "IL";
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
        return "MR";
    }
    else if (state==DR)
    {
        return "DR";
    }
    else if (state==IR)
    {
        return "IR";
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
std::string CHMM::state2cigarstring(int32_t state)
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
        return "l";
    }
    else if (state==IL)
    {
        return "i";
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
        return "r";
    }
    else if (state==IR)
    {
        return "i";
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
std::string CHMM::track2cigarstring1(int32_t t, int32_t j)
{
    int32_t state = track_get_u(t);

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
    else if (state==ML || state==M || state==MR)
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
    else if (state==DL || state==D || state==DR)
    {
        return "D";
    }
    else if (state==IL || state==I || state==IR)
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
std::string CHMM::track2cigarstring2(int32_t t)
{
    int32_t state = track_get_u(t);

    if (state==ML || state==DL || state==IL)
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
    else if (state==MR || state==DR || state==IR)
    {
        return "R";
    }
    else
    {
        return " ";
    }
}
/**
 * Converts model component to string representation.
 */
std::string CHMM::component2string(int32_t component)
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
void CHMM::print_alignment()
{
    std::string pad = "\t";
    print_alignment(pad);
};

/**
 * Prints an alignment with padding.
 */
void CHMM::print_alignment(std::string& pad)
{
    if (!optimal_path_traced)
    {
        std::cerr << "path not traced\n";
    }

    std::cerr << "rflank       : " << model[LFLANK] << "\n";
    std::cerr << "repeat motif : " << model[MOTIF] << "\n";
    std::cerr << "lflank       : " << model[RFLANK] << "\n";
    std::cerr << "lflen        : " << lflen << "\n";
    std::cerr << "mlen         : " << mlen << "\n";
    std::cerr << "rflen        : " << rflen << "\n";
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
    std::cerr << "optimal path     : " << optimal_path << "\n";
    std::cerr << "optimal path ptr : " << optimal_path_ptr  << "\n";
    std::cerr << "max j: " << rlen << "\n";

    std::cerr << "probe: " << "(" << lflank_start[PROBE] << "~" << lflank_end[PROBE] << ") "
                          << "[" << motif_start[PROBE] << "~" << motif_end[PROBE] << "] "
                          << "(" << rflank_start[PROBE] << "~" << rflank_end[PROBE] << ")\n";
    std::cerr << "read : " << "(" << lflank_start[READ] << "~" << lflank_end[READ] << ") "
                          << "[" << motif_start[READ] << "~" << motif_end[READ] << "] "
                          << "(" << rflank_start[READ] << "~" << rflank_end[READ] << ")\n";
    std::cerr << "\n";
    std::cerr << "motif #           : " << motif_count << " [" << motif_start[READ] << "," << motif_end[READ] << "]\n";

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
    motif_concordance *= 100.0/motif_count;
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
        if (u==X || u==ML || u==DL || u==M || u==D || u==MR || u==DR || u==W)
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
        if (u==Y || u==ML || u==IL || u==M || u==I || u==MR || u==IR || u==Z)
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
        if (u==Y || u==ML || u==IL || u==M || u==I || u==MR || u==IR || u==Z)
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
        
        
        
};

/**
 * Prints a float matrix.
 */
void CHMM::print(float *v, size_t plen, size_t rlen)
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
void CHMM::print(int32_t *v, size_t plen, size_t rlen)
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
void CHMM::print_T()
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
void CHMM::print_U(int32_t *U, size_t plen, size_t rlen)
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
void CHMM::print_trace(int32_t state, size_t plen, size_t rlen)
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
std::string CHMM::track2string(int32_t t)
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
void CHMM::print_track(int32_t t)
{
    std::cerr << track2string(t) << "\n";
}

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