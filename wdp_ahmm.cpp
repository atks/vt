/* The MIT License

   Copyright (c) 2016m Adrian Tan <atks@umich.edu>

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

#include "wdp_ahmm.h"

#define MAXLEN 1024
#define MAXLEN_NBITS 10

#define S       0
#define M       1
#define D       2
#define I       3
#define E       5
#define N       6
#define NSTATES 6

/*for indexing single array*/
#define index(i,j) (((i)<<MAXLEN_NBITS)+(j))

/**
 * Constructor.
 */
WDP_AHMM::WDP_AHMM(bool debug)
{
    this->debug = debug;
    initialize();
};

/**
 * Constructor.
 */
WDP_AHMM::WDP_AHMM(LogTool *lt, bool debug)
{
    this->debug = debug;
    initialize();
};

/**
 * Destructor.
 */
WDP_AHMM::~WDP_AHMM()
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
void WDP_AHMM::initialize()
{
    initialize_structures();
    initialize_T();
    initialize_UV();
};

/**
 * Initializes objects for constructor.
 */
void WDP_AHMM::initialize_structures()
{
    max_len = MAXLEN;
    motif = NULL;
    mlen = 0;

    motif_discordance = new int32_t[MAXLEN];
    optimal_path = new int32_t[MAXLEN<<2];
    optimal_path_traced = false;

    typedef int32_t (WDP_AHMM::*move) (int32_t t, int32_t j);
    V = new float*[NSTATES];
    U = new int32_t*[NSTATES];
//    moves = new move*[NSTATES];
    for (size_t state=S; state<=E; ++state)
    {
        V[state] = new float[MAXLEN*MAXLEN];
        U[state] = new int32_t[MAXLEN*MAXLEN];
    }
};

/**
 * Initialize transition matrix based on parameters.
 */
void WDP_AHMM::initialize_T()
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

    T[S][M] = log10(((1-delta))/((1-eta)*(1-eta)));
    T[M][M] = log10(((1-2*delta-tau))/((1-eta)*(1-eta)));
    T[D][M] = log10(((1-epsilon-tau))/((1-eta)*(1-eta)));
    T[I][M] = T[D][M];

    T[S][D] = log10(delta/(1-eta));
    T[M][D] = log10(delta/(1-eta));
    T[D][D] = log10(epsilon/(1-eta));;

    T[M][I] = T[M][D];
    T[I][I] = T[D][D];
};

/**
 * Initializes U and V.
 */
void WDP_AHMM::initialize_UV()
{
    for (size_t i=0; i<MAXLEN; ++i)
    {
        for (size_t j=0; j<MAXLEN; ++j)
        {
            size_t c = index(i,j);
        }
    }

    V[S][index(0,0)] = 0;
//    U[S][index(0,0)] = START_TRACK;

    V[M][index(0,0)] = -INFINITY;
};

/**
 * Sets a model.
 */
void WDP_AHMM::set_model(const char* motif)
{
//    if (motif) free(motif);

    this->motif = strdup(motif);
//    mlen = strlen(model[MOTIF]);
}

/**
 * Sets delta.
 */
void WDP_AHMM::set_delta(float delta)
{
    par.delta = delta;
}

/**
 * Sets epsilon.
 */
void WDP_AHMM::set_epsilon(float epsilon)
{
    par.epsilon = epsilon;
}

/**
 * Sets tau.
 */
void WDP_AHMM::set_tau(float tau)
{
    par.tau = tau;
}

/**
 * Sets eta.
 */
void WDP_AHMM::set_eta(float eta)
{
    par.eta = eta;
}

/**
 * Sets mismatch penalty.
 */
void WDP_AHMM::set_mismatch_penalty(float mismatch_penalty)
{
    par.mismatch_penalty = mismatch_penalty;
}

/**
 * Sets debug.
 */
void WDP_AHMM::set_debug(bool debug)
{
    this->debug = debug;
}

/**
 * Get motif start position for model.
 */
int32_t WDP_AHMM::get_motif_model_spos1()
{
    return 0;
};

/**
 * Get motif end position for model.
 */
int32_t WDP_AHMM::get_motif_model_epos1()
{
    return 0;
};

/**
 * Get motif start position for read.
 */
int32_t WDP_AHMM::get_motif_read_spos1()
{
    return 0;
};

/**
 * Get motif end position for read.
 */
int32_t WDP_AHMM::get_motif_read_epos1()
{
    return 0;
};

/**
 * Get motif concordance.
 */
float WDP_AHMM::get_motif_concordance()
{
    return motif_concordance;
};

/**
 * Get exact motif count.
 */
uint32_t WDP_AHMM::get_exact_motif_count()
{
    return exact_motif_count;
};

/**
 * Get motif count.
 */
uint32_t WDP_AHMM::get_motif_count()
{
    return motif_count;
};

/**
 * Align read against model.
 */
void WDP_AHMM::align(const char* read, const char* qual)
{
//    clear_statistics();
//    optimal_path_traced = false;
//    this->read = read;
//    this->qual = qual;
//    rlen = strlen(read);
//
//    if (rlen>=MAXLEN)
//    {
//        fprintf(stderr, "[%s:%d %s] Sequence to be aligned is greater than %d currently supported, subsetting string to first 1023 characters: %d\n", __FILE__, __LINE__, __FUNCTION__, MAXLEN, rlen);
//        rlen = 1023;
//    }
//    plen = rlen;
//
//    float max = 0;
//    char maxPath = 'X';
//
//    size_t c,d,u,l;
//
//    //alignment
//    //take into consideration
//    for (size_t i=1; i<=plen; ++i)
//    {
//        //break;
//
//        for (size_t j=1; j<=rlen; ++j)
//        {
//            c = index(i,j);
//            d = index(i-1,j-1);
//            u = index(i-1,j);
//            l = index(i,j-1);
//
//            /////
//            //M//
//            /////
//            //only need to update this i>rflen
//            max_score = -INFINITY;
//            max_track = NULL_TRACK;
//            proc_comp(S, M, d, j-1, MATCH);
//            proc_comp(M, M, d, j-1, MATCH);
//            proc_comp(D, M, d, j-1, MATCH);
//            proc_comp(I, M, d, j-1, MATCH);
//            V[M][c] = max_score;
//            U[M][c] = max_track;
//            if (debug) std::cerr << "\tset M " << max_score << " - " << track2string(max_track) << "\n";
//
//            /////
//            //D//
//            /////
//            max_score = -INFINITY;
//            max_track = NULL_TRACK;
//            proc_comp(S, D, u, j, MODEL);
//            proc_comp(M, D, u, j, MODEL);
//            proc_comp(D, D, u, j, MODEL);
//            V[D][c] = max_score;
//            U[D][c] = max_track;
//            if (debug) std::cerr << "\tset D " << max_score << " - " << track2string(max_track) << "\n";
//
//            /////
//            //I//
//            /////
//            max_score = -INFINITY;
//            max_track = NULL_TRACK;
//            proc_comp(S, I, l, j-1, READ);
//            proc_comp(M, I, l, j-1, READ);
//            proc_comp(I, I, l, j-1, READ);
//            V[I][c] = max_score;
//            U[I][c] = max_track;
//            if (debug) std::cerr << "\tset I " << max_score << " - " << track2string(max_track) << "\n";
//        }
//    }
//
//    if (debug)
//    {
//        std::cerr << "\n   =V[S]=\n";
//        print(V[S], plen+1, rlen+1);
//        std::cerr << "\n   =U[S]=\n";
//        print_U(U[S], plen+1, rlen+1);
//
//        std::cerr << "\n   =V[M]=\n";
//        print(V[M], plen+1, rlen+1);
//        std::cerr << "\n   =U[M]=\n";
//        print_U(U[M], plen+1, rlen+1);
//        std::cerr << "\n   =V[D]=\n";
//        print(V[D], plen+1, rlen+1);
//        std::cerr << "\n   =U[D]=\n";
//        print_U(U[D], plen+1, rlen+1);
//        std::cerr << "\n   =V[I]=\n";
//        print(V[I], plen+1, rlen+1);
//        std::cerr << "\n   =U[I]=\n";
//        print_U(U[I], plen+1, rlen+1);
//
//        std::cerr << "\n";
//        std::cerr << "\n";
//
//        print_T();
//    }
//
//    trace_path();
//
//    exact_motif_count = motif_count;
//    motif_concordance = 0;
//    for (int32_t k=1; k<=motif_count; ++k)
//    {
//        if (motif_discordance[k])
//        {
//            --exact_motif_count;
//        }
//
//        if (mlen>=motif_discordance[k])
//        {
//            motif_concordance += (float)(mlen-motif_discordance[k]) / mlen;
//        }
//    }
//    motif_concordance /= motif_count;
};

/**
 * Trace path after alignment.
 */
void WDP_AHMM::trace_path()
{
//    //search for a complete path in M
//    size_t c;
//    optimal_score = -INFINITY;
//    optimal_track = NULL_TRACK;
//    optimal_state = TBD;
//    optimal_probe_len = 0;
//    trf_score = 0;
//    for (size_t i=lflen; i<=plen; ++i)
//    {
//        c = index(i,rlen);
//
//        if (V[M][c]>=optimal_score)
//        {
//            optimal_score = V[M][c];
//            optimal_track = U[M][c];
//            optimal_state = M;
//            optimal_probe_len = i;
//        }
//
//        if (V[D][c]>=optimal_score)
//        {
//            optimal_score = V[D][c];
//            optimal_track = U[D][c];
//            optimal_state = D;
//            optimal_probe_len = i;
//        }
//    }
//
//    //trace path
//    optimal_path_ptr = optimal_path+(MAXLEN<<2)-1;
//    int32_t i = optimal_probe_len, j = rlen;
//    int32_t last_t = make_track(optimal_state, UNMODELED, 0, 0); //dummy end track for E
//    optimal_path_len = 0;
//    int32_t u;
//    int32_t des_t, src_t = make_track(E, UNMODELED, 0, 0);
//
//    do
//    {
//        u = track_get_u(last_t);
//        last_t = U[u][index(i,j)];
//        *optimal_path_ptr = track_set_u(last_t, u);
//
//        des_t = *optimal_path_ptr;
//        collect_statistics(src_t, des_t, j);
//        if (debug) std::cerr << track2string(src_t) << " (" << i << "," << j << ") => " << track2string(des_t) << " :  " << track2string(last_t) << "\n";
//        src_t = des_t;
//
//        if (u==M)
//        {
//            --i; --j;
//        }
//        else if (u==D)
//        {
//            --i;
//        }
//        else if (u==I || u==Z)
//        {
//            --j;
//        }
//
//        --optimal_path_ptr;
//        ++optimal_path_len;
//
//    } while (track_get_u(last_t)!=S);
//
//    collect_statistics(src_t, last_t, j);
//
//    ++optimal_path_ptr;
//    optimal_path_traced = true;
};

/**
 * Collect alignment summary statistics.
 */
void WDP_AHMM::collect_statistics(int32_t src_t, int32_t des_t, int32_t j)
{
//    //std::cerr << "\t " << track2string(src_t) << " (" << j << ") => " << track2string(des_t) << "\n";
//
//    int32_t src_u = track_get_u(src_t);
//    int32_t des_u = track_get_u(des_t);
//
//    if (src_u==E)
//    {
//        if (des_u==M || des_u==D || des_u==I)
//        {
//            rflank_start[MODEL] = NAN;
//            rflank_start[READ] = NAN;
//            rflank_end[MODEL] = INFINITY;
//            rflank_end[READ] = INFINITY;
//
//            motif_end[MODEL] = track_get_c(des_t);
//            motif_count = track_get_c(des_t);
//            last_motif_pos = track_get_p(des_t);
//            motif_end[READ] = j;
//
//            //initialize array for tracking inexact repeats
//            for (int32_t k=1; k<=motif_count; ++k)
//            {
//                motif_discordance[k] = 0;
//            }
//
//            if (des_u==D || track_get_base(des_t)!=read[j-1])
//            {
//                ++motif_discordance[motif_count];
//            }
//        }
//    }
//
//    if (des_u==M)
//    {
//        if (track_get_base(des_t)!=read[j-1])
//        {
//            trf_score -= 7;
//            ++motif_discordance[track_get_c(des_t)];
//        }
//        else
//        {
//            trf_score += 2;
//        }
//    }
//
//    if (des_u==D || des_u==I)
//    {
//        trf_score -= 7;
//        ++motif_discordance[track_get_c(des_t)];
//    }
//
//    frac_no_repeats = motif_count - (mlen-last_motif_pos)/((float)mlen);
};

/**
 * Clear alignment statistics.
 */
void WDP_AHMM::clear_statistics()
{
//    motif_start[MODEL] = NAN;
//    motif_start[READ] = NAN;
//    motif_end[MODEL] = NAN;
//    motif_end[READ] = NAN;
//    motif_count = NAN;
//    exact_motif_count = NAN;
//    motif_m = NAN;
//    motif_xid = NAN;
//    motif_concordance = NAN;
}

/**
 * Update alignment statistics after collection.
 */
void WDP_AHMM::update_statistics()
{
    motif_concordance = (float)motif_m/(motif_m+motif_xid);
}

/**
 * Compute log10 emission odds based on equal error probability distribution.
 */
float WDP_AHMM::log10_emission_odds(char probe_base, char read_base, uint32_t pl, float mismatch_penalty)
{
//    if (read_base=='N' || probe_base=='N')
//    {
//        return -INFINITY;  //is this appropriate in this case?
//    }

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
float WDP_AHMM::log10_emission_odds(char probe_base, char read_base, uint32_t pl)
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
std::string WDP_AHMM::state2string(int32_t state)
{
//    if (state==S)
//    {
//        return "S";
//    }
//    else if (state==M)
//    {
//        return "M";
//    }
//    else if (state==D)
//    {
//        return "D";
//    }
//    else if (state==I)
//    {
//        return "I";
//    }
//    else if (state==E)
//    {
//        return "E";
//    }
//    else if (state==N)
//    {
//        return "N";
//    }
//    else if (state==TBD)
//    {
//        return "*";
//    }
//    else
//    {
//        return "!";
//    }
    return "";
}

/**
 * Converts state to cigar string representation.
 */
std::string WDP_AHMM::state2cigarstring(int32_t state)
{
    if (state==S)
    {
        return "S";
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
 * Converts state to cigar string representation.
 */
std::string WDP_AHMM::track2cigarstring1(int32_t t, int32_t j)
{
//    int32_t state = track_get_u(t);
//
//    if (state==S)
//    {
//        return "S";
//    }
//    else if (state==M)
//    {
//        if (track_get_base(t)==read[j-1])
//        {
//            return "M";
//        }
//        else
//        {
//            return "*";
//        }
//    }
//    else if (state==D)
//    {
//        return "D";
//    }
//    else if (state==I)
//    {
//        return "I";
//    }
//    else if (state==Z)
//    {
//        return "Z";
//    }
//    else if (state==E)
//    {
//        return "E";
//    }
//    else if (state==N)
//    {
//        return "N";
//    }
//    else if (state==TBD)
//    {
//        return "*";
//    }
//    else
//    {
//        return "!";
//    }
    return "";
}

/**
 * Converts track to cigar string representation.
 */
std::string WDP_AHMM::track2cigarstring2(int32_t t)
{
//    int32_t state = track_get_u(t);
//
//    if (state==M)
//    {
//        return (track_get_c(t)%2==0?"+":"o");
//    }
//    else if (state==D)
//    {
//        return (track_get_c(t)%2==0?"+":"o");
//    }
//    else if (state==I)
//    {
//        return (track_get_c(t)%2==0?"+":"o");
//    }
//    else
//    {
//        return " ";
//    }
    return "";
}

/**
 * Converts model component to string representation.
 */
std::string WDP_AHMM::component2string(int32_t component)
{
//    if (component==MOTIF)
//    {
//        return "m";
//    }
//    else if (component==UNMODELED)
//    {
//        return "!";
//    }
//    else if (component==READ)
//    {
//        return "s";
//    }
//    else if (component==UNCERTAIN)
//    {
//        return "?";
//    }
//    else
//    {
//        return "!";
//    }
    return "";
}

/**
 * Prints an alignment.
 */
void WDP_AHMM::print_alignment()
{
    std::string pad = "\t";
    print_alignment(pad);
};

/**
 * Prints an alignment with padding.
 */
void WDP_AHMM::print_alignment(std::string& pad)
{
    if (!optimal_path_traced)
    {
        std::cerr << "path not traced\n";
    }

    std::cerr << "=================================\n";
    std::cerr << "WDP_AHMM\n";
    std::cerr << "*********************************\n";
    std::cerr << "repeat motif : " << motif << "\n";
    std::cerr << "mlen         : " << mlen << "\n";
//    std::cerr << "plen         : " << plen << "\n";
//    std::cerr << "\n";
    std::cerr << "read         : " << read << "\n";
    std::cerr << "rlen         : " << rlen << "\n";
//    std::cerr << "\n";
//    std::cerr << "optimal score: " << optimal_score << "\n";
//    std::cerr << "optimal state: " << state2string(optimal_state) << "\n";
//    std::cerr << "optimal track: " << track2string(optimal_track) << "\n";
//    std::cerr << "optimal probe len: " << optimal_probe_len << "\n";
//    std::cerr << "optimal path length : " << optimal_path_len << "\n";
//    std::cerr << "max j: " << rlen << "\n";
//    std::cerr << "mismatch penalty: " << par.mismatch_penalty << "\n";
//    std::cerr << "\n";
//
//    std::cerr << "model: " << "(" << lflank_start[MODEL] << "~" << lflank_end[MODEL] << ") "
//                          << "[" << motif_start[MODEL] << "~" << motif_end[MODEL] << "]\n";
//    std::cerr << "read : " << "(" << lflank_start[READ] << "~" << lflank_end[READ] << ") "
//                          << "[" << motif_start[READ] << "~" << motif_end[READ] << "]"
//                          << "[" << rflank_start[READ] << "~" << rflank_end[READ] << "]\n";
//    std::cerr << "\n";
//    std::cerr << "motif #                     : " << motif_count << " [" << motif_start[READ] << "," << motif_end[READ] << "]\n";
//    std::cerr << "motif concordance           : " << motif_concordance << "% (" << exact_motif_count << "/" << motif_count << ")\n";
//    std::cerr << "last motif position         : " << last_motif_pos << "\n";
//    std::cerr << "motif discordance           : ";
//    for (int32_t k=1; k<=motif_count; ++k)
//    {
//        std::cerr << motif_discordance[k] << (k==motif_count?"\n":"|");
//    }
//    std::cerr << "fractional no. repeat units : " << frac_no_repeats << "\n";
//    std::cerr << "repeat tract length         : " << rlen << "\n";
//    std::cerr << "TRF Score                   : " << trf_score << "\n";
//
//    std::cerr << "\n";
//
//    //print path
//    int32_t* path;
//    path = optimal_path_ptr;
//    std::cerr << "Model:  ";
//    int32_t t = NULL_TRACK;
//    int32_t j = 0;
//    while (path<optimal_path+(MAXLEN<<2))
//    {
//        int32_t u = track_get_u(*path);
//        if (u==M || u==D)
//        {
//            std::cerr << track_get_base(*path);
//        }
//        else
//        {
//            std::cerr << '-';
//        }
//        ++path;
//    }
//    std::cerr << " \n";
//
//    std::cerr << "       S";
//    path = optimal_path_ptr;
//    j=1;
//    while (path<optimal_path+(MAXLEN<<2))
//    {
//        std::cerr << track2cigarstring1(*path,j);
//        int32_t u = track_get_u(*path);
//        if (u==M || u==I || u==Z)
//        {
//            ++j;
//        }
//        ++path;
//    }
//    std::cerr << "E\n";
//
//    path = optimal_path_ptr;
//    std::cerr << "        ";
//    while (path<optimal_path+(MAXLEN<<2))
//    {
//        std::cerr << track2cigarstring2(*path);
//        ++path;
//    }
//    std::cerr << " \n";
//
//    path = optimal_path_ptr;
//    j=1;
//    std::cerr << "Read:   ";
//    while (path<optimal_path+(MAXLEN<<2))
//    {
//        int32_t u = track_get_u(*path);
//        if (u==M || u==I || u==Z)
//        {
//            std::cerr << read[j-1];
//            ++j;
//        }
//        else
//        {
//            std::cerr << '-';
//        }
//        ++path;
//    }
//    std::cerr << " \n";
//    std::cerr << "=================================\n";
};

/**
 * Prints a float matrix.
 */
void WDP_AHMM::print(float *v, size_t plen, size_t rlen)
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
void WDP_AHMM::print(int32_t *v, size_t plen, size_t rlen)
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
void WDP_AHMM::print_T()
{
//    std::cerr << "\t";
//    for (size_t j=S; j<=Z; ++j)
//    {
//        std::cerr << std::setw(8) << std::setprecision(2) << std::fixed << state2string(j);
//    }
//    std::cerr << "\n";
//
//    for (size_t i=S; i<=Z; ++i)
//    {
//        std::cerr << "\t";
//        for (size_t j=S; j<=Z; ++j)
//        {
//            if (j)
//            {
//                std::cerr << std::setw(8) << std::setprecision(2) << std::fixed << T[i][j];
//            }
//            else
//            {
//                std::cerr << state2string(i) << std::setw(8) << std::setprecision(2) << std::fixed << T[i][j];
//            }
//        }
//        std::cerr << "\n";
//    }
//    std::cerr << "\n";
};

/**
 * Prints U.
 */
void WDP_AHMM::print_U(int32_t *U, size_t plen, size_t rlen)
{
//    std::cerr << std::setprecision(1) << std::fixed;
//    std::string state;
//    for (size_t i=0; i<plen; ++i)
//    {
//        for (size_t j=0; j<rlen; ++j)
//        {
//            int32_t t = U[index(i,j)];
//            state = state2string(track_get_u(t));
//            std::cerr << (state.size()==1 ? "   " : "  ")
//                      << state << "|"
//                      << component2string(track_get_d(t)) << "|"
//                      << track_get_c(t) << "|"
//                      << track_get_p(t) << (j==rlen-1?"\n":"   ");
//        }
//    }
};

/**
 * Prints U and V.
 */
void WDP_AHMM::print_trace(int32_t state, size_t plen, size_t rlen)
{
//    std::cerr << std::setprecision(1) << std::fixed;
//    int32_t *u = U[state];
//    float *v = V[state];
//    std::string s;
//    for (size_t i=0; i<plen; ++i)
//    {
//        for (size_t j=0; j<rlen; ++j)
//        {
//            int32_t t = u[index(i,j)];
//            s = state2string(track_get_u(t));
//            std::cerr << (s.size()==1 ? "   " : "  ")
//                      << s << "|"
//                      << component2string(track_get_d(t)) << "|"
//                      << track_get_c(t) << "|"
//                      << track_get_p(t) << "|"
//                      << v[index(i,j)];
//        }
//
//        std::cerr << "\n";
//    }
};

/**
 * Prints track.
 */
void WDP_AHMM::print_track(int32_t t)
{
//    std::cerr << track2string(t) << "\n";
}

#undef MAXLEN
#undef MAXLEN_NBITS
#undef S
#undef M
#undef I
#undef D
#undef E
#undef N
#undef NSTATES

#undef index
