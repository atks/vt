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
#define MAXLEN 250
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

#define LFLANK 0
#define MOTIF  1
#define RFLANK 2
#define UNMODELED 3
#define UNCERTAIN 4

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
void CHMM::initialize(const char* lflank, const char* ru, const char* rflank)
{
    this->lflank = lflank;
    this->ru = ru;
    this->rflank = rflank;

    delta = 0.001;
    epsilon = 0.5;
    tau = 0.1;
    eta = 0.001;

    logOneSixteenth = log10(1.0/16.0);

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
    T[M][M] = log10((tau*(1-2*delta-tau))/(eta*eta*(1-eta)*(1-eta)));
    T[D][M] = log10(((1-epsilon-tau))/((1-eta)*(1-eta)));
    T[I][M] = T[I][M];

    T[S][D] = log10((tau*delta)/(eta*eta*(1-eta)));
    T[X][D] = T[S][D];
    T[Y][D] = T[S][D];
    T[ML][D] = log10((tau*delta)/((1-eta)));
    T[M][D] = log10(delta/(1-eta));

    T[S][I] = log10((tau*delta)/(eta*eta*eta*(1-eta)));
    T[X][I] = T[S][I];
    T[Y][I] = T[S][I];
    T[ML][I] = log10((tau*delta)/(eta*(1-eta)));
    T[M][I] = log10(delta/(1-eta));

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
    T[DR][DR] = log10(epsilon/(1-eta));

    T[MR][IR] = T[MR][DR];
    T[IR][IR] = T[MR][IR];

    T[S][W] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta*eta));
    T[X][W] = T[S][W];
    T[Y][W] = T[S][W];
    T[ML][W] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta));
    T[M][W] = log10((tau*tau)/(eta*eta*eta*eta));
    T[D][W] = T[M][W];
    T[I][W] = log10((tau*tau)/(eta*eta*eta));
    T[MR][W] = log10(tau/eta);
    T[W][W] = 0;

    T[S][Z] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta*eta));
    T[X][Z] = T[S][Z];
    T[Y][Z] = T[S][Z];
    T[ML][Z] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta));
    T[M][Z] = log10((tau*tau)/(eta*eta*eta*eta));
    T[D][Z] = T[M][W];
    T[I][Z] = log10((tau*tau)/(eta*eta*eta));
    T[MR][Z] = log10(tau/eta);
    T[W][Z] = 0;
    T[Z][Z] = 0;


    V = new double*[NSTATES];
    U = new int32_t*[NSTATES];

    //the best alignment V_ for subsequence (i,j)
    for (size_t state=X; state<=Z; ++state)
    {
        V[state] = new double[MAXLEN*MAXLEN];
        U[state] = new int32_t[MAXLEN*MAXLEN];
    }

    //used for back tracking, this points to the state prior to the alignment for subsequence (i,j)
    //that ends with the corresponding state


    int32_t t=0;
    for (size_t i=0; i<MAXLEN; ++i)
    {
        for (size_t j=0; j<MAXLEN; ++j)
        {
            size_t c = i*MAXLEN+j;

            //X
            if (j) //(i,j)
            {
                V[X][c] = -DBL_MAX;
                U[X][c] = make_track(t,N,UNMODELED,0,0);
            }
            else
            {
                V[X][c] = 0;
                if (i) // (i,0)
                {
                    //t=0;
                    if (i==1)
                    {
                        t =  make_track(t,S,LFLANK,0,i);
                    }
                    else
                    {
                        t =  make_track(t,X,LFLANK,0,i);
                    }

                    U[X][c] = t;
                }
                else // (0,0)
                {
                    U[X][c] = make_track(t,N,UNMODELED,0,0);
                }
            }

            //Y
            V[Y][c] = 0;
            if (i)
            {
                if (j) // (i,j)
                {
                    U[Y][c] = j==1? make_track(t,X,UNMODELED,0,0) : make_track(t,Y,UNMODELED,0,0);
                }
                else // (i,0)
                {
                    V[Y][c] = -DBL_MAX;
                    U[Y][c] = make_track(t,N,UNMODELED,0,0);
                }
            }
            else
            {
                if (j) // (0,j)
                {
                    U[Y][c] = j==1? make_track(t,S,UNMODELED,0,0) : make_track(t,Y,UNMODELED,0,0);
                }
                else // (0,0)
                {
                    U[Y][c] = make_track(t,N,UNMODELED,0,0);
                }
            }

            //M
            V[M][c] = -DBL_MAX;
            if (!i || !j)
            {
                U[M][c] = make_track(t,N,UNMODELED,0,0);
            }
            else
            {
                U[M][c] = make_track(t,N,UNCERTAIN,0,0);
            }




//            V_I[c] = -DBL_MAX;
//            V_D[c] = -DBL_MAX;
//            V_W[c] = -DBL_MAX;
//            V_Z[c] = -DBL_MAX;
//
//            if (j)
//            {
//                U_M[c] = make_track(t,N,UNMODELED,0,0);
//                U_I[c] = make_track(t,N,UNMODELED,0,0);
//                U_D[c] = make_track(t,N,UNMODELED,0,0);
//                U_W[c] = make_track(t,N,UNMODELED,0,0);
//                U_Z[c] = make_track(t,N,UNMODELED,0,0);
//            }
//            else
//            {
//                U_M[c] = make_track(t,N,UNMODELED,0,0);
//                U_I[c] = make_track(t,N,UNMODELED,0,0);
//                U_D[c] = make_track(t,N,UNMODELED,0,0);
//                U_W[c] = make_track(t,N,UNMODELED,0,0);
//                U_Z[c] = make_track(t,N,UNMODELED,0,0);
//            }
        }
    }


    logEta = log10(eta);
    logTau = log10(tau);

    V[X][0*MAXLEN+0] = 0;
    V[Y][0*MAXLEN+0] = 0;
    V[M][0*MAXLEN+0] = 0;
    V[W][0*MAXLEN+0] = 0;
    V[Z][0*MAXLEN+0] = 0;
//    U_X[0*MAXLEN+0] = 'N';
//    U_X[1*MAXLEN+0] = 'S';
//    U_Y[0*MAXLEN+0] = 'N';
//    U_Y[0*MAXLEN+1] = 'S';
//    U_M[0*MAXLEN+0] = 'N';
//    U_M[1*MAXLEN+1] = 'S';

};


/**
 * 
 */
void CHMM::proc_comp(int32_t A, int32_t B, int32_t i, bool match)
{
    //1. check if move is valid
    
    if (V[A][i]!=-DBL_MAX && T[A][B]!=-DBL_MAX)
    {        
//        if ()
//        {
//
//        }
    }
    else
    {
        //ignore this transition
    }
}

/**
 * Align y against x.
 */
void CHMM::align(const char* read, const char* qual, bool debug)
{
    this->read = read; //read
    this->qual = qual;
    rlen = strlen(read);
    lflen = strlen(lflank);
    rflen = strlen(rflank);
    plen = lflen + rlen + rflen;

    if (rlen>MAXLEN)
    {
        fprintf(stderr, "[%s:%d %s] Sequence to be aligned is greater than %d currently supported: %d\n", __FILE__, __LINE__, __FUNCTION__, MAXLEN, rlen);
        exit(1);
    }

    if (1)
    {
        std::cerr << "rflank: " << rflank << "\n";
        std::cerr << "ru: " << ru << "\n";
        std::cerr << "lflank: " << lflank << "\n";
        std::cerr << "lflen: " << lflen << "\n";
        std::cerr << "rflen: " << rflen << "\n";
        std::cerr << "plen: " << plen << "\n";

        std::cerr << "read: " << read << "\n";
        std::cerr << "rlen: " << rlen << "\n";
    }

    double max = 0;
    char maxPath = 'X';

    size_t c,d,u,l;
    double x_ml, y_ml, ml_ml, il_ml, dl_ml;

    //alignment
    //take into consideration
    for (size_t i=1; i<=plen; ++i)
    {
        for (size_t j=1; j<=rlen; ++j)
        {
            c = i*MAXLEN+j;
            d = (i-1)*MAXLEN+(j-1);
            u = (i-1)*MAXLEN+j;
            l = i*MAXLEN+(j-1);

            ////////////////////////////
            //X matrices are invariant//
            ////////////////////////////

            ////////////////////////////
            //Y matrices are invariant//
            ////////////////////////////

            //////
            //ML//
            //////
            proc_comp(X, ML, d, true);


/**
 * Checks move from state A to B, with direction d and 
 */


//            x_ml = V_ML[d] + T[X][M];
//            y_ml = V_ML[d] + T[Y][M];
//            ml_ml = V_ML[d] + T[M][M];
//            il_ml = V_IL[d] + T[I][M];
//            dl_ml = V_DL[d] + T[D][M];


//            if (ym>max) //special case
//            {
//                max = ym;
//                maxPath = 'Y';
//            }
//            if (mm>max)
//            {
//                max = mm;
//                maxPath = (i==1&&j==1) ? 'S' : 'M';
//            }
//            if (im>max)
//            {
//                max = im;
//                maxPath = 'I';
//            }
//            if (dm>max)
//            {
//                max = dm;
//                maxPath = 'D';
//            }
//
//            V_M[i*MAXLEN+j] = max + log10_emission_odds(x[i-1], y[j-1], lt->pl2prob((uint32_t) qual[j-1]-33));
//            U_M[i*MAXLEN+j] = maxPath;

//            //D
//            double md = V_M[(i-1)*MAXLEN+j] + T[M][D];
//            double dd = V_D[(i-1)*MAXLEN+j] + T[D][D];
//
//            max = md;
//            maxPath = 'M';
//
//            if (dd>max)
//            {
//                max = dd;
//                maxPath = 'D';
//            }
//
//            V_D[i*MAXLEN+j] = max;
//            U_D[i*MAXLEN+j] = maxPath;
//
//            //I
//            double mi = V_M[i*MAXLEN+(j-1)] + T[M][I];
//            double ii = V_I[i*MAXLEN+(j-1)] + T[I][I];
//
//            max = mi;
//            maxPath = 'M';
//
//            if (ii>max)
//            {
//                max = ii;
//                maxPath = 'I';
//            }
//
//            V_I[i*MAXLEN+j] = max;
//            U_I[i*MAXLEN+j] = maxPath;
//
//            //W
//            double mw = V_M[(i-1)*MAXLEN+j] + T[M][W];
//            double ww = V_W[(i-1)*MAXLEN+j] + T[W][W];
//
//            max = mw;
//            maxPath = 'M';
//
//            if (ww>max)
//            {
//                max = ww;
//                maxPath = 'W';
//            }
//
//            V_W[i*MAXLEN+j] = max;
//            U_W[i*MAXLEN+j] = maxPath;
//
//            //Z
//            double mz = V_M[i*MAXLEN+(j-1)] + T[M][Z];
//            double wz = V_W[i*MAXLEN+(j-1)] + T[W][Z];
//            double zz = V_Z[i*MAXLEN+(j-1)] + T[Z][Z];
//
//            max = mz;
//            maxPath = 'M';
//
//            if (wz>max)
//            {
//                max = wz;
//                maxPath = 'W';
//            }
//            if (zz>max)
//            {
//                max = zz;
//                maxPath = 'Z';
//            }
//
//            V_Z[i*MAXLEN+j] = max;
//            U_Z[i*MAXLEN+j] = maxPath;
        }

//        V_M[rlen*MAXLEN+rlen] += logTau-logEta;
    }

    if (1)
    {
        std::cerr << "\n=V[X]=\n";
        print(V[X], plen+1, rlen+1);
        std::cerr << "\n=U[X]=\n";
        print_U(U[X], plen+1, rlen+1);
        std::cerr << "\n=V[Y]=\n";
        print(V[Y], plen+1, rlen+1);
        std::cerr << "\n=U[Y]=\n";
        print_U(U[Y], plen+1, rlen+1);

        std::cerr << "\n=V[ML]=\n";
        print(V[ML], plen+1, rlen+1);
        std::cerr << "\n=U[ML]=\n";
        print_U(U[ML], plen+1, rlen+1);
        std::cerr << "\n=V[DL]=\n";
        print(V[DL], plen+1, rlen+1);
        std::cerr << "\n=U[DL]=\n";
        print_U(U[DL], plen+1, rlen+1);
        std::cerr << "\n=V[IL]=\n";
        print(V[IL], plen+1, rlen+1);
        std::cerr << "\n=U[IL]=\n";
        print_U(U[IL], plen+1, rlen+1);

        std::cerr << "\n=V[M]=\n";
        print(V[M], plen+1, rlen+1);
        std::cerr << "\n=U[M]=\n";
        print_U(U[M], plen+1, rlen+1);
        std::cerr << "\n=V[D]=\n";
        print(V[D], plen+1, rlen+1);
        std::cerr << "\n=U[D]=\n";
        print_U(U[D], plen+1, rlen+1);
        std::cerr << "\n=V[I]=\n";
        print(V[I], plen+1, rlen+1);
        std::cerr << "\n=U[I]=\n";
        print_U(U[I], plen+1, rlen+1);

        std::cerr << "\n=V[MR]=\n";
        print(V[MR], plen+1, rlen+1);
        std::cerr << "\n=U[MR]=\n";
        print_U(U[MR], plen+1, rlen+1);
        std::cerr << "\n=V[DR]=\n";
        print(V[DR], plen+1, rlen+1);
        std::cerr << "\n=U[DR]=\n";
        print_U(U[DR], plen+1, rlen+1);
        std::cerr << "\n=V[IR]=\n";
        print(V[IR], plen+1, rlen+1);
        std::cerr << "\n=U[IR]=\n";
        print_U(U[IR], plen+1, rlen+1);

        std::cerr << "\n=V[W]=\n";
        print(V[W], plen+1, rlen+1);
        std::cerr << "\n=U[W]=\n";
        print_U(U[W], plen+1, rlen+1);
        std::cerr << "\n=V[Z]=\n";
        print(V[Z], plen+1, rlen+1);
        std::cerr << "\n=U[Z]=\n";
        print_U(U[Z], plen+1, rlen+1);

        std::cerr << "\n";
    }


    trace_path();
};

/**
 * Trace path after alignment.
 */
void CHMM::trace_path()
{
//    double globalMax = V[M[rlen*MAXLEN+rlen];
//    char globalMaxPath = 'M';
//    if (V_W[rlen*MAXLEN+rlen]>globalMax)
//    {
//        globalMax = V_W[rlen*MAXLEN+rlen];
//        globalMaxPath = 'W';
//    }
//    if (V_Z[rlen*MAXLEN+rlen]>globalMax)
//    {
//        globalMax = V_Z[rlen*MAXLEN+rlen];
//        globalMaxPath = 'Z';
//    }
//
//    matchStartX = -1;
//    matchEndX = -1;
//    matchStartY = -1;
//    matchEndY = -1;
//    noBasesAligned = 0;
//    if (globalMaxPath == 'M')
//    {
//        matchEndX = rlen;
//        matchEndY = rlen;
//        ++noBasesAligned;
//    }
//
//    maxLogOdds = globalMax;
//
//    ss.str("");
//    ss << globalMaxPath;
//
//    //recursively trace U_
//    trace_U_(globalMaxPath, rlen, rlen);
//
//    U_ = reverse(ss.str());
//
//    bool inI = false;
//    bool inD = false;
//
//    //records appearance of Indels in probe and read
//    indelStartsInX.clear();
//    indelEndsInX.clear();
//    indelStartsInY.clear();
//    indelEndsInY.clear();
//    indelStartsInPath.clear();
//    indelEndsInPath.clear();
//    indelStatusInPath.clear();
//
//    matchedBases = 0;
//    mismatchedBases = 0;
//
//    uint32_t x_index=0, y_index=0;
//    for (uint32_t i=0; i<U_.size(); ++i)
//    {
//        if (!inI && U_[i]=='I')
//        {
//            inI = true;
//            indelStartsInPath.push_back(i);
//            indelStatusInPath.push_back('I');
//            indelStartsInX.push_back(x_index);
//            indelStartsInY.push_back(y_index);
//        }
//
//        if (inI && U_[i]!='I')
//        {
//            inI = false;
//            indelEndsInPath.push_back(i);
//            indelEndsInX.push_back(x_index);
//            indelEndsInY.push_back(y_index);
//        }
//
//        if (!inD && U_[i]=='D')
//        {
//            inD = true;
//            indelStartsInPath.push_back(i);
//            indelStatusInPath.push_back('D');
//            indelStartsInX.push_back(x_index);
//            indelStartsInY.push_back(y_index);
//        }
//
//        if (inD && U_[i]!='D')
//        {
//            inD = false;
//            indelEndsInPath.push_back(i);
//            indelEndsInX.push_back(x_index);
//            indelEndsInY.push_back(y_index);
//        }
//
//        if (U_[i]=='M')
//        {
//            if (x[x_index]== y[y_index])
//            {
//                ++matchedBases;
//            }
//            else
//            {
//                ++mismatchedBases;
//            }
//
//            ++x_index;
//            ++y_index;
//        }
//
//        if (U_[i]=='I' || U_[i]=='Y' || U_[i]=='Z')
//        {
//            ++y_index;
//        }
//
//        if (U_[i]=='D' || U_[i]=='X' || U_[i]=='W')
//        {
//            ++x_index;
//        }
//    }
//
//    left_align();
}

/**
 * Compute log10 emission odds based on equal error probability distribution.
 */
double CHMM::log10_emission_odds(char readBase, char probeBase, uint32_t pl)
{
    //4 encodes for N
    if (readBase=='N' || probeBase=='N')
    {
        //silent match
        return -DBL_MAX;
    }

    if (readBase!=probeBase)
    {
        return lt->pl2log10_ed3(pl);
    }
    else
    {
        return lt->pl2log10_1me(pl);
    }
};

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
//    std::stringstream xAligned;
//    std::stringstream yAligned;
//    std::stringstream qualAligned;
//    uint32_t xIndex=0, yIndex=0;
//    for (uint32_t i=0; i<U_.size(); ++i)
//    {
//        char state = U_.at(i);
//        if (state=='S' || state=='E')
//        {
//            xAligned << state;
//            yAligned << state;
//            qualAligned << state;
//        }
//        else if (state=='X'||state=='W')
//        {
//            xAligned << x[xIndex++];
//            yAligned << '.';
//            qualAligned << '.';
//        }
//        else if (state=='M')
//        {
//            xAligned << x[xIndex++];
//            yAligned << y[yIndex];
//            qualAligned << qual[yIndex++];
//        }
//        else if (state=='D')
//        {
//            xAligned << x[xIndex++];
//            yAligned << '-';
//            qualAligned << '-';
//        }
//        else if (state=='I')
//        {
//            xAligned << '-';
//            yAligned << y[yIndex];
//            qualAligned << qual[yIndex++];
//        }
//        else if (state=='Y'||state=='Z')
//        {
//            xAligned << '.';
//            yAligned << y[yIndex];
//            qualAligned << qual[yIndex++];
//        }
//    }
//
//    std::cerr << pad << "X    " << xAligned.str() << "E\n";
//    std::cerr << pad << "Path " << U_ << "E\n";
//    std::cerr << pad << "Y    " << yAligned.str() << "E\n";
//    std::cerr << pad << "Qual " << qualAligned.str() << "E\n\n";
//
//    std::cerr << pad << "Alignment in Probe : [" << matchStartX << "," << matchEndX<< "]\n";
//    std::cerr << pad << "Alignment in Read  : [" << matchStartY << "," << matchEndY<< "]\n";
//    std::cerr << pad << "# Aligned Bases : " << noBasesAligned << "\n";
//    std::cerr << pad << "# Matched Bases : " << matchedBases << "\n";
//    std::cerr << pad << "# Mismatched Bases : " << mismatchedBases << "\n";
//    std::cerr << pad << "Max Log odds    : " << maxLogOdds << "\n";
};

/**
 * Prints a double matrix.
 */
void CHMM::print(double *v, size_t plen, size_t rlen)
{
    for (size_t i=0; i<plen; ++i)
    {
        for (size_t j=0; j<rlen; ++j)
        {
            std::cerr << (v[i*MAXLEN+j]==-DBL_MAX?NAN:v[i*MAXLEN+j]) << "\t";
        }

        std::cerr << "\n";
    }
};

/**
 * Prints a char matrix.
 */
void CHMM::print(int32_t *v, size_t plen, size_t rlen)
{
    for (size_t i=0; i<plen; ++i)
    {
        for (size_t j=0; j<rlen; ++j)
        {
          std::cerr << v[i*MAXLEN+j] << "\t";
        }

        std::cerr << "\n";
    }
};

/**
 * Prints U.
 */
void CHMM::print_U(int32_t *U, size_t plen, size_t rlen)
{
    for (size_t i=0; i<plen; ++i)
    {
        for (size_t j=0; j<rlen; ++j)
        {
            int32_t t = U[i*MAXLEN+j];
            std::cerr << state2string(track_get_u(t)) << "|"
                      << component2string(track_get_d(t)) << "|"
                      << track_get_c(t) << "|"
                      << track_get_p(t) << (j==rlen-1?"\n":"   ");
        }
    }
};


/**
 * Prints track.
 */
void CHMM::print_track(int32_t t)
{
    std::cerr << state2string(track_get_u(t)) << "|"
          << component2string(track_get_d(t)) << "|"
          << track_get_c(t) << "|"
          << track_get_p(t) << "\n";
}

#define track_get_u(t) (((t)&0xFF000000)>>24)
#define track_get_d(t) (((t)&0x00FF0000)>>16)
#define track_get_c(t) (((t)&0x0000FF00)>>8)
#define track_get_p(t) (((t)&0x000000FF))

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

#undef LFLANK
#undef MOTIF
#undef RFLANK