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
#define LM 3
#define LI 4
#define LD 5
#define M  6
#define I  7
#define D  8
#define RM 9
#define RI 10
#define RD 11
#define W  12
#define Z  13
#define E  14

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
    delete V_X;
    delete V_Y;
    delete V_LM;
    delete V_LI;
    delete V_LD;
    delete V_M;
    delete V_I;
    delete V_D;
    delete V_RM;
    delete V_RI;
    delete V_RD;
    delete V_W;
    delete V_Z;

    delete positionLM;
    delete positionLI;
    delete positionLD;
    delete positionM;
    delete positionI;
    delete positionD;
    delete positionRM;
    delete positionRI;
    delete positionRD;

    delete U_X;
    delete U_Y;
    delete U_LM;
    delete U_LI;
    delete U_LD;
    delete U_M;
    delete U_I;
    delete U_D;
    delete U_RM;
    delete U_RI;
    delete U_RD;
    delete U_W;
    delete U_Z;
};

/**
 * Initializes object, helper function for constructor.
 */
void CHMM::initialize(const char* lflank, const char* ru, const char* rflank)
{
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

    T[S][LM] = log10((1-tau)/(eta*(1-eta)*(1-eta)));
    T[X][LM] = T[S][LM];
    T[Y][LM] = T[S][LM];
    T[LM][LM] = log10(((1-2*delta-tau))/((1-eta)*(1-eta)));
    T[LD][LM] = log10(((1-epsilon))/((1-eta)*(1-eta)));
    T[LI][LM] = T[LD][LM];

    T[LM][LD] = log10((delta)/((1-eta)));
    T[LD][LD] = log10((epsilon)/((1-eta)));

    T[LM][LI] = T[LM][LD];
    T[LI][LI] = T[LD][LD];

    T[S][M] = log10((tau*(1-2*delta-tau))/(eta*eta*eta*(1-eta)*(1-eta)));
    T[X][M] = T[S][M];
    T[Y][M] = T[S][M];
    T[LM][M] = log10((tau*(1-2*delta-tau))/(eta*eta*eta*(1-eta)*(1-eta)));
    T[M][M] = log10((tau*(1-2*delta-tau))/(eta*eta*(1-eta)*(1-eta)));
    T[D][M] = log10(((1-epsilon-tau))/((1-eta)*(1-eta)));
    T[I][M] = T[I][M];

    T[S][D] = log10((tau*delta)/(eta*eta*(1-eta)));
    T[X][D] = T[S][D];
    T[Y][D] = T[S][D];
    T[LM][D] = log10((tau*delta)/((1-eta)));
    T[M][D] = log10(delta/(1-eta));

    T[S][I] = log10((tau*delta)/(eta*eta*eta*(1-eta)));
    T[X][I] = T[S][I];
    T[Y][I] = T[S][I];
    T[LM][I] = log10((tau*delta)/(eta*(1-eta)));
    T[M][I] = log10(delta/(1-eta));

    T[S][RM] = log10((tau*tau*(1-tau))/(eta*eta*eta*eta*eta*(1-eta)*(1-eta)));
    T[X][RM] = T[S][RM];
    T[Y][RM] = T[S][RM];
    T[LM][RM] = log10((tau*tau*(1-tau))/(eta*eta*eta*eta*(1-eta)*(1-eta)));
    T[M][RM] = log10((tau*(1-tau))/(eta*eta*eta*(1-eta)*(1-eta)));
    T[D][RM] = T[M][RM];
    T[I][RM] = log10((tau*(1-tau))/(eta*eta*(1-eta)*(1-eta)));
    T[RM][RM] = log10(((1-2*delta-tau))/((1-eta)*(1-eta)));
    T[RD][RM] = log10(((1-epsilon))/((1-eta)*(1-eta)));
    T[RI][RM] = T[RD][RM];

    T[RM][RD] = log10(delta/(1-eta));
    T[RD][RD] = log10(epsilon/(1-eta));

    T[RM][RI] = T[RM][RD];
    T[RI][RI] = T[RM][RI];

    T[S][W] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta*eta));
    T[X][W] = T[S][W];
    T[Y][W] = T[S][W];
    T[LM][W] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta));
    T[M][W] = log10((tau*tau)/(eta*eta*eta*eta));
    T[D][W] = T[M][W];
    T[I][W] = log10((tau*tau)/(eta*eta*eta));
    T[RM][W] = log10(tau/eta);
    T[W][W] = 0;

    T[S][Z] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta*eta));
    T[X][Z] = T[S][Z];
    T[Y][Z] = T[S][Z];
    T[LM][Z] = log10((tau*tau*tau)/(eta*eta*eta*eta*eta));
    T[M][Z] = log10((tau*tau)/(eta*eta*eta*eta));
    T[D][Z] = T[M][W];
    T[I][Z] = log10((tau*tau)/(eta*eta*eta));
    T[RM][Z] = log10(tau/eta);
    T[W][Z] = 0;
    T[Z][Z] = 0;

    //the best alignment V_ for subsequence (i,j)
    V_X  = new double[MAXLEN*MAXLEN];
    V_Y  = new double[MAXLEN*MAXLEN];
    V_LM = new double[MAXLEN*MAXLEN];
    V_LI = new double[MAXLEN*MAXLEN];
    V_LD = new double[MAXLEN*MAXLEN];
    V_M  = new double[MAXLEN*MAXLEN];
    V_I  = new double[MAXLEN*MAXLEN];
    V_D  = new double[MAXLEN*MAXLEN];
    V_RM = new double[MAXLEN*MAXLEN];
    V_RI = new double[MAXLEN*MAXLEN];
    V_RD = new double[MAXLEN*MAXLEN];
    V_W  = new double[MAXLEN*MAXLEN];
    V_Z  = new double[MAXLEN*MAXLEN];

    //the matching position of the model sequence for the best alignment for subsequence (i,j)
    positionLM = new int32_t[MAXLEN*MAXLEN];
    positionLI = new int32_t[MAXLEN*MAXLEN];
    positionLD = new int32_t[MAXLEN*MAXLEN];
    positionM  = new int32_t[MAXLEN*MAXLEN];
    positionI  = new int32_t[MAXLEN*MAXLEN];
    positionD  = new int32_t[MAXLEN*MAXLEN];
    positionRM = new int32_t[MAXLEN*MAXLEN];
    positionRI = new int32_t[MAXLEN*MAXLEN];
    positionRD = new int32_t[MAXLEN*MAXLEN];

    //used for back tracking, this points to the state prior to the alignment for subsequence (i,j)
    //that ends with the corresponding state
    U_X  = new char[MAXLEN*MAXLEN];
    U_Y  = new char[MAXLEN*MAXLEN];
    U_LM = new char[MAXLEN*MAXLEN];
    U_LI = new char[MAXLEN*MAXLEN];
    U_LD = new char[MAXLEN*MAXLEN];
    U_M  = new char[MAXLEN*MAXLEN];
    U_I  = new char[MAXLEN*MAXLEN];
    U_D  = new char[MAXLEN*MAXLEN];
    U_RM = new char[MAXLEN*MAXLEN];
    U_RI = new char[MAXLEN*MAXLEN];
    U_RD = new char[MAXLEN*MAXLEN];
    U_W  = new char[MAXLEN*MAXLEN];
    U_Z  = new char[MAXLEN*MAXLEN];

    for (size_t i=0; i<MAXLEN; ++i)
    {
        for (size_t j=0; j<MAXLEN; ++j)
        {
            size_t c = i*MAXLEN+j;
           
            //X
            //note that for j>0, the values are invalid in X since Y can never preceed X
            if (j) //(i,j)
            {
                V_X[c] = 0;
                U_X[c] = 'N';
            }
            else
            {
                V_X[c] = -DBL_MAX;
                if (i) // (i,0)
                {
                    U_X[c] = i==1? 'S' : 'X';
                }  
                else // (0,0)
                {
                    U_X[c] = 'N';   
                } 
            }    
           
            //Y    
            V_Y[c] = 0;
            if (i)
            {
                if (j) // (i,j)
                {
                    U_Y[c] = j==1? 'X' : 'Y';
                }  
                else // (i,0)
                {
                    U_Y[c] = 'N';   
                } 
            }
            else
            {
                if (j) // (0,j)
                {
                    U_Y[c] = j==1? 'S' : 'Y';
                }  
                else // (0,0)
                {
                    U_Y[c] = 'N';   
                }
            } 
            
            V_M[c] = -DBL_MAX;
            V_I[c] = -DBL_MAX;
            V_D[c] = -DBL_MAX;
            V_W[c] = -DBL_MAX;
            V_Z[c] = -DBL_MAX;

            if (j)
            {
                U_X[c] = 'Y';
                U_Y[c] = 'Y';
                U_M[c] = 'Y';
                U_I[c] = 'Y';
                U_D[c] = 'Y';
                U_W[c] = 'Y';
                U_Z[c] = 'Y';
            }
            else
            {
                U_X[c] = 'X';
                U_Y[c] = 'X';
                U_M[c] = 'X';
                U_I[c] = 'X';
                U_D[c] = 'X';
                U_W[c] = 'X';
                U_Z[c] = 'X';
            }
        }
    }
    

    logEta = log10(eta);
    logTau = log10(tau);

    V_X[0*MAXLEN+0] = 0;
    V_Y[0*MAXLEN+0] = 0;
    V_M[0*MAXLEN+0] = 0;
    V_W[0*MAXLEN+0] = 0;
    V_Z[0*MAXLEN+0] = 0;
    U_X[0*MAXLEN+0] = 'N';
    U_X[1*MAXLEN+0] = 'S';
    U_Y[0*MAXLEN+0] = 'N';
    U_Y[0*MAXLEN+1] = 'S';
    U_M[0*MAXLEN+0] = 'N';
    U_M[1*MAXLEN+1] = 'S';

};

/**
 * Align y against x.
 */
void CHMM::align(const char* read, const char* qual, bool debug)
{
    this->read = read; //read
    this->qual = qual;
    rlen = strlen(read);
    plen = lflen + rlen + rflen;

    if (rlen>MAXLEN)
    {
        fprintf(stderr, "[%s:%d %s] Sequence to be aligned is greater than %d currently supported: %d\n", __FILE__, __LINE__, __FUNCTION__, MAXLEN, rlen);
        exit(1);
    }

    double max = 0;
    char maxPath = 'X';

    size_t c,d,u,l;

    //alignment
    //take into consideration
    for (size_t i=1; i<=rlen; ++i)
    {
        for (size_t j=1; j<=rlen; ++j)
        {
            c = i*MAXLEN+j;
            d = (i-1)*MAXLEN+(j-1);
            u = (i-1)*MAXLEN+j;
            l = i*MAXLEN+(j-1);

            //X matrices are invariant

            //Y matrices are invariant

            ////////////////
            //LM
            ////////////////
            
            //pick up base
            double xlm = T[X][M];
            
            double ylm = T[Y][M];
            double lmlm = V_LM[d] + T[M][M];
            double lilm = V_LI[d] + T[I][M];
            double ldlm = V_LD[d] + T[D][M];

            max = xlm;
            maxPath = 'X';
//
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
//           // V_M[i*MAXLEN+j] = max + log10_emission_odds(x[i-1], y[j-1], lt->pl2prob((uint32_t) qual[j-1]-33));
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

        V_M[rlen*MAXLEN+rlen] += logTau-logEta;
    }

    if (debug)
    {
        std::cerr << "\n=X=\n";
        print(V_X, rlen+1);
        std::cerr << "\n=Y=\n";
        print(V_Y, rlen+1);
        std::cerr << "\n=M=\n";
        print(V_M, rlen+1);
        std::cerr << "\n=D=\n";
        print(V_D, rlen+1);
        std::cerr << "\n=I=\n";
        print(V_I, rlen+1);
        std::cerr << "\n=W=\n";
        print(V_W, rlen+1);
        std::cerr << "\n=Z=\n";
        print(V_Z, rlen+1);
        std::cerr << "\n=Path X=\n";
        print(U_X, rlen+1);
        std::cerr << "\n=Path Y=\n";
        print(U_Y, rlen+1);
        std::cerr << "\n=Path M=\n";
        print(U_M, rlen+1);
        std::cerr << "\n=Path D=\n";
        print(U_D, rlen+1);
        std::cerr << "\n=Path I=\n";
        print(U_I, rlen+1);
        std::cerr << "\n=Path W=\n";
        print(U_W, rlen+1);
        std::cerr << "\n=Path Z=\n";
        print(U_Z, rlen+1);
    }

    trace_path();
};

/**
 * Trace U_ after alignment.
 */
void CHMM::trace_path()
{
//    double globalMax = V_M[rlen*MAXLEN+rlen];
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
double CHMM::log10_emission_odds(char readBase, char probeBase, double e)
{
    //4 encodes for N
    if (readBase=='N' || probeBase=='N')
    {
        //silent match
        return 0;
    }

    if (readBase!=probeBase)
    {
        return log10(e/3)-logOneSixteenth;
    }
    else
    {
        return log10(1-e)-logOneSixteenth;
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
void CHMM::print(double *v, uint32_t rlen)
{
    for (uint32_t i=0; i<rlen; ++i)
    {
        for (uint32_t j=0; j<rlen; ++j)
        {
            std::cerr << (v[i*MAXLEN+j]==-DBL_MAX?-1000:v[i*MAXLEN+j]) << "\t";
        }

        std::cerr << "\n";
    }
};

/**
 * Prints a char matrix.
 */
void CHMM::print(char *v, uint32_t rlen)
{
    for (uint32_t i=0; i<rlen; ++i)
    {
        for (uint32_t j=0; j<rlen; ++j)
        {
          std::cerr << v[i*MAXLEN+j] << "\t";
        }

        std::cerr << "\n";
    }
};

#undef MAXLEN
#undef S
#undef X
#undef Y
#undef LM
#undef LI
#undef LD
#undef M
#undef I
#undef D
#undef RM
#undef RI
#undef RD
#undef W
#undef Z
#undef E