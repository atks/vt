/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#include "lhmm1.h"

#include "lhmm1.h"

LHMM1::LHMM1()
{
    initialize(500);
};

void LHMM1::initialize(int32_t maxLength)
{
    this->maxLength = maxLength;
//  delta = 0.01;
//  epsilon = 0.1;
//  tau = 0.1;
//  eta = 0.001;

//  delta = 0.01;
//  epsilon = 0.1;
//  tau = 0.1;
//  eta = 0.001;

//  last
//  delta = 0.001;
//  epsilon = 0.5;
//  tau = 0.1;
//  eta = 0.01;

//  delta = 0.01;
//  epsilon = 0.5;
//  tau = 0.1;
//  eta = 0.001;

    delta = 0.001;
    epsilon = 0.5;
    tau = 0.1;
    eta = 0.001;

    logOneSixteenth = log10(1.0/16.0);

    //for matched portion
    tMM = log10(1-2*delta-tau);
    tMI = log10(delta);
    tMD = log10(delta);
    tII = log10(epsilon);
    tDD = log10(epsilon);
    tDM = log10(1-epsilon);
    tIM = log10(1-epsilon);

    txx = log10((1-eta)/(1-eta));
    txy = log10((eta*(1-eta))/(eta*(1-eta)));
    tyy = log10((1-eta)/(1-eta));

    txm = log10((eta*eta)/(eta*eta*(1-eta)*(1-eta)));
    tym = log10((eta)/(eta*(1-eta)*(1-eta)));
    tsm = log10((eta*eta)/(eta*eta*eta*(1-eta)*(1-eta)));
    tmm = log10((1-2*delta-tau)/((1-eta)*(1-eta)));
    tim = log10((1-epsilon)/(1-eta));
    tdm = log10((1-epsilon)/(1-eta));

    tmd = log10(delta/(1-eta));
    tdd = log10(epsilon/(1-eta));

    tmi = log10(delta/(1-eta));
    tii = log10(epsilon/(1-eta));

    tmw = log10((tau*(1-eta))/(eta*(1-eta)));
    tww = log10((1-eta)/(1-eta));

    tmz = log10((tau*eta*(1-eta))/(eta*eta*(1-eta)));
    twz = log10((eta*(1-eta))/(eta*(1-eta)));
    tzz = log10((1-eta)/(1-eta));

    //assume alignments can't possibly be maxLength bases or more
    //todo:: allow for resizing if you do encounter fragments which are really long
    std::vector<double> temp(maxLength+1, -DBL_MAX);
    std::vector<char> tempX(maxLength+1, 'Y');
    tempX[0] = 'X';
    std::vector<char> tempY(maxLength+1, 'Y');
    tempY[0] = 'X';

    X.clear();
    Y.clear();
    M.clear();
    I.clear();
    D.clear();
    W.clear();
    Z.clear();
    pathX.clear();
    pathY.clear();
    pathM.clear();
    pathD.clear();
    pathI.clear();
    pathW.clear();
    pathZ.clear();

    for (uint32_t i=0; i<maxLength+1; ++i)
    {
        X.push_back(temp);
        Y.push_back(temp);
        M.push_back(temp);
        I.push_back(temp);
        D.push_back(temp);
        W.push_back(temp);
        Z.push_back(temp);
        pathX.push_back(tempX);
        pathY.push_back(tempY);
        pathM.push_back(tempX);
        pathD.push_back(tempX);
        pathI.push_back(tempX);
        pathW.push_back(tempX);
        pathZ.push_back(tempX);
    }

    logEta = log10(eta);
    logTau = log10(tau);

    X[0][0] = 0;
    Y[0][0] = 0;
    M[0][0] = 0;
    W[0][0] = 0;
    Z[0][0] = 0;
    pathX[0][0] = 'N';
    pathX[1][0] = 'S';
    pathY[0][0] = 'N';
    pathY[0][1] = 'S';
    pathM[0][0] = 'N';
    pathM[1][1] = 'S';

    for (uint32_t k=1; k<maxLength+1; ++k)
    {
        X[k][0] = X[k-1][0] + txx;
        X[0][k] = -DBL_MAX;
        Y[k][0] = -DBL_MAX;
        Y[0][k] = Y[0][k-1] + tyy;
        W[k][0] = W[k-1][0] + tww;
        W[0][k] = -DBL_MAX;
        Z[k][0] = -DBL_MAX;
        Z[0][k] = Z[0][k-1] + tzz;
    }

    X[0][0] = -DBL_MAX;
    Y[0][0] = -DBL_MAX;
    W[0][0] = -DBL_MAX;
    Z[0][0] = -DBL_MAX;

};

bool LHMM1::containsIndel()
{
    return (indelStatusInPath.size()!=0);
}

/**
convert PLs to probabilities.
*/
double LHMM1::pl2prob(uint32_t PL)
{
    if (PL>=PLs.size())
    {
        if (PL > 3236)
            PL = 3236;

        for (uint32_t i=PLs.size(); i<=PL; ++i)
        {
            PLs.push_back(pow(10, -((double) i)/10.0));
        }
    }

    return PLs[PL];
}

/**
 *Updates matchStart, matchEnd, globalMaxPath and path.
 */
void LHMM1::tracePath()
{
    double globalMax = M[xlen][ylen];
    char globalMaxPath = 'M';
    if (W[xlen][ylen]>globalMax)
    {
        globalMax = W[xlen][ylen];
        globalMaxPath = 'W';
    }
    if (Z[xlen][ylen]>globalMax)
    {
        globalMax = Z[xlen][ylen];
        globalMaxPath = 'Z';
    }

    matchStartX = -1;
    matchEndX = -1;
    matchStartY = -1;
    matchEndY = -1;
    noBasesAligned = 0;
    if (globalMaxPath == 'M')
    {
        matchEndX = xlen;
        matchEndY = ylen;
        ++noBasesAligned;
    }

    maxLogOdds = globalMax;

    //std::cerr << "\tGlobal max Path : " << globalMaxPath << "\n";
    //std::cerr << "\tGlobal max log odds : " << globalMax << "\n";

    std::stringstream ss;
    ss << globalMaxPath;

    //recursively trace path
    tracePath(ss, globalMaxPath, xlen, ylen);

    path = reverse(ss.str());

    bool inI = false;
    bool inD = false;
    //std::cerr << "*************************************\n";
    //std::cerr << "PATH :" << path << "\n";

    ////////////////////////////////////////////////
    //records appearance of Indels in probe and read
    ////////////////////////////////////////////////
    indelStartsInX.clear();
    indelEndsInX.clear();
    indelStartsInY.clear();
    indelEndsInY.clear();
    indelStartsInPath.clear();
    indelEndsInPath.clear();
    indelStatusInPath.clear();

    uint32_t x=0, y=0;
    for (uint32_t i=0; i<path.size(); ++i)
    {
        if (!inI && path[i]=='I')
        {
            inI = true;
            indelStartsInPath.push_back(i);
            indelStatusInPath.push_back('I');
            indelStartsInX.push_back(x);
            indelStartsInY.push_back(y);
        }

        if (inI && path[i]!='I')
        {
            inI = false;
            indelEndsInPath.push_back(i);
            indelEndsInX.push_back(x);
            indelEndsInY.push_back(y);
        }

        if (!inD && path[i]=='D')
        {
            inD = true;
            indelStartsInPath.push_back(i);
            indelStatusInPath.push_back('D');
            indelStartsInX.push_back(x);
            indelStartsInY.push_back(y);
        }

        if (inD && path[i]!='D')
        {
            inD = false;
            indelEndsInPath.push_back(i);
            indelEndsInX.push_back(x);
            indelEndsInY.push_back(y);
        }

        if (path[i]=='M')
        {
            ++x;
            ++y;
        }

        if (path[i]=='I' || path[i]=='Y' || path[i]=='Z')
        {
            ++y;
        }

        if (path[i]=='D' || path[i]=='X' || path[i]=='W')
        {
            ++x;
        }
    }

//    std::cerr << "\t# INS in Path : ";
//    for (uint32_t i=0; i<insertionStartsInPath.size(); ++i)
//    {
//         std::cerr << "(" << insertionStartsInPath[i] << "," << insertionEndsInPath[i] << ") ";
//    }
//    std::cerr << "\n";
//
//    std::cerr << "\t# DEL in Path : ";
//    for (uint32_t i=0; i<deletionStartsInPath.size(); ++i)
//    {
//         std::cerr << "(" << deletionStartsInPath[i] << "," << deletionEndsInPath[i] << ") ";
//    }
//    std::cerr << "\n";
//    std::cerr << "*************************************\n";
    //left align path
    left_align();

//    std::cerr << "*************************************\n";
//    printAlignment();
//    std::cerr << "*************************************\n";

    // std::cerr << "\tPath: " << path << "\n";
}

void LHMM1::tracePath(std::stringstream& ss, char state, uint32_t i, uint32_t j)
{
    //std::cout << state << "\n";
    if (i>0 && j>0)
    {
        if (state=='X')
        {
            //std::cerr << pathX[i][j] << " " << i << " " << j << " " << matchStartY << "\n";
            ss << pathX[i][j];
            tracePath(ss, pathX[i][j], i-1, j);
        }
        else if (state=='Y')
        {
            //std::cerr << pathY[i][j] << " " << i << " " << j << " " << matchStartY << "\n";
            ss << pathY[i][j];
            tracePath(ss, pathY[i][j], i, j-1);
        }
        else if (state=='M')
        {
//            if (matchEndX==-1)
//            {
//                matchEndX = i;
//                matchEndY = j;
//
//                //std::cerr << "set matchEndX "  << matchEndX  << " " << matchEndY  <<  "\n";
//            }
//
            if (matchStartX==-1 && (pathM[i][j] =='X' || pathM[i][j]=='Y'))
            {
               matchStartX = i;
               matchStartY = j;

               //std::cerr << "set match starts "  << matchStartX  << " " << matchStartY  <<  "\n";
            }


            //std::cout << pathM[i][j] << " " << i << " " << j << " " << matchEndY << "\n";
            ss << pathM[i][j];
            tracePath(ss, pathM[i][j], i-1, j-1);
            ++noBasesAligned;
        }
        else if (state=='I')
        {
            //std::cout << pathI[i][j] << " " << i << " " << j << "\n";
            ss << pathI[i][j];
            tracePath(ss, pathI[i][j], i, j-1);
        }
        else if (state=='D')
        {
            //std::cout << pathD[i][j] << " " << i << " " << j << "\n";
            ss << pathD[i][j];
            tracePath(ss, pathD[i][j], i-1, j);
            ++noBasesAligned;
        }
        else if (state=='W')
        {
            if (matchEndX==-1 && pathW[i][j] =='M')
            {
                matchEndX = i-1;
                matchEndY = j;
            }


            //std::cout << pathW[i][j] << " " << i << " " << j << "\n";
            ss << pathW[i][j];
            tracePath(ss, pathW[i][j], i-1, j);
        }
        else if (state=='Z')
        {
            if (matchEndX==-1 && pathZ[i][j] =='M')
            {
                matchEndX = i;
                matchEndY = j-1;
            }

            //std::cout << pathZ[i][j] << " " << i << " " << j << "\n";
            ss << pathZ[i][j];
            tracePath(ss, pathZ[i][j], i, j-1);
        }
        else if (state=='S')
        {
            //std::cout << "\n";
            if (matchStartX==-1)
            {
               matchStartX = i+1;
            }

            if (matchStartY==-1)
            {
               matchStartY = j+1;
            }

        }
    }
    else if (i==0 && j>0)
    {
        if (matchStartY==-1)
        {
           matchStartY = j+1;
        }

        // std::cout << pathY[i][j] << " " << i << " " << j << "\n";
        ss << pathY[i][j];
        tracePath(ss, pathY[i][j], i, j-1);
    }
    else if (i>0 && j==0)
    {
        if (matchStartY==-1)
        {
           matchStartY = j+1;
        }
        // std::cout << pathX[i][j] << " " << i << " " << j << "\n";
        ss << pathX[i][j];
        tracePath(ss, pathX[i][j], i-1, j);
    }
    else
    {
        //std::cout << "\n";
        //ss << pathX[i][j];
    }
}

/**
 *Left align indels in an alignment
 */
void LHMM1::left_align()
{
//.......... GATAAGTGAG GAGGAAAAGC GA---GGAGA TCTTATTTGACAACTGCE
//YYYYYYYYYY MMMMMMMMMM MMMMMMMMMM MMIIIMMMMM MMMMMWWWWWWWWWWWWE
//ACTGCTTCTA GATAAGTGAG GAGGAAAAGC GATAAGGAGA TCTTA............E
//# INS in Probe : (22,22)
//# INS in Read : (32,35)

    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='I')
        {
            while (indelStartsInX[i]>1 && indelStartsInY[i]>1)
            {
//                printAlignment();
//
//                std::cerr << "X[" << indelStartsInX[i]-1 << "] "  << x[indelStartsInX[i]-1] << "\n";
//                std::cerr << "Y[" << indelStartsInY[i]-1 << "] "  << y[indelStartsInY[i]-1] << "\n";
//                std::cerr << "Y[" << indelEndsInY[i]-1 << "] "  << y[indelEndsInY[i]-1] << "\n";
//
                if (//(indelStartsInY[i]!=indelEndsInY[i]-1) && //single indel need not be shifted
                    x[indelStartsInX[i]-1]==y[indelStartsInY[i]-1] &&
                    y[indelStartsInY[i]-1]==y[indelEndsInY[i]-1])
                {
//                    std::cerr << "CHECK BEFORE:" << path << "\n";
//                    std::cerr << "P[" << indelStartsInPath[i]-1 << "] "  << path[indelStartsInPath[i]-1] << "\n";
//                    std::cerr << "P[" << indelEndsInPath[i]-1 << "] "  << path[indelEndsInPath[i]-1] << "\n";

                    path[indelStartsInPath[i]-1] = 'I';
                    path[indelEndsInPath[i]-1] = 'M';

                    --indelStartsInX[i];
                    --indelEndsInX[i];
                    --indelStartsInY[i];
                    --indelEndsInY[i];
                    --indelStartsInPath[i];
                    --indelEndsInPath[i];

//                    std::cerr << "SHIFT 1 to the left\n";

                }
                else
                {
                    break;
                }
            }
        }
        else
        {
            //  ==================
            //  2) 6:81319178:AATTT:A:-4
            //  ==================
            //  ref probe     ATCTATTAAAGATATATGTCAATTTATTAGTGCATATCAAAGAGCAAGT (20/49)
            //  read sequence GATATATGTCAATTAGTGCATATCAAAGAGCAAGTTGACATGTTTTTTCAATATATTCATTTCTCTAATTTATCTC
            //  ==================
            //  X    SATCTATTAAAGATATATGTCAATTTATTAGTGCATATCAAAGAGCAAGT.........................................E
            //  Path SXXXXXXXXXXMMMMMMMMMMMMMDDDDMMMMMMMMMMMMMMMMMMMMMMZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZE
            //  Y    S..........GATATATGTCAAT----TAGTGCATATCAAAGAGCAAGTTGACATGTTTTTTCAATATATTCATTTCTCTAATTTATCTCE
            //  Match Start     : 0
            //  Match End       : 35
            //  # Aligned Bases : 39
            //  # Matched Bases : 35
            //  # Mismatched Bases : 0
            //  # INS in Probe :
            //  # INS in Read :
            //  # INS in Path :
            //  # DEL in Probe : (23,27)
            //  # DEL in Read : (13,13)
            //  # DEL in Path : (24,28)

            while (indelStartsInX[i]>1 && indelStartsInY[i]>1)
            {
//                printAlignment();
//                std::cerr << "X[" << indelStartsInX[i]-1 << "] "  << x[indelStartsInX[i]-1] << "\n";
//                std::cerr << "Y[" << indelStartsInY[i]-1 << "] "  << y[indelStartsInY[i]-1] << "\n";
//                std::cerr << "Y[" << indelEndsInY[i]-1 << "] "  << y[indelEndsInY[i]-1] << "\n";

                if (//(indelStartsInX[i]!=indelEndsInX[i]-1) && //single indel need not be shifted
                    x[indelStartsInX[i]-1]==y[indelStartsInY[i]-1] &&
                    x[indelStartsInX[i]-1]==x[indelEndsInX[i]-1])
                {
//                    std::cerr << "CHECK BEFORE:" << path << "\n";
//                    std::cerr << "P[" << indelStartsInPath[i]-1 << "] "  << path[indelStartsInPath[i]-1] << "\n";
//                    std::cerr << "P[" << indelEndsInPath[i]-1 << "] "  << path[indelEndsInPath[i]-1] << "\n";

                    path[indelStartsInPath[i]-1] = 'D';
                    path[indelEndsInPath[i]-1] = 'M';

                    --indelStartsInX[i];
                    --indelEndsInX[i];
                    --indelStartsInY[i];
                    --indelEndsInY[i];
                    --indelStartsInPath[i];
                    --indelEndsInPath[i];

//                    std::cerr << "SHIFT 1 to the left\n";

                }
                else
                {
                    break;
                }
            }

        }
    }
}

/**
Align and compute genotype likelihood.
*/
void LHMM1::align(double& llk, const char* _x, const char* _y, const char* _qual, bool debug)
{
    x = _x;
    y = _y;
    qual = _qual;

    //adds a starting character at the fron of each string that must be matched
    xlen = strlen(x);
    ylen = strlen(y);
    
    int32_t mlen = xlen>ylen?xlen:ylen;
    if (mlen>this->maxLength)
    {    
        initialize(mlen+10);
    }
    
    double max = 0;
    char maxPath = 'X';

    //std::cerr << "x: " << x << " " << xlen << "\n" ;
    //std::cerr << "y: " << y << " " << ylen << "\n" ;

    //construct possible solutions
    for (uint32_t i=1; i<=xlen; ++i)
    {
        for (uint32_t j=1; j<=ylen; ++j)
        {
            //std::cerr << i << " " << j  << "\n" ;

            //X
            double xx = X[i-1][j] + txx;

            max = xx;
            maxPath = 'X';

            X[i][j] = max;
            pathX[i][j] = maxPath;

            //Y
            double xy = X[i][j-1] + txy;
            double yy = Y[i][j-1] + tyy;

            max = xy;
            maxPath = 'X';

            if (yy>max)
            {
                max = yy;
                maxPath = 'Y';
            }

            Y[i][j] = max;
            pathY[i][j] = maxPath;

            //M
            double xm = X[i-1][j-1] + txm;
            double ym = Y[i-1][j-1] + tym;
            double mm = M[i-1][j-1] + ((i==1&&j==1) ? tsm : tmm);
            double im = I[i-1][j-1] + tim;
            double dm = D[i-1][j-1] + tdm;

            max = xm;
            maxPath = 'X';

            if (ym>max) //special case
            {
                max = ym;
                maxPath = 'Y';
            }
            if (mm>max)
            {
                max = mm;
                maxPath = (i==1&&j==1) ? 'S' : 'M';
            }
            if (im>max)
            {
                max = im;
                maxPath = 'I';
            }
            if (dm>max)
            {
                max = dm;
                maxPath = 'D';
            }

            M[i][j] = max + logEmissionOdds(x[i-1], y[j-1], pl2prob((uint32_t) qual[j-1]-33));
            pathM[i][j] = maxPath;

            //D
            double md = M[i-1][j] + tmd;
            double dd = D[i-1][j] + tdd;

            max = md;
            maxPath = 'M';

            if (dd>max)
            {
                max = dd;
                maxPath = 'D';
            }

            D[i][j] = max;
            pathD[i][j] = maxPath;

            //I
            double mi = M[i][j-1] + tmi;
            double ii = I[i][j-1] + tii;

            max = mi;
            maxPath = 'M';

            if (ii>max)
            {
                max = ii;
                maxPath = 'I';
            }

            I[i][j] = max;
            pathI[i][j] = maxPath;

            //W
            double mw = M[i-1][j] + tmw;
            double ww = W[i-1][j] + tww;

            max = mw;
            maxPath = 'M';

            if (ww>max)
            {
                max = ww;
                maxPath = 'W';
            }

            W[i][j] = max;
            pathW[i][j] = maxPath;

            //Z
            double mz = M[i][j-1] + tmz;
            double wz = W[i][j-1] + twz;
            double zz = Z[i][j-1] + tzz;

            max = mz;
            maxPath = 'M';

            if (wz>max)
            {
                max = wz;
                maxPath = 'W';
            }
            if (zz>max)
            {
                max = zz;
                maxPath = 'Z';
            }

            Z[i][j] = max;
            pathZ[i][j] = maxPath;
        }

        M[xlen][ylen] += logTau-logEta;
    }

    if (debug)
    {
        std::cerr << "\n=X=\n";
        printVector(X, xlen+1, ylen+1);
        std::cerr << "\n=Y=\n";
        printVector(Y, xlen+1, ylen+1);
        std::cerr << "\n=M=\n";
        printVector(M, xlen+1, ylen+1);
        std::cerr << "\n=D=\n";
        printVector(D, xlen+1, ylen+1);
        std::cerr << "\n=I=\n";
        printVector(I, xlen+1, ylen+1);
        std::cerr << "\n=W=\n";
        printVector(W, xlen+1, ylen+1);
        std::cerr << "\n=Z=\n";
        printVector(Z, xlen+1, ylen+1);
        std::cerr << "\n=Path X=\n";
        printVector(pathX, xlen+1, ylen+1);
        std::cerr << "\n=Path Y=\n";
        printVector(pathY, xlen+1, ylen+1);
        std::cerr << "\n=Path M=\n";
        printVector(pathM, xlen+1, ylen+1);
        std::cerr << "\n=Path D=\n";
        printVector(pathD, xlen+1, ylen+1);
        std::cerr << "\n=Path I=\n";
        printVector(pathI, xlen+1, ylen+1);
        std::cerr << "\n=Path W=\n";
        printVector(pathW, xlen+1, ylen+1);
        std::cerr << "\n=Path Z=\n";
        printVector(pathZ, xlen+1, ylen+1);
    }

    tracePath();
};

/**
*Compute log likelihood of matched portion of alignment
*/
void LHMM1::computeLogLikelihood(double& llk, std::string& _path, const char* qual)
{

//  std::cerr <<" basequals " << qual << "\n";
//  std::cerr <<" read      " << y << "\n";
//  std::cerr <<" path      " << path << "\n";
    matchedBases = 0;
    mismatchedBases = 0;

    llk = 0;
    double q = 0;
    char lastState = 'S';
    uint32_t xIndex=0, yIndex=0;
    for (uint32_t i=0; i<_path.size(); ++i)
    {
        char state = _path.at(i);
        //std::cerr << state << " " << xIndex << " " << yIndex << "\n";

        if (state=='M')
        {
            q = pl2prob((uint32_t) qual[yIndex]-33);
            llk = lt.log10prod(llk, logEmission(x[xIndex], y[yIndex], q));
            //compute for perfect fit

            if (x[xIndex]== y[yIndex])
            {
                ++matchedBases;
            }
            else
            {
                ++mismatchedBases;
            }

            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMM);
            }
            else if (lastState=='D')
            {
                llk = lt.log10prod(llk, tDM);
            }
            else if (lastState=='I')
            {
                llk = lt.log10prod(llk, tIM);
            }

            //std::cerr << "M) "<< lastState << " "  << llk << " " << x[xIndex] << " " <<  y[yIndex] << " " << pl2prob((uint32_t) qual[yIndex]-33) << " "<< qual[yIndex-1]-33 << "\n";

            ++xIndex;
            ++yIndex;
        }
        else if (state=='D')
        {
            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMD);
            }
            else if (lastState=='D')
            {
                llk = lt.log10prod(llk, tDD);
            }

            //std::cerr << "D) " << lastState << " " << llk << " " << x[xIndex] << " " <<  y[yIndex] << "\n";

            ++xIndex;
        }
        else if (state=='I')
        {
            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMI);
            }
            else if (lastState=='I')
            {
                llk = lt.log10prod(llk, tII);
            }

            //std::cerr << "I) " << lastState << " " << llk << " " << x[xIndex] << " " <<  y[yIndex] << "\n";

            ++yIndex;
        }
        else if (state=='X')
        {
            ++xIndex;
        }
        else if (state=='Y')
        {
            ++yIndex;
        }
        else if (state=='W')
        {
            ++xIndex;
        }
        else if (state=='Z')
        {
            ++yIndex;
        }

        lastState = state;
    }
};

/**
 * Compute log likelihood of matched portion of alignment
 */
void LHMM1::computeLogLikelihood(double& llk, double& perfectllk, std::string& _path, const char* qual)
{

//  std::cerr <<" basequals " << qual << "\n";
//  std::cerr <<" read      " << y << "\n";
//  std::cerr <<" path      " << path << "\n";
    matchedBases = 0;
    mismatchedBases = 0;

    llk = 0;
    perfectllk = 0;
    double q = 0;
    char lastState = 'S';
    uint32_t xIndex=0, yIndex=0;
    for (uint32_t i=0; i<_path.size(); ++i)
    {

        char state = _path.at(i);
        //std::cerr << state << " " << xIndex << " " << yIndex << "\n";

        if (state=='M')
        {
            q = pl2prob((uint32_t) qual[yIndex]-33);
            llk = lt.log10prod(llk, logEmission(x[xIndex], y[yIndex], q));
            //compute for perfect fit
            perfectllk = lt.log10prod(perfectllk, logEmission(x[xIndex], x[xIndex], q));

            if (x[xIndex]== y[yIndex])
            {
                ++matchedBases;
            }
            else
            {
                ++mismatchedBases;
            }

            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMM);
            }
            else if (lastState=='D')
            {
                llk = lt.log10prod(llk, tDM);
            }
            else if (lastState=='I')
            {
                llk = lt.log10prod(llk, tIM);
            }

            perfectllk = lt.log10prod(llk, tMM);

            //std::cerr << "M) "<< lastState << " "  << llk << " " << x[xIndex] << " " <<  y[yIndex] << " " << pl2prob((uint32_t) qual[yIndex]-33) << " "<< qual[yIndex-1]-33 << "\n";

            ++xIndex;
            ++yIndex;
        }
        else if (state=='D')
        {
            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMD);
            }
            else if (lastState=='D')
            {
                llk = lt.log10prod(llk, tDD);
            }

           // std::cerr << "D) " << lastState << " " << llk << " " << x[xIndex] << " " <<  y[yIndex] << "\n";

            ++xIndex;
        }
        else if (state=='I')
        {
            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMI);
            }
            else if (lastState=='I')
            {
                llk = lt.log10prod(llk, tII);
            }

            // std::cerr << "I) " << lastState << " " << llk << " " << x[xIndex] << " " <<  y[yIndex] << "\n";

            ++yIndex;
        }
        else if (state=='X')
        {
            ++xIndex;
        }
        else if (state=='Y')
        {
            ++yIndex;
        }
        else if (state=='W')
        {
            ++xIndex;
        }
        else if (state=='Z')
        {
            ++yIndex;
        }

        lastState = state;
    }
};

/**
*Compute log likelihood of matched portion of alignment
*/
void LHMM1::computeLogLikelihood(double& llk, const char* qual)
{
    tracePath();
    llk = 0;
    char lastState = 'S';
    uint32_t xIndex=0, yIndex=0;

    for (uint32_t i=0; i<path.size(); ++i)
    {
        //std::cerr << xIndex << " " << yIndex << "\n";

        char state = path.at(i);

        if (state=='M')
        {
            llk = lt.log10prod(llk, logEmission(x[xIndex], y[yIndex], pl2prob((uint32_t) qual[yIndex-1]-33)));
            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMM);
            }
            else if (lastState=='D')
            {
                llk = lt.log10prod(llk, tDM);
            }
            else if (lastState=='I')
            {
                llk = lt.log10prod(llk, tIM);
            }

            ++xIndex;
            ++yIndex;
        }
        else if (state=='D')
        {
            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMD);
            }
            else if (lastState=='D')
            {
                llk = lt.log10prod(llk, tDD);
            }
            ++xIndex;
        }
        else if (state=='I')
        {
            if (lastState=='M')
            {
                llk = lt.log10prod(llk, tMI);
            }
            else if (lastState=='I')
            {
                llk = lt.log10prod(llk, tII);
            }
            ++yIndex;
        }
        else if (state=='X')
        {
            ++xIndex;
        }
        else if (state=='Y')
        {
            ++yIndex;
        }
        else if (state=='W')
        {
            ++xIndex;
        }
        else if (state=='Z')
        {
            ++yIndex;
        }

        lastState = state;
    }
}

std::string& LHMM1::getPath()
{
    return path;
}

double LHMM1::logEmissionOdds(char readBase, char probeBase, double e)
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

double LHMM1::logEmission(char readBase, char probeBase, double e)
{
    //4 encodes for N
    if (readBase=='N' || probeBase=='N')
    {
        //silent match
        return 0;
    }

    if (readBase!=probeBase)
    {
        return log10(e/3);
    }
    else
    {
        return log10(1-e);
    }
};

double LHMM1::emission(char readBase, char probeBase, double e)
{
    //4 encodes for N
    if (readBase=='N' || probeBase=='N')
    {
        //silent match
        return 0;
    }

    if (readBase!=probeBase)
    {
        return e/3;
    }
    else
    {
        return 1-e;
    }
};

std::string LHMM1::reverse(std::string s)
{
    std::string rs;
    for (std::string::reverse_iterator rit=s.rbegin() ; rit < s.rend(); rit++ )
    {
        rs.push_back(*rit);
    }

    return rs;
};

void LHMM1::printVector(std::vector<std::vector<double> >& v, uint32_t xLen, uint32_t yLen)
{
    for (uint32_t i=0; i<xLen; ++i)
    {
        for (uint32_t j=0; j<yLen; ++j)
        {
          std::cerr << (v[i][j]==-DBL_MAX?-1000:v[i][j]) << "\t";
        }

        std::cerr << "\n";
    }
};

void LHMM1::printVector(std::vector<std::vector<char> >& v, uint32_t xLen, uint32_t yLen)
{
    for (uint32_t i=0; i<xLen; ++i)
    {
        for (uint32_t j=0; j<yLen; ++j)
        {
          std::cerr << v[i][j] << "\t";
        }

        std::cerr << "\n";
    }
};

void LHMM1::printVector(std::vector<std::vector<double> >& v)
{
    for (uint32_t i=0; i<v.size(); ++i)
    {
        for (uint32_t j=0; j<v[i].size(); ++j)
        {
          std::cerr << v[i][j] << "\t";
        }

        std::cerr << "\n";
    }
};

void LHMM1::printAlignment(std::string& pad, std::stringstream& log)
{
    std::stringstream xAligned;
    std::stringstream yAligned;
    std::stringstream qualAligned;
    uint32_t xIndex=0, yIndex=0;
    for (uint32_t i=0; i<path.size(); ++i)
    {
        char state = path.at(i);
        if (state=='S' || state=='E')
        {
            xAligned << state;
            yAligned << state;
            qualAligned << state;
        }
        else if (state=='X'||state=='W')
        {
            xAligned << x[xIndex++];
            yAligned << '.';
            qualAligned << '.';
        }
        else if (state=='M')
        {
            xAligned << x[xIndex++];
            yAligned << y[yIndex];
            qualAligned << qual[yIndex++];
        }
        else if (state=='D')
        {
            xAligned << x[xIndex++];
            yAligned << '-';
            qualAligned << '-';
        }
        else if (state=='I')
        {
            xAligned << '-';
            yAligned << y[yIndex];
            qualAligned << qual[yIndex++];
        }
        else if (state=='Y'||state=='Z')
        {
            xAligned << '.';
            yAligned << y[yIndex];
            qualAligned << qual[yIndex++];
        }
    }

    log << pad << "X    " << xAligned.str() << "E\n";
    log << pad << "Path " << path << "E\n";
    log << pad << "Y    " << yAligned.str() << "E\n";
    log << pad << "Qual " << qualAligned.str() << "E\n\n";

    log << pad << "Alignment in Probe : [" << matchStartX << "," << matchEndX<< "]\n";
    log << pad << "Alignment in Read  : [" << matchStartY << "," << matchEndY<< "]\n";
    log << pad << "# Aligned Bases : " << noBasesAligned << "\n";
    log << pad << "# Matched Bases : " << matchedBases << "\n";
    log << pad << "# Mismatched Bases : " << mismatchedBases << "\n";

    log << pad << "# INS in Probe : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='I')
        {
            log << "[" << indelStartsInX[i] << "," << indelEndsInX[i] <<"] ";
        }
    }
    log << "\n";

    log << pad << "# DELS in Probe : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='D')
        {
            log << "[" << indelStartsInX[i] << "," << indelEndsInX[i] <<"] ";
        }
    }
    log << "\n";

    log << pad << "# INS in Read : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='I')
        {
            log << "[" << indelStartsInY[i] << "," << indelEndsInY[i] <<"] ";
        }
    }
    log << "\n";

    log << pad << "# DELS in Read : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='D')
        {
            log << "[" << indelStartsInY[i] << "," << indelEndsInY[i] <<"] ";
        }
    }
    log << "\n";

    log << pad << "Max Log odds    : " << maxLogOdds << "\n";
};

void LHMM1::printAlignment(std::string& pad)
{
    //std::cerr << "\tAlignments\n" ;

    //trace back
    //std::cerr << "\tBefore alignment" << "\n";
    //std::cerr << "\tX    " << x << "\n";
    //std::cerr << "\tY    " << y << "\n";
    //std::cerr << "\n";
    //std::cerr << "\tAfter alignment" << "\n";
    //use path to align x and y
    std::stringstream xAligned;
    std::stringstream yAligned;
    std::stringstream qualAligned;
    uint32_t xIndex=0, yIndex=0;
    for (uint32_t i=0; i<path.size(); ++i)
    {
        char state = path.at(i);
        if (state=='S' || state=='E')
        {
            xAligned << state;
            yAligned << state;
            qualAligned << state;
        }
        else if (state=='X'||state=='W')
        {
            xAligned << x[xIndex++];
            yAligned << '.';
            qualAligned << '.';
        }
        else if (state=='M')
        {
            xAligned << x[xIndex++];
            yAligned << y[yIndex];
            qualAligned << qual[yIndex++];
        }
        else if (state=='D')
        {
            xAligned << x[xIndex++];
            yAligned << '-';
            qualAligned << '-';
        }
        else if (state=='I')
        {
            xAligned << '-';
            yAligned << y[yIndex];
            qualAligned << qual[yIndex++];
        }
        else if (state=='Y'||state=='Z')
        {
            xAligned << '.';
            yAligned << y[yIndex];
            qualAligned << qual[yIndex++];
        }
    }

    std::cerr << pad << "X    " << xAligned.str() << "E\n";
    std::cerr << pad << "Path " << path << "E\n";
    std::cerr << pad << "Y    " << yAligned.str() << "E\n";
    std::cerr << pad << "Qual " << qualAligned.str() << "E\n\n";

    std::cerr << pad << "Alignment in Probe : [" << matchStartX << "," << matchEndX<< "]\n";
    std::cerr << pad << "Alignment in Read  : [" << matchStartY << "," << matchEndY<< "]\n";
    std::cerr << pad << "# Aligned Bases : " << noBasesAligned << "\n";
    std::cerr << pad << "# Matched Bases : " << matchedBases << "\n";
    std::cerr << pad << "# Mismatched Bases : " << mismatchedBases << "\n";

    std::cerr << pad << "# INS in Probe : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='I')
        {
            std::cerr << "[" << indelStartsInX[i] << "," << indelEndsInX[i] <<"] ";
        }
    }
    std::cerr << "\n";

    std::cerr << pad << "# DELS in Probe : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='D')
        {
            std::cerr << "[" << indelStartsInX[i] << "," << indelEndsInX[i] <<"] ";
        }
    }
    std::cerr << "\n";

    std::cerr << pad << "# INS in Read : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='I')
        {
            std::cerr << "[" << indelStartsInY[i] << "," << indelEndsInY[i] <<"] ";
        }
    }
    std::cerr << "\n";

    std::cerr << pad << "# DELS in Read : ";
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i)
    {
        if (indelStatusInPath[i]=='D')
        {
            std::cerr << "[" << indelStartsInY[i] << "," << indelEndsInY[i] <<"] ";
        }
    }
    std::cerr << "\n";

    std::cerr << pad << "Max Log odds    : " << maxLogOdds << "\n";
};

void LHMM1::printAlignment()
{
    std::string pad = "\t";
    printAlignment(pad);
};

double LHMM1::score(char a, char b)
{
    if (a==b)
        return 0;
    else
        return -1;
};

/**
 * Checks if deletion exists in alignment.
 */
bool LHMM1::deletion_start_exists(uint32_t pos, uint32_t& rpos)
{
    rpos = 0;
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i) 
    {
        if (indelStatusInPath[i]=='D' &&
            indelStartsInX[i]==pos)
        {
            rpos = indelStartsInY[i];
            return true;
        }    
    }
    
    return false;
}

/**
 * Checks if insertion exists in alignment.
 */
bool LHMM1::insertion_start_exists(uint32_t pos, uint32_t& rpos)
{
    rpos = 0;
    for (uint32_t i=0; i<indelStatusInPath.size(); ++i) 
    {
        if (indelStatusInPath[i]=='I' &&
            indelStartsInX[i]==pos)
        {
            rpos = indelStartsInY[i];
            return true;
        }    
    }
    
    return false;
}