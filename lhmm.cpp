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

#include "lhmm.h"

LHMM::LHMM()
{
    delta = 0.001;
    epsilon = 0.5;
    tau = 0.1;
    eta = 0.001;

    logOneSixteenth = log10(1.0/16.0);

    transition[X][X] = 0; //log10((1-eta)/(1-eta));
    transition[X][Y] = 0; //log10((eta*(1-eta))/(eta*(1-eta)));
    transition[Y][Y] = 0; //log10((1-eta)/(1-eta));

    transition[X][M] = log10(1/((1-eta)*(1-eta))); //log10((eta*eta)/(eta*eta*(1-eta)*(1-eta)));
    transition[Y][M] = transition[X][M]; //log10((eta)/(eta*(1-eta)*(1-eta)));
    transition[S][M] = log10(1/(eta*(1-eta)*(1-eta))); //log10((eta*eta)/(eta*eta*eta*(1-eta)*(1-eta)));
    transition[M][M] = log10((1-2*delta-tau)/((1-eta)*(1-eta)));
    transition[I][M] = log10((1-epsilon)/(1-eta));
    transition[D][M] = transition[I][M]; //log10((1-epsilon)/(1-eta));

    transition[M][D] = log10(delta/(1-eta));
    transition[D][D] = log10(epsilon/(1-eta));

    transition[M][I] = transition[M][D]; //log10(delta/(1-eta));
    transition[I][I] = transition[D][D]; //log10(epsilon/(1-eta));

    transition[M][W] = log10(tau/eta); //log10((tau*(1-eta))/(eta*(1-eta)));
    transition[W][W] = 0; //log10((1-eta)/(1-eta));

    transition[M][Z] = transition[M][W]; //log10((tau*eta*(1-eta))/(eta*eta*(1-eta)));
    transition[W][Z] = 0; //log10((eta*(1-eta))/(eta*(1-eta)));
    transition[Z][Z] = 0; //log10((1-eta)/(1-eta));

    scoreX = new double[MAXLEN*MAXLEN];
    scoreY = new double[MAXLEN*MAXLEN];
    scoreM = new double[MAXLEN*MAXLEN];
    scoreI = new double[MAXLEN*MAXLEN];
    scoreD = new double[MAXLEN*MAXLEN];
    scoreW = new double[MAXLEN*MAXLEN];
    scoreZ = new double[MAXLEN*MAXLEN];
    
    pathX = new char[MAXLEN*MAXLEN];
    pathY = new char[MAXLEN*MAXLEN];
    pathM = new char[MAXLEN*MAXLEN];
    pathI = new char[MAXLEN*MAXLEN];
    pathD = new char[MAXLEN*MAXLEN];
    pathW = new char[MAXLEN*MAXLEN];
    pathZ = new char[MAXLEN*MAXLEN];
    
    //assume alignments can't possibly be maxLength bases or more
    for (int32_t i=0; i<MAXLEN; ++i)
    {
        for (int32_t j=0; j<MAXLEN; ++j)
        {
            scoreX[i*MAXLEN+j] = -DBL_MAX;
            scoreY[i*MAXLEN+j] = -DBL_MAX;
            scoreM[i*MAXLEN+j] = -DBL_MAX;
            scoreI[i*MAXLEN+j] = -DBL_MAX;
            scoreD[i*MAXLEN+j] = -DBL_MAX;
            scoreW[i*MAXLEN+j] = -DBL_MAX;
            scoreZ[i*MAXLEN+j] = -DBL_MAX;
            
            if (j)
            {
                pathX[i*MAXLEN+j] = 'Y';
                pathY[i*MAXLEN+j] = 'Y';
                pathM[i*MAXLEN+j] = 'Y';
                pathI[i*MAXLEN+j] = 'Y';
                pathD[i*MAXLEN+j] = 'Y';
                pathW[i*MAXLEN+j] = 'Y';
                pathZ[i*MAXLEN+j] = 'Y';
            }
            else
            {
                pathX[i*MAXLEN+j] = 'X';
                pathY[i*MAXLEN+j] = 'X';
                pathM[i*MAXLEN+j] = 'X';
                pathI[i*MAXLEN+j] = 'X';
                pathD[i*MAXLEN+j] = 'X';
                pathW[i*MAXLEN+j] = 'X';
                pathZ[i*MAXLEN+j] = 'X';
            }
        }
    }

    logEta = log10(eta);
    logTau = log10(tau);

    scoreX[0*MAXLEN+0] = 0;
    scoreY[0*MAXLEN+0] = 0;
    scoreM[0*MAXLEN+0] = 0;
    scoreW[0*MAXLEN+0] = 0;
    scoreZ[0*MAXLEN+0] = 0;
    pathX[0*MAXLEN+0] = 'N';
    pathX[1*MAXLEN+0] = 'S';
    pathY[0*MAXLEN+0] = 'N';
    pathY[0*MAXLEN+1] = 'S';
    pathM[0*MAXLEN+0] = 'N';
    pathM[1*MAXLEN+1] = 'S';

    for (uint32_t k=1; k<MAXLEN; ++k)
    {
        scoreX[k*MAXLEN+0] = scoreX[(k-1)*MAXLEN+0] + transition[X][X];
        scoreX[0*MAXLEN+k] = -DBL_MAX;
        scoreY[k*MAXLEN+0] = -DBL_MAX;
        scoreY[0*MAXLEN+k] = scoreY[0*MAXLEN+(k-1)] + transition[Y][Y];
        scoreW[k*MAXLEN+0] = scoreW[(k-1)*MAXLEN+0] + transition[W][W];
        scoreW[0*MAXLEN+k] = -DBL_MAX;
        scoreZ[k*MAXLEN+0] = -DBL_MAX;
        scoreZ[0*MAXLEN+k] = scoreZ[0*MAXLEN+(k-1)] + transition[Z][Z];
    }

    scoreX[0*MAXLEN+0] = -DBL_MAX;
    scoreY[0*MAXLEN+0] = -DBL_MAX;
    scoreW[0*MAXLEN+0] = -DBL_MAX;
    scoreZ[0*MAXLEN+0] = -DBL_MAX;

};

/**
Align and compute genotype likelihood.
*/
void LHMM::align(double& llk, const char* x, const char* y, const char* qual, bool debug)
{
    //std::cerr << "Running this\n" ;
    this->x = x;
    this->y = y;
    this->qual = qual;

    //adds a starting character at the fron of each string that must be matched
    xlen = strlen(x);
    ylen = strlen(y);
    double max = 0;
    char maxPath = 'X';

    //construct possible solutions
    for (uint32_t i=1; i<=xlen; ++i)
    {
        for (uint32_t j=1; j<=ylen; ++j)
        {
            //X
            double xx = scoreX[(i-1)*MAXLEN+j] + transition[X][X];

            max = xx;
            maxPath = 'X';

            scoreX[i*MAXLEN+j] = max;
            pathX[i*MAXLEN+j] = maxPath;

            //Y
            double xy = scoreX[i*MAXLEN+(j-1)] + transition[X][Y];
            double yy = scoreY[i*MAXLEN+(j-1)] + transition[Y][Y];

            max = xy;
            maxPath = 'X';

            if (yy>max)
            {
                max = yy;
                maxPath = 'Y';
            }

            scoreY[i*MAXLEN+j] = max;
            pathY[i*MAXLEN+j] = maxPath;

            //M
            double xm = scoreX[(i-1)*MAXLEN+(j-1)] + transition[X][M];
            double ym = scoreY[(i-1)*MAXLEN+(j-1)] + transition[Y][M];
            double mm = scoreM[(i-1)*MAXLEN+(j-1)] + ((i==1&&j==1) ? transition[S][M] : transition[M][M]);
            double im = scoreI[(i-1)*MAXLEN+(j-1)] + transition[I][M];
            double dm = scoreD[(i-1)*MAXLEN+(j-1)] + transition[D][M];

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

            scoreM[i*MAXLEN+j] = max + logEmissionOdds(x[i-1], y[j-1], lt.pl2prob((uint32_t) qual[j-1]-33));
            pathM[i*MAXLEN+j] = maxPath;

            //D
            double md = scoreM[(i-1)*MAXLEN+j] + transition[M][D];
            double dd = scoreD[(i-1)*MAXLEN+j] + transition[D][D];

            max = md;
            maxPath = 'M';

            if (dd>max)
            {
                max = dd;
                maxPath = 'D';
            }

            scoreD[i*MAXLEN+j] = max;
            pathD[i*MAXLEN+j] = maxPath;

            //I
            double mi = scoreM[i*MAXLEN+(j-1)] + transition[M][I];
            double ii = scoreI[i*MAXLEN+(j-1)] + transition[I][I];

            max = mi;
            maxPath = 'M';

            if (ii>max)
            {
                max = ii;
                maxPath = 'I';
            }

            scoreI[i*MAXLEN+j] = max;
            pathI[i*MAXLEN+j] = maxPath;

            //W
            double mw = scoreM[(i-1)*MAXLEN+j] + transition[M][W];
            double ww = scoreW[(i-1)*MAXLEN+j] + transition[W][W];

            max = mw;
            maxPath = 'M';

            if (ww>max)
            {
                max = ww;
                maxPath = 'W';
            }

            scoreW[i*MAXLEN+j] = max;
            pathW[i*MAXLEN+j] = maxPath;

            //Z
            double mz = scoreM[i*MAXLEN+(j-1)] + transition[M][Z];
            double wz = scoreW[i*MAXLEN+(j-1)] + transition[W][Z];
            double zz = scoreZ[i*MAXLEN+(j-1)] + transition[Z][Z];

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

            scoreZ[i*MAXLEN+j] = max;
            pathZ[i*MAXLEN+j] = maxPath;
        }

        scoreM[xlen*MAXLEN+ylen] += logTau-logEta;
    }

    if (debug)
    {
//        std::cerr << "\n=X=\n";
//        printVector(scoreX, xlen+1, ylen+1);
//        std::cerr << "\n=Y=\n";
//        printVector(scoreY, xlen+1, ylen+1);
//        std::cerr << "\n=M=\n";
//        printVector(scoreM, xlen+1, ylen+1);
//        std::cerr << "\n=D=\n";
//        printVector(scoreD, xlen+1, ylen+1);
//        std::cerr << "\n=I=\n";
//        printVector(scoreI, xlen+1, ylen+1);
//        std::cerr << "\n=W=\n";
//        printVector(scoreW, xlen+1, ylen+1);
//        std::cerr << "\n=Z=\n";
//        printVector(scoreZ, xlen+1, ylen+1);
//        std::cerr << "\n=Path X=\n";
//        printVector(pathX, xlen+1, ylen+1);
//        std::cerr << "\n=Path Y=\n";
//        printVector(pathY, xlen+1, ylen+1);
//        std::cerr << "\n=Path M=\n";
//        printVector(pathM, xlen+1, ylen+1);
//        std::cerr << "\n=Path D=\n";
//        printVector(pathD, xlen+1, ylen+1);
//        std::cerr << "\n=Path I=\n";
//        printVector(pathI, xlen+1, ylen+1);
//        std::cerr << "\n=Path W=\n";
//        printVector(pathW, xlen+1, ylen+1);
//        std::cerr << "\n=Path Z=\n";
//        printVector(pathZ, xlen+1, ylen+1);
    }

    tracePath();
};

double LHMM::logEmissionOdds(char readBase, char probeBase, double e)
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

//bool LHMM::containsIndel()
//{
//    return (indelStatusInPath.size()!=0);
//}



/**
 * Updates matchStart, matchEnd, globalMaxPath and path.
 */
void LHMM::tracePath()
{
    double globalMax = scoreM[xlen*MAXLEN+ylen];
    char globalMaxPath = 'M';
    if (scoreW[xlen*MAXLEN+ylen]>globalMax)
    {
        globalMax = scoreW[xlen*MAXLEN+ylen];
        globalMaxPath = 'W';
    }
    if (scoreZ[xlen*MAXLEN+ylen]>globalMax)
    {
        globalMax = scoreZ[xlen*MAXLEN+ylen];
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

void LHMM::tracePath(std::stringstream& ss, char state, uint32_t i, uint32_t j)
{
    //std::cout << state << "\n";
    if (i>0 && j>0)
    {
        if (state=='X')
        {
            //std::cerr << pathX[i*MAXLEN+j] << " " << i << " " << j << " " << matchStartY << "\n";
            ss << pathX[i*MAXLEN+j];
            tracePath(ss, pathX[i*MAXLEN+j], i-1, j);
        }
        else if (state=='Y')
        {
            //std::cerr << pathY[i*MAXLEN+j] << " " << i << " " << j << " " << matchStartY << "\n";
            ss << pathY[i*MAXLEN+j];
            tracePath(ss, pathY[i*MAXLEN+j], i, j-1);
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
            if (matchStartX==-1 && (pathM[i*MAXLEN+j] =='X' || pathM[i*MAXLEN+j]=='Y'))
            {
               matchStartX = i;
               matchStartY = j;

               //std::cerr << "set match starts "  << matchStartX  << " " << matchStartY  <<  "\n";
            }


            //std::cout << pathM[i*MAXLEN+j] << " " << i << " " << j << " " << matchEndY << "\n";
            ss << pathM[i*MAXLEN+j];
            tracePath(ss, pathM[i*MAXLEN+j], i-1, j-1);
            ++noBasesAligned;
        }
        else if (state=='I')
        {
            //std::cout << pathI[i*MAXLEN+j] << " " << i << " " << j << "\n";
            ss << pathI[i*MAXLEN+j];
            tracePath(ss, pathI[i*MAXLEN+j], i, j-1);
        }
        else if (state=='D')
        {
            //std::cout << pathD[i*MAXLEN+j] << " " << i << " " << j << "\n";
            ss << pathD[i*MAXLEN+j];
            tracePath(ss, pathD[i*MAXLEN+j], i-1, j);
            ++noBasesAligned;
        }
        else if (state=='W')
        {
            if (matchEndX==-1 && pathW[i*MAXLEN+j] =='M')
            {
                matchEndX = i-1;
                matchEndY = j;
            }


            //std::cout << pathW[i*MAXLEN+j] << " " << i << " " << j << "\n";
            ss << pathW[i*MAXLEN+j];
            tracePath(ss, pathW[i*MAXLEN+j], i-1, j);
        }
        else if (state=='Z')
        {
            if (matchEndX==-1 && pathZ[i*MAXLEN+j] =='M')
            {
                matchEndX = i;
                matchEndY = j-1;
            }

            //std::cout << pathZ[i*MAXLEN+j] << " " << i << " " << j << "\n";
            ss << pathZ[i*MAXLEN+j];
            tracePath(ss, pathZ[i*MAXLEN+j], i, j-1);
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

        // std::cout << pathY[i*MAXLEN+j] << " " << i << " " << j << "\n";
        ss << pathY[i*MAXLEN+j];
        tracePath(ss, pathY[i*MAXLEN+j], i, j-1);
    }
    else if (i>0 && j==0)
    {
        if (matchStartY==-1)
        {
           matchStartY = j+1;
        }
        // std::cout << pathX[i*MAXLEN+j] << " " << i << " " << j << "\n";
        ss << pathX[i*MAXLEN+j];
        tracePath(ss, pathX[i*MAXLEN+j], i-1, j);
    }
    else
    {
        //std::cout << "\n";
        //ss << pathX[i*MAXLEN+j];
    }
}

/**
 *Left align indels in an alignment
 */
void LHMM::left_align()
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
//                std::cerr << "scoreX[" << indelStartsInX[i]-1 << "] "  << x[indelStartsInX[i]-1] << "\n";
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
//                std::cerr << "scoreX[" << indelStartsInX[i]-1 << "] "  << x[indelStartsInX[i]-1] << "\n";
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

std::string& LHMM::getPath()
{
    return path;
}



double LHMM::logEmission(char readBase, char probeBase, double e)
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

double LHMM::emission(char readBase, char probeBase, double e)
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

std::string LHMM::reverse(std::string s)
{
    std::string rs;
    for (std::string::reverse_iterator rit=s.rbegin() ; rit < s.rend(); rit++ )
    {
        rs.push_back(*rit);
    }

    return rs;
};

void LHMM::printVector(double v[][MAXLEN], uint32_t xLen, uint32_t yLen)
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

void LHMM::printVector(char v[][MAXLEN], uint32_t xLen, uint32_t yLen)
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

void LHMM::printVector(double v[][MAXLEN])
{
    for (uint32_t i=0; i<MAXLEN; ++i)
    {
        for (uint32_t j=0; j<MAXLEN; ++j)
        {
          std::cerr << v[i][j] << "\t";
        }

        std::cerr << "\n";
    }
};

void LHMM::printAlignment(std::string& pad, std::stringstream& log)
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

void LHMM::printAlignment(std::string& pad)
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

void LHMM::printAlignment()
{
    std::string pad = "\t";
    printAlignment(pad);
};

double LHMM::score(char a, char b)
{
    if (a==b)
        return 0;
    else
        return -1;
};

/**
 * Checks if deletion exists in alignment.
 */
bool LHMM::deletion_start_exists(uint32_t pos, uint32_t& rpos)
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
bool LHMM::insertion_start_exists(uint32_t pos, uint32_t& rpos)
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