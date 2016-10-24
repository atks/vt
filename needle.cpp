/* The MIT License

   Copyright (c) 2015 Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>

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

#include "needle.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

NeedlemanWunsch::NeedlemanWunsch(bool debug)
{
    this->debug = debug;
}

void NeedlemanWunsch::align(const char* ref, const char* read)
{
    this->ref = ref;
    this->read = read;
    this->len_ref = strlen(ref);
    this->len_read = strlen(read);

    matrix.clear();
    matrix.resize((len_ref + 1) * (len_read + 1), EMPTY);

    // fill first column and row
    scores.resize(len_read + 1);
    for (unsigned i = 0; i <= len_ref; ++i)
        matrix.at(i * (len_read + 1)) = CIGAR_D;
    for (unsigned i = 0; i <= len_read; ++i)
    {
        matrix.at(i) = CIGAR_I;
        scores.at(i) = i * params.score_gap;
    }
    matrix.at(0) = CIGAR_M;

    // fill remainder of the matrix
    int offset = 0;
    int score_n = 0;
    for (unsigned i = 1; i <= len_ref; ++i)
    {
        offset += len_read + 1;
        score_n += params.score_gap;
        int best_score = 0;
        for (unsigned j = 1; j <= len_read; ++j)
        {
            int score_nw = scores.at(j - 1);
            int score_w = scores.at(j);

            int score_diag = score_nw;
            Traceback best_dir;
            if (ref[i - 1] == read[j - 1])
            {
                best_dir = CIGAR_M;
                score_diag += params.score_match;
            }
            else
            {
                best_dir = CIGAR_X;
                score_diag += params.score_mismatch;
            }

            best_score = score_diag;

            if (score_w + params.score_gap > best_score)
            {
                best_score = score_w + params.score_gap;
                best_dir = CIGAR_D;
            }
            else if (score_n + params.score_gap > best_score)
            {
                best_score = score_n + params.score_gap;
                best_dir = CIGAR_I;
            }

            scores.at(j - 1) = score_n;
            score_n = best_score;
            matrix.at(offset + j) = best_dir;
        }
        scores.back() = best_score;
    }

}

void NeedlemanWunsch::trace_path()
{
    int i = len_ref;
    int j = len_read;
    int k = (len_read + 1) * (len_ref + 1) - 1;

    trace.clear();

    while (i>0 || j>0)
    {
        trace.push_back(matrix.at(k));
        switch ((int32_t) matrix.at(k))
        {
            case CIGAR_X:
            case CIGAR_M:
                --i;
                --j;
                k -= len_read + 2;
                break;
            case CIGAR_I:
                --j;
                --k;
                break;
            case CIGAR_D:
                --i;
                k -= len_read + 1;
                break;
        }
    }

    std::reverse(trace.begin(), trace.end());
}

void NeedlemanWunsch::print_alignment(std::string const & pad)
{
    // print reference
    std::cerr << pad;
    for (unsigned i = 0, j = 0; i < trace.size(); ++i)
    {
        if (trace.at(i) == CIGAR_X || trace[i] == CIGAR_M || trace.at(i) == CIGAR_D)
            std::cerr << ref[j++];
        else
            std::cerr << "-";
    }
    std::cerr << "\n";

    // print trace
    std::cerr << pad;
    for (unsigned i = 0, j = 0, k = 0; i < trace.size(); ++i)
    {
        if (trace[i] == CIGAR_X || trace[i] == CIGAR_M)
            std::cerr << ((ref[j++] == read[k++]) ? '|' : ' ');
        else if (trace[i] == CIGAR_D)
        {
            j++;
            std::cerr << " ";
        }
        else  // trace[i] == CIGAR_I
        {
            k++;
            std::cerr << " ";
        }
    }
    std::cerr << "\n";

    // print read
    std::cerr << pad;
    for (unsigned i = 0, j = 0; i < trace.size(); ++i)
    {
        if (trace[i] == CIGAR_X || trace[i] == CIGAR_M || trace[i] == CIGAR_I)
            std::cerr << read[j++];
        else
            std::cerr << "-";
    }
    std::cerr << "\n";
}
