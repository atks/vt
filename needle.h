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
#ifndef NEEDLE_H
#define NEEDLE_H

#include "utils.h"

class NWParameters
{
    public:

    int score_match;
    int score_mismatch;
    int score_gap;

    NWParameters() : score_match(0), score_mismatch(-1), score_gap(-1)
    {}
};

class NeedlemanWunsch
{
    public:

    // in CIGAR: "?", "D", "M", "X", "I"
    enum Traceback { EMPTY, CIGAR_D, CIGAR_M, CIGAR_X, CIGAR_I };

    bool debug;
    const char* ref;
    const char* read;

    int len_ref;
    int len_read;

    std::vector<int> scores;
    std::vector<Traceback> matrix;
    std::vector<Traceback> trace;

    NWParameters params;

    /**
     * Constructor.
     */
    NeedlemanWunsch(bool debug=false);

    void set_read(const char * read) 
    {
        this->read = read;
    }

    /**
     * Align and compute genotype likelihood.
     */
    void align(const char* ref, const char* read);

    /**
     * Trace path after alignment and write to trace.
     */
    void trace_path();

    /**
     * Prints an alignment.
     */
    void print_alignment()
    {
        print_alignment("");
    }

    /**
     * Prints an alignment with padding.
     */
    void print_alignment(std::string const & pad);
};

#endif  // ifndef NEEDLE_H
