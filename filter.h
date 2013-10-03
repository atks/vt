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

#ifndef FILTER_H
#define FILTER_H

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/join.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/tokenizer.hpp"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "hts_utils.h"
#include "synced_reader.h"
#include "ordered_reader.h"
#include "ordered_writer.h"
#include "variant_manip.h"

#define LT 0
#define LE 1
#define EQ 2
#define GT 3
#define GE 4

/**
 * Class for filtering VCF records.
 */
class Filter
{
    public:
        
    std::string tag;
    int32_t comparison;
    float value;

    Filter() {};
            
    Filter(std::string tag, int32_t comparison, float value);
    
    bool apply(bcf_hdr_t *h, bcf1_t *v);
    
    void parse(std::string filter);
};

#endif