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

#ifndef PEEK_H
#define PEEK_H

#include <cstdio>
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "program.h"
#include "hts_utils.h"
#include "bcf_ordered_reader.h"
#include "variant_manip.h"
#include "filter.h"

#define MONOMORPHIC 0
#define BIALLELIC 1
#define TRIALLELIC 2
#define TETRAALLELIC 3
#define GE_PENTAALLELIC 4
#define GE_TRIALLELIC 5
#define GE_TETRAALLELIC 6
#define POLYMORPHIC 7

#define NO_ALLELE_CATEGORIES 8

void peek(int argc, char ** argv);
    
#endif
