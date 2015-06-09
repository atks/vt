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

#ifndef UTILS_H
#define UTILS_H

#include <sys/stat.h>
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char* delims, std::string& str, uint32_t limit=UINT_MAX, bool clear=true, bool collapse=true);

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char* delims, const char* str, uint32_t limit=UINT_MAX, bool clear=true, bool collapse=true);


/**
 * Casts a string into int32.  Returns true if successful.
 */
bool str2int32(std::string& s, int32_t& i);

/**
 * Casts a string into uint32.  Returns true if successful.
 */
bool str2uint32(std::string& s, uint32_t& i);

/**
 * Appends cuurent working directoy to a path.
 * Returns true if successful.
 */
bool append_cwd(std::string& path);

#endif
