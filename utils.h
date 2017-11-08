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

#include <unistd.h>
#include <climits>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>
#include <list>
#include <queue>
#include "hts_utils.h"

#define fix_neg_zero(f) ((f)==0?0:(f))

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char* delims, std::string& str, uint32_t limit=UINT_MAX, bool clear=true, bool collapse=true);

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char* delims, const char* str, uint32_t limit=UINT_MAX, bool clear=true, bool collapse=true);

/**
 * Joins a vector of strings into a string - PERL style
 */
std::string join(std::vector<std::string>& vec, std::string delim);

/**
 * Joins a map with with string keys into a string - PERL style
 */
std::string join(std::map<std::string, int32_t>& map, std::string delim);

/**
 * Casts a string into int32.  Returns true if successful.
 */
bool str2int32(std::string& s, int32_t& i);
    
/**
 * Casts a string into int32.  Aborts if unsuccessful.
 */
int32_t str2int32(std::string& s);    

/**
 * Casts a string into uint32.  Returns true if successful.
 */
bool str2uint32(std::string& s, uint32_t& i);

/**
 * Casts a string into double.  Returns true if successful.
 */
bool str2double(std::string& s, double& d);

/**
 * Returns a string in lower case.
 */
std::string to_lower(std::string& s);
    
/**
 * Appends current working directoy to a path.
 * Returns true if successful.
 */
bool append_cwd(std::string& path);

#endif