/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#ifndef PREGEX_H
#define PREGEX_H

#include "pcre2.h"
#include <cstdlib>
#include <cstring>
#include <string>

/**
 * Class for PERL regular expressions.
 * Wrapper for pcre2.
 */
class PERLregex
{
    public:
    std::string regex;
    pcre2_code *re;
    PCRE2_SPTR pattern;     /* PCRE2_SPTR is a pointer to unsigned code units of */
    PCRE2_SPTR subject;     /* the appropriate width (8, 16, or 32 bits). */
    PCRE2_SPTR name_table;

    int32_t crlf_is_newline;
    int32_t errornumber;
    int32_t find_all;
    int32_t i;
    int32_t namecount;
    int32_t name_entry_size;
    int32_t rc;
    int32_t utf8;

    uint32_t option_bits;
    uint32_t newline;

    PCRE2_SIZE erroroffset;
    PCRE2_SIZE *ovector;

    size_t subject_length;
    pcre2_match_data *match_data;

    /**
     * Constructor.
     */
    PERLregex();

    /**
     * Destructor.
     */
    ~PERLregex();
    
    /**
     * Sets the regular expression amd compiles it for matching later.
     */
    void set(std::string& regex);

    /**
     * Sets the regular expression amd compiles it for matching later.
     */
    void set(char* regex);
        
    /**
     * Matches a text against a regular expression that has been compiled in set().
     */
    bool match(std::string& text);
        
    /**
     * Matches a text against a regular expression that has been compiled in set().
     */
    bool match(char* text);    
};

#endif
