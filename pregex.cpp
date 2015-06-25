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

#include "pregex.h"

/**
 * Constructor.
 */
PERLregex::PERLregex()
{
   regex = "";
   re = NULL;
   match_data = NULL;
};

/**
 * Destructor.
 */
PERLregex::~PERLregex()
{
    if (re) pcre2_code_free(re);
    if (match_data) pcre2_match_data_free(match_data);
};

/**
 * Sets the regular expression amd compiles it for matching later.
 */
void PERLregex::set(std::string& regex)
{
    set(const_cast<char*>(regex.c_str()));
};

/**
 * Sets the regular expression amd compiles it for matching later.
 */
void PERLregex::set(char* regex)
{
    if (re) pcre2_code_free(re);
    if (match_data) pcre2_match_data_free(match_data);

    this->regex = regex;
    pattern = (PCRE2_SPTR) regex;

    re = pcre2_compile(
                pattern,               /* the pattern */
                PCRE2_ZERO_TERMINATED, /* indicates pattern is zero-terminated */
                0,                     /* default options */
                &errornumber,          /* for error number */
                &erroroffset,          /* for error offset */
                NULL);

    if (re == NULL)
    {
        PCRE2_UCHAR buffer[256];
        pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
        fprintf(stderr, "[E:%s] Regular expression compilation failed : %s at position %d\n", __FUNCTION__, buffer, (int32_t) erroroffset);
        exit(1);
    }

    match_data = pcre2_match_data_create_from_pattern(re, NULL);
};

/**
 * Matches a text against a regular expression that has been compiled in set().
 */
bool PERLregex::match(std::string& text)
{
    return match(const_cast<char*>(text.c_str()));
};

/**
 * Matches a text against a regular expression that has been compiled in set().
 */
bool PERLregex::match(char* text)
{
    subject = (PCRE2_SPTR) text;
    subject_length = strlen((char *)subject);

    rc = pcre2_match(
                re,                   /* the compiled pattern */
                subject,              /* the subject string */
                subject_length,       /* the length of the subject */
                0,                    /* start at offset 0 in the subject */
                0,                    /* default options */
                match_data,           /* block for storing the result */
                NULL);                /* use default match context */

    return rc>0;
};