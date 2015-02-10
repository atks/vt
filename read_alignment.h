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

#ifndef READ_ALIGNMENT_H
#define READ_ALIGNMENT_H

#include "hts_utils.h"
#include "utils.h"
#include "variant_manip.h"
#include <vector>

#define M 0 //match
#define X 1 //mismath
#define I 2 //insertion
#define D 3 //deletion
#define J 4 //right overhang soft clip
#define K 5 //left overhang soft clip


/**
 * A feature structure.
 */
typedef struct
{
    int32_t type;
    int32_t start1;
    int32_t end1;
       
} feature_t;


/**
 * A read object that allows represents the alignment.
 */
class ReadAlignment
{
    public:

    char* seq;
    char* qual;
    uint32_t* cigar;
    char* md;
    
    uint32_t* cptr;
    char* mptr;
    
    /**
     * Constructor
     */
    ReadAlignment();

    /**
     * Sets the read object with a new read alignment.
     */
    void set(char* seq, char* qual, uint32_t* cigar, char* md);
       
    /**
     * Gets the next feature.
     */
    void get_next_feature(feature_t* f);

    private:
};

#endif