/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#ifndef REPEAT_REGION_H
#define REPEAT_REGION_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <queue>
#include <list>
#include "hts_utils.h"

#define REFERENCE     0
#define ALLELE_EXACT  1
#define ALLELE_FUZZY  2

class RepeatRegion
{
    public:
    uint32_t beg1;
    uint32_t end1;
    std::string ref;

    /**
     * Constructor.
     */
    RepeatRegion() {};

    /**
     * Constructor.
     */
    RepeatRegion(uint32_t beg1, char* ref)
    {
        this->beg1 = beg1;
        this->ref.assign(ref);
        this->end1 = beg1 + this->ref.size() - 1;
    };

    /**
     * Initialize ReferenceRegion.
     */
    void initialize(uint32_t beg1, char* ref)
    {
        this->beg1 = beg1;
        this->ref.assign(ref);
        this->end1 = beg1 + this->ref.size() - 1;
    };
};

#endif