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

#include "repeat_region.h"

/**
 * Constructor.
 */
RepeatRegion::RepeatRegion() {};

/**
 * Constructor.
 */
RepeatRegion::RepeatRegion(uint32_t beg1, char* ref)
{
    initialize(beg1, ref);
};

/**
 * Initialize RepeatRegion.
 */
void RepeatRegion::initialize(uint32_t beg1, char* ref)
{
    this->beg1 = beg1;
    this->ref.assign(ref);
    this->end1 = beg1 + this->ref.size() - 1;
};

/**
 * Clears RepeatRegion.
 */
void RepeatRegion::clear()
{
    this->beg1 = 0;
    this->end1 = 0;
    this->ref.clear();
};