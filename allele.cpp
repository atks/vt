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

#include "allele.h"

Allele::Allele(int32_t type, int32_t diff, int32_t alen, int32_t dlen, int32_t tlen, int32_t mlen, int32_t ts)
{
    this->type = type;
    this->diff = diff;
    this->alen = alen;
    this->dlen = dlen;
    this->tlen = tlen;
    this->mlen = mlen;
    this->ts = ts;
    this->tv = mlen-ts;
    this->ins = dlen>0?1:0;
    this->del = dlen<0?1:0;
};

Allele::Allele()
{
    clear();
};

Allele::~Allele() {};

void Allele::clear()
{
    type = VT_REF;
    diff = 0;
    alen = 0;
    dlen = 0;
    tlen = 0;
    mlen = 0;
    ts = 0;
    tv = 0;
    ins = 0;
};

void Allele::print()
{
    std::cerr << "\ttype: " << type << "\n";
    std::cerr << "\tdiff: " << diff << "\n";
    std::cerr << "\talen: " << alen << "\n";
    std::cerr << "\tdlen: " << dlen << "\n";
    std::cerr << "\ttlen: " << tlen << "\n";
};