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


/**
 * Constructor for explicit sequence variant.
 */
Allele::Allele(int32_t type, int32_t diff, int32_t alen, int32_t dlen, int32_t mlen, int32_t ts, int32_t tv)
{
    this->type = type;
    this->diff = diff;
    this->alen = alen;
    this->dlen = dlen;
    this->mlen = mlen;
    this->ts = ts;
    this->tv = tv;
    this->ins = dlen>0?1:0;
    this->del = dlen<0?1:0;
};

/**
 * Constructor for VNTR.
 */
Allele::Allele(int32_t type)
{
    this->type = type;
    this->diff = 0;
    this->alen = 0;
    this->dlen = 0;
    this->mlen = 0;
    this->ts = 0;
    this->tv = 0;
    this->ins = 0;
    this->del = 0;
    this->sv_type = "";
}

/**
 * Constructor for SV.
 */
Allele::Allele(int32_t type, std::string& sv_type)
{
    this->type = type;
    this->diff = 0;
    this->alen = 0;
    this->dlen = 0;
    this->mlen = 0;
    this->ts = 0;
    this->tv = 0;
    this->ins = 0;
    this->del = 0;
    this->sv_type = reduce_sv_type(sv_type);
};

/**
 * Destructor.
 */
Allele::~Allele() {};

/**
 * Special dictionary for some reserve types.
 * CN\d+ be CNV and VN\d+ be VNTR.
 */
std::string Allele::reduce_sv_type(std::string& sv_type)
{
    size_t len = sv_type.size();
    if (len>=5)
    {
        if (sv_type.at(0)=='<' && sv_type.at(1)=='C' && sv_type.at(2)=='N' && sv_type.at(len-1)=='>' )
        {
            for (size_t i=3; i<len-1; ++i)
            {
                if (sv_type.at(i)<'0' || sv_type.at(i)>'9')
                {
                    return sv_type;
                }
            }

            return "<CNV>";
        }
        else if (sv_type.at(0)=='<' && sv_type.at(1)=='V' && sv_type.at(2)=='N' && sv_type.at(len-1)=='>' )
        {
            for (size_t i=3; i<len-1; ++i)
            {
                if (sv_type.at(i)<'0' || sv_type.at(i)>'9')
                {
                    return sv_type;
                }
            }

            return "<VNTR>";
        }
    }

    return sv_type;
};

/**
 * Clear variables.
 */
void Allele::clear()
{
    type = VT_REF;
    diff = 0;
    alen = 0;
    dlen = 0;
    mlen = 0;
    ts = 0;
    tv = 0;
    ins = 0;
    sv_type = "";
};

/**
 * Print allele.
 */
void Allele::print()
{
    std::cerr << "\ttype: " << type << "\n";
    std::cerr << "\tdiff: " << diff << "\n";
    std::cerr << "\talen: " << alen << "\n";
    std::cerr << "\tdlen: " << dlen << "\n";
    std::cerr << "\tsv_type: " << sv_type << "\n";
};