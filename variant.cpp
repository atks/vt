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

#include "variant.h"

/**
 * Constructor.
 */
Variant::Variant(bcf1_t* v)
{
    this->v = v;
}

/**
 * Constructor.
 */
Variant::Variant()
{
    type = VT_REF;
    alleles.clear();
}

/**
 * Destructor.
 */
Variant::~Variant() {};

/**
 * Clears variant information.
 */
void Variant::clear()
{
    type = VT_REF;
    alleles.clear();
    vntr.clear();
};

/**
 * Prints variant information.
 */
void Variant::print()
{
    std::cerr << "type : " << vtype2string(type) << "\n";
    std::cerr << "motif: " << vntr.motif << "\n";
    std::cerr << "rlen : " << vntr.motif.size() << "\n";

    for (int32_t i=0; i<alleles.size(); ++i)
    {
        std::cerr << "\tallele: " << i << "\n";
        std::cerr << "\t  type: " << vtype2string(alleles[i].type) << "\n";
        std::cerr << "\t  diff: " << alleles[i].diff << "\n";
        std::cerr << "\t  alen: " << alleles[i].alen << "\n";
        std::cerr << "\t  dlen: " << alleles[i].dlen << "\n";
    }
};

/**
 * Gets a string representation of the underlying VNTR.
 */
void Variant::get_vntr_string(kstring_t* s)
{
    s->l = 0;
    kputs(chrom.c_str(), s);
    kputc(':', s);
    kputw(vntr.rbeg1, s);
    kputc(':', s);
    kputs(vntr.repeat_tract.c_str(), s);
    kputc(':', s);
    kputs("<VNTR>", s);
};

/**
 * Gets a string representation of the underlying VNTR.
 */
void Variant::get_fuzzy_vntr_string(kstring_t* s)
{
    s->l = 0;
    kputs(chrom.c_str(), s);
    kputc(':', s);
    kputw(vntr.fuzzy_rbeg1, s);
    kputc(':', s);
    kputs(vntr.fuzzy_repeat_tract.c_str(), s);
    kputc(':', s);
    kputs("<VNTR>", s);
};

/**
 * Converts VTYPE to string.
 */
std::string Variant::vtype2string(int32_t VTYPE)
{
    std::string s;

    if (!VTYPE)
    {
        s += (s.size()==0) ? "" : "/";
        s += "REF";
    }

    if (VTYPE & VT_SNP)
    {
        s += (s.size()==0) ? "" : "/";
        s += "SNP";
    }

    if (VTYPE & VT_MNP)
    {
        s += (s.size()==0) ? "" : "/";
        s += "MNP";
    }

    if (VTYPE & VT_INDEL)
    {
        s += (s.size()==0) ? "" : "/";
        s += "INDEL";
    }

    if (VTYPE & VT_CLUMPED)
    {
        s += (s.size()==0) ? "" : "/";
        s += "CLUMPED";
    }

    if (VTYPE & VT_VNTR)
    {
        s += (s.size()==0) ? "" : "/";
        s += "VNTR";
    }

    if (VTYPE & VT_SV)
    {
        s += (s.size()==0) ? "" : "/";
        s += "SV";
    }

    return s;
}