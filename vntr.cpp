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

#include "vntr.h"

/**
 * Constructor.
 */
VNTR::VNTR(){};

/**
 * Clear object.
 */
void VNTR::clear()
{
    rid = -1;
    motif.clear();
    ru.clear();
    mlen = 0;

    exact_repeat_tract.clear();
    exact_beg1 = 0;
    exact_end1 = 0;
    exact_lflank.clear();
    exact_rflank.clear();
    exact_score = 0;
    exact_rl = 0;
    exact_no_exact_ru = 0;
    exact_total_no_ru = 0;

    fuzzy_repeat_tract.clear();
    fuzzy_beg1 = 0;
    fuzzy_end1 = 0;
    fuzzy_lflank.clear();
    fuzzy_rflank.clear();
    fuzzy_score = 0;
    fuzzy_rl = 0;
    fuzzy_no_exact_ru = 0;
    fuzzy_total_no_ru = 0;

    is_large_repeat_tract = false;
}

/**
 * Checks for equality.
 */
bool VNTR::equals(VNTR& vntr)
{
    return (rid==vntr.rid &&
            exact_beg1==vntr.exact_beg1 &&
            exact_end1==vntr.exact_end1 &&
            motif==vntr.motif);
}

/**
 * Return the string of unique bases in a motif.
 */
std::string VNTR::get_basis(std::string& motif)
{
    return get_basis(const_cast<char*>(motif.c_str()), motif.size());
}

/**
 * Return the string of unique bases in a motif.
 */
std::string VNTR::get_basis(char* motif, uint32_t n)
{
    bool bases[4] = {false, false, false, false};

//    std::cerr << "input " << motif << " " << n << "\n";

    for (uint32_t i=0; i<n; ++i)
    {
        char base = motif[i];

        if (base=='A')
        {
            bases[0] = true;
        }
        else if (base=='C')
        {
            bases[1] = true;
        }
        else if (base=='G')
        {
            bases[2] = true;
        }
        else if (base=='T')
        {
            bases[3] = true;
        }
    }

    std::string basis;
    if (bases[0]) basis.append(1, 'A');
    if (bases[1]) basis.append(1, 'C');
    if (bases[2]) basis.append(1, 'G');
    if (bases[3]) basis.append(1, 'T');

//    std::cerr << "return " << basis << "\n";

    return basis;
}

/**
 * Print object.
 */
void VNTR::print()
{
    std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    std::cerr << "VNTR Summary\n";
    std::cerr << "rid          : " << rid << "\n";
    std::cerr << "motif        : " << motif << "\n";
    std::cerr << "ru           : " << ru << "\n";
    std::cerr << "basis        : " << basis << "\n";
    std::cerr << "\n";
    std::cerr << "Exact\n";
    std::cerr << "repeat_tract                    : " << exact_repeat_tract << "\n";
    std::cerr << "position                        : [" << exact_beg1 << "," << exact_end1 << "]\n";
    std::cerr << "reference repeat unit length    : " << exact_rl << "\n";
    std::cerr << "longest allele length           : " << exact_ll << "\n";
    std::cerr << "score               : " << exact_score << "\n";
    std::cerr << "repeat units                    : " << exact_rl << "\n";
    std::cerr << "exact repeat units              : " << exact_no_exact_ru << "\n";
    std::cerr << "total no. of repeat units       : " << exact_total_no_ru << "\n";
    std::cerr << "\n";
    std::cerr << "Fuzzy\n";
    std::cerr << "repeat_tract                    : " << fuzzy_repeat_tract << "\n";
    std::cerr << "position                        : [" << fuzzy_beg1 << "," << fuzzy_end1 << "]\n";
    std::cerr << "reference repeat unit length    : " << fuzzy_rl << "\n";
    std::cerr << "longest allele length           : " << fuzzy_ll << "\n";
    std::cerr << "score               : " << fuzzy_score << "\n";
    std::cerr << "repeat units                    : " << fuzzy_rl << "\n";
    std::cerr << "exact repeat units              : " << fuzzy_no_exact_ru << "\n";
    std::cerr << "total no. of repeat units       : " << fuzzy_total_no_ru << "\n";
    std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
};