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
    repeat_tract.clear();
    rbeg1 = 0;
    rend1 = 0;
    lflank.clear();
    rflank.clear();

    motif.clear();
    ru.clear();
    mlen = 0;

    motif_score = 0;
    motif_concordance = 0;
    rl = 0;
    no_exact_ru = 0;
    total_no_ru = 0;
}

/**
 * Get VNTR representation in string format.
 */
void VNTR::get_vntr_allele_string(std::string& var)
{
    var.clear();
    var.append(repeat_tract.c_str());
    var.append(1, ':');
    var.append("<VNTR>");
}

/**
 * Print object.
 */
void VNTR::print()
{
    std::cerr << "++++++++++++++++++++++++++++\n";
    std::cerr << "VNTR Summary\n";
    std::cerr << "repeat_tract : " << repeat_tract << "\n";
    std::cerr << "pos1         : " << rbeg1 << "\n";
    std::cerr << "motif        : " << motif << "\n";
    std::cerr << "ru           : " << ru << "\n";
    std::cerr << "rl           : " << rl << "\n";
    std::cerr << "++++++++++++++++++++++++++++\n";
};