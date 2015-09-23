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

    is_large_repeat_tract = false;
}

/**
 * Checks for equality.
 */
bool VNTR::equals(VNTR& vntr)
{
    return (rid==vntr.rid &&
            rbeg1==vntr.rbeg1 &&
            rend1==vntr.rend1 &&
            motif==vntr.motif);
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
    std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    std::cerr << "VNTR Summary\n";
    std::cerr << "rid          : " << rid << "\n";
    std::cerr << "motif        : " << motif << "\n";
    std::cerr << "ru           : " << ru << "\n";
    std::cerr << "\n";
    std::cerr << "Exact\n";
    std::cerr << "repeat_tract                    : " << repeat_tract << "\n";
    std::cerr << "position                        : [" << rbeg1 << "," << rend1 << "]\n";
    std::cerr << "reference repeat unit length    : " << rl << "\n";
    std::cerr << "longest allele length           : " << ll << "\n";
    std::cerr << "motif_concordance               : " << motif_concordance << "\n";
    std::cerr << "repeat units                    : " << rl << "\n";
    std::cerr << "exact repeat units              : " << no_exact_ru << "\n";
    std::cerr << "total no. of repeat units       : " << total_no_ru << "\n";
    std::cerr << "\n";
    std::cerr << "Fuzzy\n";
    std::cerr << "repeat_tract                    : " << fuzzy_repeat_tract << "\n";
    std::cerr << "position                        : [" << fuzzy_rbeg1 << "," << fuzzy_rend1 << "]\n";
    std::cerr << "reference repeat unit length    : " << fuzzy_rl << "\n";
    std::cerr << "longest allele length           : " << fuzzy_ll << "\n";
    std::cerr << "motif_concordance               : " << fuzzy_motif_concordance << "\n";
    std::cerr << "repeat units                    : " << fuzzy_rl << "\n";
    std::cerr << "exact repeat units              : " << fuzzy_no_exact_ru << "\n";
    std::cerr << "total no. of repeat units       : " << fuzzy_total_no_ru << "\n";
    std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";


    std::string fuzzy_repeat_tract;   //repeat tract
    int32_t fuzzy_rbeg1;              //beginning of repeat tract
    int32_t fuzzy_rend1;              //end of repeat tract
    float fuzzy_rl;                   //number of repeat units on repeat tract
    float fuzzy_ll;                   //number of repeat units on longest allele
    float fuzzy_motif_concordance;    //motif concordance from hmm
    int32_t fuzzy_no_exact_ru;        //number exact repeat units from hmm
    int32_t fuzzy_total_no_ru;        //total no of repeat units from hmm
    std::string fuzzy_lflank;         //left flank
    std::string fuzzy_rflank;         //right flank
};