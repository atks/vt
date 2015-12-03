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
    motif_score = 0;
    mlen = 0;

    exact_repeat_tract.clear();
    exact_rbeg1 = 0;
    exact_rend1 = 0;
    exact_lflank.clear();
    exact_rflank.clear();
    exact_motif_concordance = 0;
    exact_rl = 0;
    exact_no_exact_ru = 0;
    exact_total_no_ru = 0;

    fuzzy_repeat_tract.clear();
    fuzzy_rbeg1 = 0;
    fuzzy_rend1 = 0;
    fuzzy_lflank.clear();
    fuzzy_rflank.clear();
    fuzzy_motif_concordance = 0;
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
            exact_rbeg1==vntr.exact_rbeg1 &&
            exact_rend1==vntr.exact_rend1 &&
            motif==vntr.motif);
}

/**
 * Return the string of unique bases in a motif.
 */
std::string VNTR::get_basis(std::string& motif)
{
    bool bases[4] = {false, false, false, false};
    
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        char base = motif.at(i);
        
        if (base<=67)
        {
            if (base==65)
            {
                if (!bases[0]) 
                {
                    bases[0] = true;
                }
            }  
            else
            {
                if (!bases[1]) 
                {
                    bases[1] = true;
                }
            }  
        }   
        else
        {
            if (base==71)
            {
                if (!bases[2]) 
                {
                    bases[2] = true;
                }
            }  
            else
            {
                if (!bases[3]) 
                {
                    bases[3] = true;
                }
            } 
        } 
    }
    
    std::string basis;
    if (bases[0]) basis.append(1, 'A');
    if (bases[1]) basis.append(1, 'C');
    if (bases[2]) basis.append(1, 'G');
    if (bases[3]) basis.append(1, 'T');
        
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
    std::cerr << "position                        : [" << exact_rbeg1 << "," << exact_rend1 << "]\n";
    std::cerr << "reference repeat unit length    : " << exact_rl << "\n";
    std::cerr << "longest allele length           : " << exact_ll << "\n";
    std::cerr << "motif_concordance               : " << exact_motif_concordance << "\n";
    std::cerr << "repeat units                    : " << exact_rl << "\n";
    std::cerr << "exact repeat units              : " << exact_no_exact_ru << "\n";
    std::cerr << "total no. of repeat units       : " << exact_total_no_ru << "\n";
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