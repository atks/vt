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

#ifndef PED_H
#define PED_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <string>
#include <iostream>
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "hts_utils.h"
#include "utils.h"

#define PED_UNKNOWN_SEX -1
#define PED_MALE 0
#define PED_FEMALE  1

class Trio
{
    public:
    std::string pedigree;
    std::string father;
    std::string mother;
    std::string child;
    int32_t father_index;
    int32_t mother_index;
    int32_t child_index;
    int32_t child_sex;    
    
    Trio(std::string pedigree, std::string father, std::string mother, std::string child, int32_t child_sex);

    private:
};

class Pedigree
{
    public:
    std::vector<Trio> trios;
    std::string ped_file;
    
    /**
     * Constructs and initialize a Pedigree object.
     */
    Pedigree(std::string& ped_file);

    private:
};

#endif