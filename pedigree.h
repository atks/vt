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

#include "utils.h"
#include "hts_utils.h"

#define PED_MALE     0
#define PED_FEMALE   1
#define PED_OTHER    2

class PEDRecord
{
    public:
    std::string pedigree;
    std::vector<std::string> individual;
    std::string father;
    std::string mother;
    int32_t individual_sex;

    /**
     * Constructs and initialize a PEDRecord object.
     */
    PEDRecord(std::string pedigree,
              std::vector<std::string>& individual,
              std::string father,
              std::string mother,
              int32_t individual_sex);

    /**
     * Returns true if record is representative of a trio, false otherwise.
     */
    bool is_trio();

    /**
     * Returns true if record contains duplicates, false otherwise.
     */
    bool is_duplicated();

    /**
     * Prints the contents of this record.
     */
    void print();

    private:
};

class Pedigree
{
    public:
    std::vector<PEDRecord> recs;
    std::string ped_file;

    /**
     * Constructs and initialize a Pedigree object.
     */
    Pedigree(std::string& ped_file);

    private:
};

#endif