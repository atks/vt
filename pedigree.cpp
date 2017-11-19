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

#include "pedigree.h"

/**
 * Constructs and initialize a PEDRecord object.
 */
PEDRecord::PEDRecord(std::string pedigree, std::vector<std::string>& individual, 
                     std::string father, std::string mother, int32_t individual_sex)
{
    this->pedigree = pedigree;
    this->individual = individual;
    this->father = father;
    this->mother = mother;
    this->individual_sex = individual_sex;
};

/**
 * Returns true if record is representative of a trio, false otherwise.
 */
bool PEDRecord::is_trio()
{
    return this->father!="" && this->mother!= "";
};

/**
 * Returns true if record contains duplicates, false otherwise.
 */
bool PEDRecord::is_duplicated()
{
    return this->individual.size()>1;
};

/**
 * Returns true if individual is male.
 */
bool PEDRecord::is_male()
{
    return individual_sex==PED_MALE;
};

/**
 * Returns true if individual is male.
 */
bool PEDRecord::is_female()
{
    return individual_sex==PED_FEMALE;
};

/**
 * Prints the contents of this record.
 */
void PEDRecord::print()
{
    fprintf(stderr, "================================\n");
    fprintf(stderr, "Family ID      : %s\n", pedigree.c_str());
    fprintf(stderr, "\n");
    fprintf(stderr, "Individual ID  : %s (%zd)\n", individual[0].c_str(), individual.size());
    fprintf(stderr, "Father ID      : %s\n", father.c_str());
    fprintf(stderr, "Mother ID      : %s\n", mother.c_str());
    fprintf(stderr, "Sex            : %s\n", (individual_sex==1?"male":(individual_sex==2?"female":"other")));
    fprintf(stderr, "================================\n");
};
    
/**
 * Constructs and initialize a Pedigree object.
 */
Pedigree::Pedigree(std::string& ped_file)
{
    htsFile *hts = hts_open(ped_file.c_str(), "r");
    kstring_t s = {0,0,0};
    std::vector<std::string> vec;
    while (hts_getline(hts, '\n', &s)>=0)
    {
        if (s.s[0] == '#')
            continue;

        std::string line(s.s);
        split(vec, "\t ", line);

        std::string pedigree = vec[0];
        std::string individual = vec[1];
        std::string father = vec[2];
        std::string mother = vec[3];
        int32_t individual_sex = PED_OTHER;
        
        vec[4] = to_lower(vec[4]);
        
        if (vec[4]=="male" || vec[4]=="1")
        {
            individual_sex = PED_MALE;
        }
        else if (vec[4] == "female"|| vec[4]=="2")
        {
            individual_sex = PED_FEMALE;
        }
        else if (vec[4] == "other")
        {
            individual_sex = PED_OTHER;
        }
        
        vec.resize(0);
        split(vec, ",", individual);
        
        PEDRecord rec(pedigree, vec, father, mother, individual_sex);
        recs.push_back(rec);
    }
    hts_close(hts);
    if (s.m) free(s.s);
};