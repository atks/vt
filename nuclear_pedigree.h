/* The MIT License

   Copyright (c) 2015 Hyun Min Kang <atks@umich.edu>

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

#ifndef NUCLEAR_PEDIGREE_H
#define NUCLEAR_PEDIGREE_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "hts_utils.h"
#include "utils.h"

/**
 * A class for storing nuclear family allowing duplicates
 *
 */
class NuclearFamilySample;
class NuclearFamilyPerson;
class NuclearFamily;
class NuclearPedigree;

class NuclearFamilySample {  // Represents each sequenced/genotyped sample.
 public:
  std::string sampleID; // sample ID
  int index;      // sample index in the BCF file
  NuclearFamilySample(const char* id) : sampleID(id), index(-1) {};
};

class NuclearFamilyPerson {  // represents each individuals
 public:
  std::vector<NuclearFamilySample*> samples; // could be multiple samples if duplicated
  int sex;  // 0 : unknown, 1 : male, 2 : female
  NuclearFamilyPerson(int _sex) : sex(_sex) {}
  const char* getUID() { return samples.front()->sampleID.c_str(); }
  bool hasSample(const char* smID) {
    for(size_t i=0; i < samples.size(); ++i) {
      if ( samples[i]->sampleID.compare(smID) == 0 )
	return true;
    }
    return false;
  }
  bool addSample(NuclearFamilySample* pSample) { samples.push_back(pSample); }
  int removeSamplesWithoutIndex(std::map<std::string,NuclearFamilySample*>& smIDmap) {
    int nRemoved = 0;
    for(int i=(int)samples.size()-1; i >= 0; --i) {
      if ( samples[i]->index < 0 ) {
	smIDmap.erase(samples[i]->sampleID);
	delete(samples[i]);               // call desctuctor
	samples.erase(samples.begin()+i); // remove the pointer from the vector
	++nRemoved;
      }
    }
    return nRemoved;
  }
};

class NuclearFamily {
 public:
  std::string famID;
  NuclearFamilyPerson* pDad;
  NuclearFamilyPerson* pMom;
  std::vector<NuclearFamilyPerson*> pKids;

  NuclearFamily(const char* fid) : famID(fid), pDad(NULL), pMom(NULL) {}
};

class NuclearPedigree
{
 public:
  std::map<std::string, NuclearFamily*> famIDmap;
  std::map<std::string, NuclearFamilySample*> smIDmap;

  NuclearPedigree(const char* pedFile);
  void addPerson(const char* famID, const std::vector<std::string>& smIDs, const char* dadID, const char* momID, int sex);
  bool setSampleIndex(const char* smID, int index);
  int numPeople();
  int numSamplesWithIndex();
  int removeSamplesWithoutIndex();
};

#endif
