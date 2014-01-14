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

#ifndef GENCODE_H
#define GENCODE_H

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
#include "interval_tree.h"
#include "variant_manip.h"
#include "genome_interval.h"
#include "tbx_ordered_reader.h"

#define GC_FT_EXON 0
#define GC_FT_CDS  1
#define GC_FT_START_CODON 2
#define GC_FT_STOP_CODON 3

class GENCODERecord : public Interval
{
    public:
    std::string gene;
    int32_t feature;
    std::string chrom;
    char strand;
    int32_t frame;
    int32_t exonNo;
    bool fivePrimeConservedEssentialSpliceSite;
    bool threePrimeConservedEssentialSpliceSite;
    bool containsStartCodon;
    bool containsStopCodon;
    uint32_t level;

    GENCODERecord(std::string& _chrom, uint32_t _start, uint32_t _end, char _strand,
                  std::string& _gene, int32_t _feature, int32_t _frame, int32_t _exonNo,
                  bool _fivePrimeConservedEssentialSpliceSite, bool _threePrimeConservedEssentialSpliceSite,
                  bool _containsStartCodon, bool _containsStopCodon,
                  uint32_t _level);

    /**
     * Prints this GENCODE record to STDERR.
     */
    void print();

    /**
     * Converts feature to string.
     */
    void feature2string(int32_t feature, kstring_t *s);

    private:
};

class GENCODE
{
    public:
    std::string gencode_gtf_file;
    std::string ref_fasta_file;
    faidx_t *fai;
    std::map<std::string, IntervalTree*> CHROM;
    std::stringstream token;

    /**
     * Constructs and initialize a GENCODE object.
     */
    GENCODE(std::string& gencode_gtf_file, std::string& ref_fasta_file, std::vector<GenomeInterval>& intervals);

    /**
     * Constructs a GENCODE object.
     */
    GENCODE(std::string& gencode_gtf_file, std::string& ref_fasta_file);

    /**
     * Initialize a vector of intervals.
     */
    void initialize(std::vector<GenomeInterval>& intervals);

    /**
     * Initialize a chromosome in the GENCODE tree.
     */
    void initialize(std::string& chrom);

    /**
     * Gets overlapping intervals with chrom:start1-end1.
     */
    void search(std::string& chrom, int32_t start1, int32_t end1, std::vector<Interval*>& intervals);

    /**
     * Splits a line into a map - PERL style.
     */
    void split_gtf_attribute_field(std::map<std::string, std::string>& map, std::string& str);
};

#endif