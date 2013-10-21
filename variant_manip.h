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

#ifndef VARIANT_MANIP_H
#define VARIANT_MANIP_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <sstream>
#include "htslib/faidx.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "tclap/CmdLine.h"
#include "tclap/Arg.h"
#include "hts_utils.h"

#define VT_REF   0
#define VT_SNP   1
#define VT_MNP   2
#define VT_INDEL   4
#define VT_INSERTION  8
#define VT_DELETION   16
#define VT_STR 32
#define VT_EXACT_STR 64
#define VT_INEXACT_STR 128
#define VT_COMPLEX 256
#define VT_SV 512
#define VT_CR 1024

/**
 * Methods for manipulating variants
 * 1. variant identification
 * 2. left trimming
 * 3. left alignment and right trimming
 * 4. haplotype contruction
 */
class VariantManip
{
    public:
        
    faidx_t *fai;
    
    VariantManip(std::string ref_fasta_file);

    /**
     * Detects near by STRs
     */
    bool detect_str(const char* chrom, uint32_t pos1, std::string& motif, uint32_t& tlen);

    /**
     * Converts VTYPE to string
     */
    std::string vtype2string(int32_t VTYPE);

    /**
     * Classifies Indels into the following categories:
     * 1. Homopolymers
     * 2. Dinucleotides
     * 3. Trinucleotides
     */
    int32_t classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele, std::string& motif, uint32_t& tlen);

    /**
     * Left trims a variant of unnecesary nucleotides.
     */
    void left_trim(std::vector<std::string>& alleles, uint32_t& pos1, uint32_t& left_trimmed) ;

    /**
     * Left aligns a variant.
     */
    void left_align(std::vector<std::string>& alleles, uint32_t& pos1, const char* chrom, uint32_t& leftAligned, uint32_t& right_trimmed);

    private:
};

#endif
