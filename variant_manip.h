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
#include "htslib/vcf.h"
#include "hts_utils.h"
#include "variant.h"
#include "allele.h"

//int BitCount(unsigned int u)
//{
//     unsigned int uCount;
//
//#ifndef countbit32
//#define countbit32(x) (y = u - ((u >> 1) & 033333333333) - ((u >> 2) & 011111111111); ((y + (y >> 3)) & 030707070707) % 63);
//
//#endif

#define VT_REF      0   //dlen==0 && diff==0
#define VT_SNP      1   //min(rlen,alen)==1 && diff==1
#define VT_MNP      2   //min(rlen,alen)==diff
#define VT_INDEL    4   //diff!=0 && (rlen==1 || alen==1)
#define VT_CLUMPED  8   //all others
#define VT_SV       16  //structural variant tags

/**
 * Methods for manipulating variants
 */
class VariantManip
{
    public:
    faidx_t *fai;
    bool reference_present;

    /**
     * Constructor.
     *
     * @ref_fasta_file reference sequence FASTA file.
     */
    VariantManip(std::string ref_fasta_file);

    /**
     * Constructor.
     */
    VariantManip();

    /**
     * Detects near by STRs.
     */
    bool detect_str(bcf_hdr_t *h, bcf1_t *v, Variant& variant);

    /**
     * Detects near by STRs.
     */
    bool detect_str(const char* chrom, uint32_t pos1, Variant& variant);

    /**
     * Converts VTYPE to string.
     */
    std::string vtype2string(int32_t VTYPE);

    /**
     * Converts VTYPE to string.
     */
    void vtype2string(int32_t vtype, kstring_t *s);

    /**
     * Inspects variants.
     */
    int32_t inspect(char** allele, int32_t n_allele);

    /**
     * Classifies variants.
     */
    int32_t classify_variant(bcf_hdr_t *h, bcf1_t *v, Variant& variant, bool in_situ_left_trimming = true);

    /**
     * Classifies variants.
     */
    int32_t classify_variant(bcf_hdr_t *h, bcf1_t *v);

    /**
     * Classifies variants.
     */
    int32_t classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele, Variant& variant, bool in_situ_left_trimming = true);

    /**
     * Classifies variants.
     */
    int32_t classify_variant(const char* chrom, uint32_t pos1, char** allele, int32_t n_allele);

    /**
     * Checks if a variant is normalized.
     */
    bool is_normalized(char** alleles, int32_t n_allele);

    /**
     * Right trims or left extend a variant.
     */
    void right_trim_or_left_extend(std::vector<std::string>& alleles, uint32_t& pos1, const char* chrom, uint32_t& left_extended, uint32_t& right_trimmed);

    /**
     * Left trims a variant with unnecesary nucleotides.
     */
    void left_trim(std::vector<std::string>& alleles, uint32_t& pos1, uint32_t& left_trimmed) ;

    /**
     * Generates a probing haplotype with flanks around the variant of interest.
     */
    void generate_probes(const char* chrom,
                        int32_t pos1, uint32_t probeDiff,
                        std::vector<std::string>& alleles, //store alleles
                        std::vector<std::string>& probes, //store probes
                        uint32_t min_flank_length,
                        int32_t& preambleLength); //store preamble length

    private:

    /**
     * Recursive helper method for generateProbes.
     */
    void generate_probes(const char* chrom,
                        int32_t pos1,
                        uint32_t flankLength,
                        uint32_t& currentDiff,
                        uint32_t& length,
                        uint32_t gald,
                        std::vector<uint32_t>& diff,
                        std::vector<std::string>& alleles,
                        std::vector<std::string>& probes);

};

#endif
