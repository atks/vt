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

#ifndef CANDIDATE_REGION_EXTRACTOR_H
#define CANDIDATE_REGION_EXTRACTOR_H

#include "hts_utils.h"
#include "utils.h"
#include "ahmm.h"
#include "variant_manip.h"
#include "vntr.h"
#include "reference_sequence.h"

//modes for extracting candidate repeat region
#define REFERENCE                   0
#define EXACT_LEFT_RIGHT_ALIGNMENT  1
#define FUZZY_LEFT_RIGHT_ALIGNMENT  2

/**
 * Extract a candidate repeat region without knowledge of repeat unit.
 */
class CandidateRegionExtractor
{
    public:

    //optiona
    bool debug;

    ///////
    //tools
    ///////
    VariantManip *vm;
    ReferenceSequence* rs;

    /**
     * Constructor.
     */
    CandidateRegionExtractor(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~CandidateRegionExtractor();

    /**
     * Annotates VNTR characteristics.
     * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE
     */
    void annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant, std::string mode);

    /**
     * Pick candidate region.
     *
     * @mode - REFERENCE     use refence field
     *       - ALLELE_EXACT  by exact alignment
     *       - ALLELE_FUZZY  by fuzzy alignment
     */
    void pick_candidate_region(Variant& variant,  int32_t mode, int32_t amode);

    /**
     * Chooses a phase of the motif that is appropriate for the alignment
     */
    std::string choose_repeat_unit(std::string& ref, std::string& motif);

    /**
     * Trim alleles.
     */
    void trim(int32_t& pos1, std::string& ref, std::string& alt);

    /**
     * Checks if a vntr is a homopolymer.
     */
    bool is_homopolymer(bcf_hdr_t* h, bcf1_t* v);

    /**
     * Extract region to for motif discovery.
     */
    void extract_regions_by_exact_alignment(Variant& variant);

    /**
     * Left align alleles.
     */
    void left_align(std::string& chrom, int32_t& pos1, std::string& ref, std::string& alt);

    /**
     * Right align alleles.
     */
    void right_align(std::string& chrom, int32_t& pos1, std::string& ref, std::string& alt);

    /**
     * Extract reference sequence region for motif discovery in a fuzzy fashion.
     */
    void extract_regions_by_fuzzy_alignment(Variant& variant);

    /**
     * Fuzzy left align alleles allowing for mismatches and indels defined by penalty.
     *
     * @chrom   - chromosome
     * @pos1    - 1 based position
     * @ref     - reference sequence
     * @alt     - alternative sequence
     * @penalty - mismatch/indels allowed
     *
     * Returns left aligned position.
     */
    uint32_t fuzzy_left_align(std::string& chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty);

    /**
     * Fuzzy right align alleles allowing for mismatches and indels defined by penalty.
     *
     * @chrom   - chromosome
     * @pos1    - 1 based position
     * @ref     - reference sequence
     * @alt     - alternative sequence
     * @penalty - mismatch/indels allowed
     *
     * Returns right aligned position.
     */
    uint32_t fuzzy_right_align(std::string& chrom, int32_t pos1, std::string ref, std::string alt, uint32_t penalty);

};

#endif