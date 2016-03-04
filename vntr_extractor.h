/* The MIT License

   Copyright (c) 2016 Adrian Tan <atks@umich.edu>

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

#ifndef VNTR_EXTRACTOR_H
#define VNTR_EXTRACTOR_H

#include "hts_utils.h"
#include "utils.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "Rmath/Rmath.h"
#include "candidate_region_extractor.h"
#include "candidate_motif_picker.h"
#include "flank_detector.h"

#define LAI2003      0 
#define KELKAR2008   1
#define FONDON2012   2
#define ANANDA2013   3
#define WILLEMS2014  4
#define TAN_KANG2015 5
#define EXACT_VNTR   6
#define FUZZY_VNTR   7

/**
 * For consolidating overlapping VNTRs.
 */
class VNTRExtractor
{
    public:

    bool debug;

    ///////
    //i/o//
    ///////
    std::string input_vcf_file;
    std::string output_vcf_file;
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    int32_t* gts;
    int32_t buffer_window_allowance;

    ////////////////
    //variant buffer
    ////////////////
    std::list<Variant *> variant_buffer; //front is most recent

    /////////////
    //INFO fields
    /////////////
    std::string MOTIF;
    std::string RU;
    std::string BASIS;
    std::string MLEN;
    std::string BLEN;
    std::string REPEAT_TRACT;
    std::string COMP;
    std::string ENTROPY;
    std::string ENTROPY2;
    std::string KL_DIVERGENCE;
    std::string KL_DIVERGENCE2;
    std::string RL;
    std::string LL;
    std::string RU_COUNTS;
    std::string SCORE;
    std::string TRF_SCORE;

    /////////
    //stats//
    /////////
    int32_t no_variants;
    int32_t no_added_vntrs;
   
    ///////
    //tools
    ///////
    ReferenceSequence *refseq;
 
    /**
     * Constructor.
     */
    VNTRExtractor(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& ref_fasta_file);

    /**
     * Update distribution of overlapping VNTRs
     */
    void update_overlapping_vntr_hist(int32_t no_overlapping_vntrs);

    /**
     * Inserts a Variant record.
     */
    void insert_variant_record_into_buffer(Variant* variant);

    /**
     * Flush variant buffer.
     */
    void flush_variant_buffer(Variant* var);

    /**
     * Compute purity by sequence content.
     */
    float compute_purity_by_sequence_content(char* repeat_tract, char* motif);

    /**
     * Consolidate multiallelic variant based on associated biallelic records
     * stored in vs.  Updates v which is to be the consolidated multiallelic
     * variant.
     *
     * returns true if the multiallelic variant is good to go.
     */
    bool consolidate_multiple_overlapping_vntrs(Variant* variant);

    /**
     * Detects a a consistent basis motif in a chain of overlapping VNTRs.
     */
    bool detect_consistent_motifs(Variant* variant);

    /**.
     * Flush variant buffer.
     */
    void flush_variant_buffer();

    /**
     * Creates aVariant object representing a VNTR record.
     */
    Variant* create_vntr_record(bcf_hdr_t* h, bcf1_t *v);

    /**
     * Close.
     */
    void close();

    private:
};

#endif
