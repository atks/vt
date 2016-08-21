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
#include "filter.h"

#define LAI2003         0
#define KELKAR2008      1
#define FONDON2012      2
#define ANANDA2013      3
#define WILLEMS2014     4
#define MONTGOMERY2014  5
#define TAN_KANG2015    6
#define EXACT_VNTR      7
#define FUZZY_VNTR      8

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
    int32_t no_samples;

    int32_t* gts;
    int32_t buffer_window_allowance;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    ////////////////
    //variant buffer
    ////////////////
    std::list<Variant *> vbuffer; //front is most recent

    /////////////
    //INFO fields
    /////////////
    std::string END;
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
    std::string ASSOCIATED_INDEL;

    int32_t vntr_classification;

    /////////
    //stats//
    /////////
    int32_t no_variants;
    int32_t no_indels;
    int32_t no_added_vntrs;
    int32_t no_duplicate_vntrs;

    ///////
    //tools
    ///////
    ReferenceSequence *refseq;

    /**
     * Constructor.
     */
    VNTRExtractor(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& fexp,  int32_t vntr_classification_code, std::string& ref_fasta_file);

    /**
     * Update distribution of overlapping VNTRs
     */
    void update_overlapping_vntr_hist(int32_t no_overlapping_vntrs);

    /**
     * Inserts a Variant record.
     */
    void insert(Variant* variant);

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

    /**
     * Process.
     */
    void process();

    /**.
     * Flush buffer.
     */
    void flush(Variant* var=NULL);

    /**.
     * Close files.
     */
    void close();

    /**
     * Process exiting variant.
     */
    void process_exit(Variant* var);

    /**
     * Copies exact VNTR features to finalized features.
     */
    void copy_exact_vntr_features_to_final_vntr_features(VNTR& vntr);

    /**
     * Copies fuzzy VNTR features to finalized features.
     */
    void copy_fuzzy_vntr_features_to_final_vntr_features(VNTR& vntr);

    /**
     * Creates a VNTR record based on classification schema.
     */
    void create_and_insert_vntr(Variant& nvar);

    /**
     * Converts VNTR classification to string
     */
    std::string vntr_code_2_str(int32_t code); 

    private:
};

#endif
