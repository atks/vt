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

#ifndef VNTR_CONSOLIDATOR_H
#define VNTR_CONSOLIDATOR_H

#include "hts_utils.h"
#include "utils.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "Rmath/Rmath.h"
#include "candidate_region_extractor.h"
#include "candidate_motif_picker.h"
#include "flank_detector.h"

/**
 * struct for storing sequence content.
 */
typedef struct
{
    std::string basis;
    float proportion;
} basis_proportion;

/**
 * Comparator for BCFPtr class.  Used in priority_queue; ensures that
 * records are ordered according to file order.
 */
class CompareBasisProportion
{
    public:
    bool operator()(basis_proportion& a, basis_proportion& b)
    {
        return a.proportion >= b.proportion;
    }
};

/**
 * For consolidating overlapping VNTRs.
 */
class VNTRConsolidator
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

    int32_t buffer_window_allowance;

    ////////////////
    //variant buffer
    ////////////////
    std::list<Variant *> variant_buffer; //front is most recent

    ////////////
    //filter ids
    ////////////
    char* overlap_snp;
    char* overlap_indel;
    char* overlap_vntr;
    int32_t overlap_snp_id;
    int32_t overlap_indel_id;
    int32_t overlap_vntr_id;

    /////////
    //stats//
    /////////
    int32_t no_snps;
    int32_t no_indels;
    int32_t no_vntrs;
    int32_t no_other_variants;

    int32_t no_total_variants;
    int32_t no_overlap_vntrs;
    std::vector<int32_t> overlapping_vntr_hist;
    int32_t no_dropped_vntrs;

    //vntrs that do not overlap with any other vntrs
    int32_t no_isolated_vntrs;

    //vntrs that overlap but have consistent repeat units
    int32_t no_clustered_consistent_ru_vntrs;
    int32_t no_merged_consistent_ru_vntrs;

    //vntrs that overlap but have consistent bases
    int32_t no_clustered_consistent_basis_vntrs;
    int32_t no_merged_consistent_basis_vntrs;

    //vntrs that overlap but have inconsistent repeat units and bases
    int32_t no_clustered_inconsistent_ru_basis_vntrs;

    //for storing basis information
    std::priority_queue<basis_proportion, std::vector<basis_proportion>, CompareBasisProportion> ordered_basis;

    ///////
    //tools
    ///////
    ReferenceSequence *refseq;
    CandidateMotifPicker* cmp;
    FlankDetector* fd;

    /**
     * Constructor.
     */
    VNTRConsolidator(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& ref_fasta_file, bool debug=false);

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
     * Detects VNTR overlapping class.
     */
    void detect_VNTR_overlapping_class(Variant* variant, bool& consistent_repeat_units, bool& consistent_bases);

    /**
     * Merge overlapping VNTRs with a consistent repeat unit.
     */
    void merge_consistent_ru_overlapping_VNTR(Variant* variant);

    /**
     * Merge overlapping VNTRs with a consistent basis.
     */
    void merge_consistent_basis_overlapping_VNTR(Variant* variant);
    
    /**.
     * Flush variant buffer.
     */
    void flush_variant_buffer();

    /**
     * Close.
     */
    void close();

    private:
};

#endif