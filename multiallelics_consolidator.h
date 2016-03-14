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

#ifndef MULTIALLELICS_CONSOLIDATOR_H
#define MULTIALLELICS_CONSOLIDATOR_H

#include "hts_utils.h"
#include "utils.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "Rmath/Rmath.h"
#include "candidate_region_extractor.h"
#include "candidate_motif_picker.h"
#include "flank_detector.h"
#include "filter.h"

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
class MultiallelicsConsolidator
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
    std::string INVOLVED_MULTIALLEIC_VARIANT;

    /////////
    //stats//
    /////////
    int32_t no_variants;
    int32_t no_added_multiallelics;

    ///////
    //tools
    ///////
    ReferenceSequence *rs;

    /**
     * Constructor.
     */
    MultiallelicsConsolidator(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& fexp, std::string& ref_fasta_file);

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
     * Process overlapping variant.
     */
    void process_overlap(Variant& nvar, Variant&cvar);

    /**
     * Creates or updates a multiallelic.
     */
    Variant* create_or_update_multiallelic(Variant& nvar, Variant& cvar);

    /**
     * Updates a multiallelic for printing.
     */
    void update_multiallelic_for_printing(Variant& var);

    private:
};

#endif
