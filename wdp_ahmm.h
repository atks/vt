/* The MIT License

   Copyright (c) 2017 Adrian Tan <atks@umich.edu>

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
#ifndef WDP_AHMM_H
#define WDP_AHMM_H

#include "hts_utils.h"
#include "utils.h"
#include "log_tool.h"

#define NSTATES 6

/**
 * Parameters for WDP_AHMM
 */
class WDP_AHMMParameters
{
    public: 
        
    float delta;
    float epsilon;
    float tau;
    float eta;
    float mismatch_penalty;
        
    WDP_AHMMParameters()
    {
        delta = 0.0001;
        epsilon = 0.001;
        tau = 0.01;
        eta = 0.01;
        mismatch_penalty = 1;
    };
};

/**
 * Wrap around Dynamic Programming implementation of AHMM
 */
class WDP_AHMM
{
    public:

    int32_t max_len;

    const char* read;
    const char* qual;

    //model variables
    char *motif;
    int32_t mlen;
    int32_t rlen;
    
    /*result variables*/
    int32_t motif_start[2], motif_end[2];
    int32_t motif_count, exact_motif_count, motif_m, motif_xid;
    float frac_no_repeats;
    int32_t* motif_discordance;
    float motif_concordance, maxLogOdds;
    int32_t trf_score;
    
    WDP_AHMMParameters par;

    //for track intermediate scores during Viterbi algorithm
    float max_score;
    int32_t max_track;

    //for optimal alignment
    bool optimal_path_traced;
    float optimal_score;
    int32_t optimal_state;
    int32_t optimal_track;
    int32_t optimal_probe_len;
    int32_t *optimal_path;     // for storage
    int32_t *optimal_path_ptr; //just a pointer
    int32_t optimal_path_len;

    float T[NSTATES][NSTATES];

    float **V;
    int32_t **U;

    LogTool *lt;

    bool debug;

    /**
     * Constructor.
     */
    WDP_AHMM(bool debug=false);

    /**
     * Constructor.
     */
    WDP_AHMM(LogTool *lt, bool debug=false);

    /**
     * Destructor.
     */
    ~WDP_AHMM();

    /**
     * Initializes object, helper function for constructor.
     */
    void initialize();

    /**
     * Initializes objects for constructor.
     */
    void initialize_structures();
        
    /**
     * Initialize transition matrix based on parameters.
     */
    void initialize_T();
    
    /**
     * Initializes U and V.
     */
    void initialize_UV();

    /**
     * Sets a model.
     */
    void set_model(const char* motif);

    /**
     * Sets delta.
     */
    void set_delta(float delta);

    /**
     * Sets epsilon.
     */
    void set_epsilon(float epsilon);
    
    /**
     * Sets tau.
     */
    void set_tau(float tau);
    
    /**
     * Sets eta.
     */
    void set_eta(float eta);

    /**
     * Sets mismatch penalty.
     */
    void set_mismatch_penalty(float mismatch_penalty);

    /**
     * Sets debug.
     */
    void set_debug(bool debug);
    
    /**
     * Get motif start position for model.
     */
    int32_t get_motif_model_spos1();

    /**
     * Get motif end position for model.
     */
    int32_t get_motif_model_epos1();

    /**
     * Get motif start position for read.
     */
    int32_t get_motif_read_spos1();

    /**
     * Get motif end position for read.
     */
    int32_t get_motif_read_epos1();
    
    /**
     * Get motif concordance.
     */
    float get_motif_concordance();
    
    /**
     * Get exact motif count.
     */
    uint32_t get_exact_motif_count();
    
    /**
     * Get motif count.
     */
    uint32_t get_motif_count();
    
    /**
     * Align and compute genotype likelihood.
     */
    void align(const char* y, const char* qual=NULL);

    /**
     * Trace path after alignment.
     */
    void trace_path();

    /**
     * Compute log10 emission odds based on equal error probability distribution contrasted against log10(1/16).
     */
    float log10_emission_odds(char read_base, char probe_base, uint32_t pl, float mismatch_penalty);

    /**
     * Compute log10 emission odds based on equal error probability distribution contrasted against log10(1/16).
     */
    float log10_emission_odds(char read_base, char probe_base, uint32_t pl);
    
    /**
     * Converts state to string representation.
     */
    std::string state2string(int32_t state);

    /**
     * Converts state to cigar string representation.
     */
    std::string state2cigarstring(int32_t state);

    /**
     * Converts track to cigar string representation.
     */
    std::string track2cigarstring1(int32_t t, int32_t j);

    /**
     * Converts state to cigar string representation.
     */
    std::string track2cigarstring2(int32_t t);

    /**
     * Converts model component to string representation.
     */
    std::string component2string(int32_t component);

    /**
     * Prints an alignment.
     */
    void print_alignment();

    /**
     * Prints an alignment with padding.
     */
    void print_alignment(std::string& pad);

    /**
     * Prints a float matrix.
     */
    void print(float *v, size_t plen, size_t rlen);

    /**
     * Prints a int32_t matrix.
     */
    void print(int32_t *v, size_t plen, size_t rlen);

    /**
     * Prints the transition matrix.
     */
    void print_T();

    /**
     * Prints U.
     */
    void print_U(int32_t *U, size_t plen, size_t rlen);

    /**
     * Prints U and V.
     */
    void print_trace(int32_t state, size_t plen, size_t rlen);

    /**
     * Collect alignment summary statistics.
     */
    void collect_statistics(int32_t t1, int32_t t2, int32_t j);

    /**
     * Clear alignment statistics.
     */
    void clear_statistics();

    /**
     * Update alignment statistics after collection.
     */
    void update_statistics();

    /**
     * Prints track.
     */
    void print_track(int32_t t);
};

#endif