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

#ifndef ORDERED_BCF_OVERLAP_MATCHER_H
#define ORDERED_BCF_OVERLAP_MATCHER_H

#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "filter.h"
#include "hts_utils.h"
#include "utils.h"

class OrderedBCFOverlapMatcher
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_file;   
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;    
    bcf1_t *v;
    
    GenomeInterval current_interval;
    std::list<bcf1_t*> buffer;
    bool end_of_file;
    int32_t no_regions;
   
    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    Variant variant;
    
	///////
    //stats
    ///////
    int32_t no_exact_overlaps;
    int32_t no_fuzzy_overlaps;
    int32_t no_nonoverlaps;
    int32_t no_variants;

    /**
     * Constructor.
     */
    OrderedBCFOverlapMatcher(std::string& file, std::vector<GenomeInterval>& intervals);

    /**
     * Constructor.
     */
    OrderedBCFOverlapMatcher(std::string& file, std::vector<GenomeInterval>& intervals, std::string fexp);
        
    /**
     * Destructor.
     */
    ~OrderedBCFOverlapMatcher();
    
    /**
     * Returns true if chrom:start1-end1 overlaps with a region in the file.
     */
    bool overlaps_with(std::string& chrom, int32_t start1, int32_t end1);
        
    /**
     * Returns true if chrom:start1-end1 overlaps with a region in the file and populates the overlapping variants.
     */
    bool overlaps_with(int32_t rid, int32_t start1, int32_t end1, std::vector<bcf1_t*>& overlap_vars);
    
    /**
     * Returns true if chrom:start1-end1 overlaps with a region in the file and populates the overlapping variants.
     * This ensures that all records in the reference VCF is processed to compute accurate overlap statistics.
     * Flushed variants are written to odw.
     */
    bool overlaps_with(int32_t rid, int32_t beg1, int32_t end1, std::vector<bcf1_t*>& overlap_vars, BCFOrderedWriter* odw);
            
    /**
     * Flushes remaining variants.
     */
    void flush();
    
    /**
     * Flushes remaining variants.
     * Flushed variants are written to odw.
     */
    void flush(BCFOrderedWriter* odw);
                
    /**
     * Increments the EXACT_OVERLAPS count of a variant record.
     */
    void increment_exact_overlap(bcf1_t* v);
    
    /**
     * Increments the FUZZY_OVERLAPS count of a variant record.
     */
    void increment_fuzzy_overlap(bcf1_t* v);
    
    /**
     * Updates the number of non overlapping and non overlapping variants.
     */
    void update_overlap_statistics(bcf1_t* v);
    
    /**
     * Updates the number of non overlapping and non overlapping variants.
     * This is always invoked when a variant is flushed.
     * Record is written to odw.
     */
    void update_overlap_statistics(bcf1_t* v, BCFOrderedWriter* odw);
    
    /**
     * Get number of exact overlap variants that has been printed and reset no_exact_overlaps.
     */
    int32_t get_no_exact_overlaps();
    
    /**
     * Get number of fuzzy overlap variants that has been printed and reset no_fuzzy_overlaps.
     */
    int32_t get_no_fuzzy_overlaps();

    /**
     * Get number of non-overlapping variants that has been printed and reset no_nonoverlaps.
     */
    int32_t get_no_nonoverlaps();

    /**
     * Is this record and exact match?.
     */
    bool is_exact_match(int32_t rid, int32_t beg1, int32_t end1, bcf1_t* v);
    
    private:
};
    
#endif
