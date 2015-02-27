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

#ifndef BAM_ORDERED_READER_H
#define BAM_ORDERED_READER_H

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include "genome_interval.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "hts_utils.h"
#include "interval_tree.h"
#include <list>
#include <map>
#include <queue>
#include <vector>

/**
 * A class for reading ordered BAM files.
 * 
 * Basically a record iterator that hides the
 * htslib interface from the programs in vt.
 *
 * The main impetus for this is that htslib
 * is currently incorporated as a very early
 * stage and is thus lacking many feature that
 * is useful to us at this point in time. This
 * allows us to isolate the changes required in 
 * future to simply methods in hts_utils.
 *
 * The following cases are supported.
 *
 * 1) Input is an unindexed file which is not necessarily ordered.
 * 2) Input is an indexed file
 *
 * This class hides the handling of indices from 
 * the user and also allows for the selection of 
 * records in intervals in both cases 1 and 2.
 */
 
class BAMOrderedReader
{
    public:
        
    ///////
    //i/o//
    ///////
    std::string file_name;    
    samFile *file;
    bam_hdr_t *hdr;
    hts_idx_t *idx;
    hts_itr_t *itr; 
    bam1_t *s;
    
    //for control
    htsFormat ftype;
    bool intervals_present;
    bool index_loaded;
    bool random_access_enabled;

    //common use
    kstring_t str;
        
    //list of intervals
    std::vector<GenomeInterval> intervals; 
    uint32_t interval_index;    
             
    /**
     * Initialize files, intervals and reference file. 
     *
     * @file_name        name of the input VCF file
     * @intervals        list of intervals, if empty, all records are selected.
     * @ref_fasta_file   reference FASTA file for CRAM
     */
    BAMOrderedReader(std::string file_name, std::vector<GenomeInterval>& intervals, std::string ref_fasta_file="");
    
    /**
     * Jump to interval. Returns false if not successful.
     *
     * @interval - string representation of interval.
     */
    bool jump_to_interval(GenomeInterval& interval);
      
    /**
     * Returns next vcf record.
     */
    bool read(bam1_t *s);
    
    /**
     * Closes the file.
     */
    void close();
    
    private:
    /**
     * Initialize next interval.
     * Returns false only if all intervals are accessed.
     */
    bool initialize_next_interval();
};
    
#endif