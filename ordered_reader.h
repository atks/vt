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

#ifndef ORDERED_READER_H
#define ORDERED_READER_H

#include <queue>
#include <sstream>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <map>
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "hts_utils.h"

/**
 * A class for reading ordered VCF/BCF files.
 * Hides the handling of indices from the user. 
 *
 * Supported for the following cases:
 *   
 * 1) Input is streamed, order is examined.
 * 
 * 2) Input is an indexed file, with selected regions.
 */
class OrderedReader
{
    public:
        
    ///////
    //i/o//
    /////// 
    std::string vcf_file;    
    vcfFile *vcf;
    BGZF *vcfgz; 
    bcf_hdr_t *hdr;  
    hts_idx_t *idx; 
    tbx_t *tbx; 
    hts_itr_t *itr; 
    bcf1_t *v;
    
    //for control
    int32_t vcf_ftype;
    bool need_random_access;
        
    //list of intervals
    std::vector<std::string> intervals; //contains intervals of interest, if build is 
    int32_t interval_index;    
        
    //shared objects for string manipulation
    kstring_t s;
    std::stringstream ss;
    
    /**
     * Initialize files and intervals. 
     */
    OrderedReader(std::string _input_vcf_file, std::vector<std::string>& _intervals);
      
    /**
     * Returns next vcf record.
     */
    bool read1(bcf1_t *v);
    
    /**
     * Gets sequence name of a record
     */
    const char* get_seqname(bcf1_t *v);
    
    private:
    /**
     * Initialize next interval.
     * Returns false only if all intervals are accessed.
     */
    bool initialize_next_interval();
};
    
#endif