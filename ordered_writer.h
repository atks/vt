/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>

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

#ifndef ORDERED_WRITER_H
#define ORDERED_WRITER_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <sstream>
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/tbx.h"
#include "hts_utils.h"

/**
 * A class for writing ordered VCF/BCF files.
 * 
 * In some cases, the processing of a file changes the coordinates 
 * slightly leading to slight local disorder within a contig,
 * instead of sorting the VCF wholesale, this class buffers the output
 * and sorts locally in a 10K base pair region before writing the records
 * out.
 */
class OrderedWriter
{
    public:
        
    ///////
    //i/o//
    /////// 
    std::string vcf_file;    
    vcfFile *vcf;
    bcf_hdr_t *hdr;  
        
    //buffer for containing records to be written out
    std::list<bcf1_t*> buffer; //most recent records in the front
    std::list<bcf1_t*> pool;  
    
    //shared objects for string manipulation
    kstring_t s;
    std::stringstream ss; 
    
    int32_t window;
    
    /**
     * Initialize output file.
     */
    OrderedWriter(std::string _input_vcf_file, int32_t _window=100000);
          
    /**
     * Gets record from pool, creates a new record if necessary
     */
    void set_hdr(bcf_hdr_t *_hdr);
    
    /**
     * Reads next record, hides the random access of different regions from the user.
     */
    void write_hdr();
    
    /**
     * Reads next record, hides the random access of different regions from the user.
     */
    void write1(bcf1_t *v);
    
    /**
     * Gets sequence name of a record
     */
    const char* get_seqname(bcf1_t *v);

    /**
     * Gets record from pool, creates a new record if necessary
     */
    bcf1_t* get_bcf1_from_pool();
        
    /**
     * Flush writable records from buffer.
     */
    void flush();
    
    private:
    /**
     * Returns record to pool 
     */ 
    void store_bcf1_into_pool(bcf1_t* v);
    
    /**
     * Flush writable records from buffer.
     */
    void flush(bool force);
};
    
#endif