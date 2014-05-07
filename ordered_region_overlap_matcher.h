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

#ifndef OVERLAP_REGION_MATCHER_H
#define OVERLAP_REGION_MATCHER_H

#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "program.h"
#include "hts_utils.h"
#include "tbx_ordered_reader.h"
#include "interval_tree.h"
#include <list>
#include "bed.h"

class OrderedRegionOverlapMatcher
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_file;   
    
    ///////
    //i/o//
    ///////
    TBXOrderedReader *todr;    
    
    kstring_t s;
    
    GenomeInterval current_interval;
    std::list<BEDRecord> buffer;
    bool end_of_file;
    int32_t no_regions;
    
    /**
     * Constructor.
     */
    OrderedRegionOverlapMatcher(std::string& file);

    /**
     * Destructor.
     */
    ~OrderedRegionOverlapMatcher();
    
    /**
     * Returns true if chrom:start1-end1 overlaps with a region in the file.
     */
    bool overlaps_with(std::string& chrom, int32_t start1, int32_t end1);
        
    private:
};
    
#endif
