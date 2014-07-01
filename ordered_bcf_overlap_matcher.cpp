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

#include "ordered_bcf_overlap_matcher.h"

/**
 * Constructor.
 */
OrderedBCFOverlapMatcher::OrderedBCFOverlapMatcher(std::string& file, std::vector<GenomeInterval>& intervals)
{
    odr = new BCFOrderedReader(file, intervals);
    v = bcf_init();
    no_regions = 0;
    current_interval.seq = "";
};

/**
 * Destructor.
 */
OrderedBCFOverlapMatcher::~OrderedBCFOverlapMatcher() {};

/**
 * Returns true if chrom:start1-end1 overlaps with a region in the file.
 */
bool OrderedBCFOverlapMatcher::overlaps_with(std::string& chrom, int32_t start1, int32_t end1)
{
    bool overlaps = false;

    //moves to new chromosome
    if (current_interval.seq!=chrom)
    {
        buffer.clear();
        current_interval.set(chrom);
        odr->jump_to_interval(current_interval);
        std::cerr << "Jumped to chromosome " << chrom << "\n";
        while (odr->read(v))
        {
            if (bcf_get_end_pos1(v)<start1) continue;
            overlaps = overlaps || (bcf_get_pos1(v)<=end1);
            buffer.push_back(v);
            if (bcf_get_pos1(v)>end1) break;
            
            v = bcf_init();
        }
    }
    else
    {
        //scythe preceding bed records
        std::list<bcf1_t*>::iterator i = buffer.begin();
        while (i!=buffer.end())
        {
            if (bcf_get_end_pos1(*i)<start1)
            {
                bcf_destroy(*i);
                i = buffer.erase(i);
                continue;
            }

            overlaps = (bcf_get_pos1(*i)<=end1);
            break;
        }
        
        if (!overlaps)        
        {
            if (buffer.empty())
            {    
                while (odr->read(v))
                {
                    if (bcf_get_end_pos1(*i)<start1) continue;
                    overlaps = overlaps || (bcf_get_pos1(*i)<=end1);
                    buffer.push_back(v);
                    if (bcf_get_pos1(v)>end1) break;
                    
                    v = bcf_init();
                }
            }
        }
    }

    return overlaps;
};
