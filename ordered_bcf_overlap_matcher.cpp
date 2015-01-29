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

    if (current_interval.seq!=chrom)
    {
        std::list<bcf1_t*>::iterator i = buffer.begin();
        while (i!=buffer.end())
        {
            bcf_destroy(*i);
            i = buffer.erase(i);
        }
        
        current_interval.set(chrom);
        if (!odr->jump_to_interval(current_interval))
        {
            fprintf(stderr, "[E:%s] cannot jump to %s\n", __FUNCTION__, current_interval.to_string().c_str());
            exit(1);
        }    
        std::cerr << "Jumped to chromosome " << chrom << "\n";
        
        v = bcf_init();
        
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            if (bcf_get_end_pos1(v)<start1) continue;
            overlaps = overlaps || (bcf_get_pos1(v)<=end1);
            buffer.push_back(v);
            if (bcf_get_pos1(v)>end1) break;
                
            v = bcf_init();
        }
    }
    else
    {
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
        
        v = bcf_init();
        
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            if (bcf_get_end_pos1(v)<start1) continue;
            overlaps = overlaps || (bcf_get_pos1(v)<=end1);
            buffer.push_back(v);
            if (bcf_get_pos1(v)>end1) break;
            
            v = bcf_init();
        }
    }

    return overlaps;
};

/**
 * Returns true if chrom:start1-end1 overlaps with a region in the file and populates the overlapping variants.
 */
bool OrderedBCFOverlapMatcher::overlaps_with(std::string& chrom, int32_t start1, int32_t end1, std::vector<bcf1_t*>& overlap_vars)
{
    overlap_vars.clear();
    bool overlaps = false;

    if (current_interval.seq!=chrom)
    {
        //clear records from previous chromosome
        std::list<bcf1_t*>::iterator i = buffer.begin();
        while (i!=buffer.end())
        {
            bcf_destroy(*i);
            i = buffer.erase(i);
        }
        
        //random access next chromosome
        current_interval.set(chrom);
        if (!odr->jump_to_interval(current_interval))
        {
            fprintf(stderr, "[E:%s] cannot jump to %s\n", __FUNCTION__, current_interval.to_string().c_str());
            return false;
        }    
        //std::cerr << "Jumped to chromosome " << chrom << "\n";
        
        //read new variants
        v = bcf_init();            
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            
            if (bcf_get_end_pos1(v)<start1) 
            {
                continue;
            }
            
            if (bcf_get_pos1(v)>end1)
            {
                buffer.push_back(v);
                v = NULL;
                break;
            }
            
            overlaps = true;
            buffer.push_back(v);
            overlap_vars.push_back(v);
            v = bcf_init();
        }
        
        if (v)
        {
            bcf_destroy(v);
        }
    }
    else
    {   
        bool need_to_read = true;
        
        //scythe records that occur prior to chrom:start1-end1
        std::list<bcf1_t*>::iterator i = buffer.begin();
        while (i!=buffer.end())
        {
            if (bcf_get_end_pos1(*i)<start1)
            {
                bcf_destroy(*i);
                i = buffer.erase(i);
                continue;
            }

            if (bcf_get_pos1(*i)>end1)
            {
                need_to_read = false;
                break;
            }
            
            overlaps = true;
            overlap_vars.push_back(*i);
            
            ++i;
        }
        
        //read new variants
        if (need_to_read)
        {
            v = bcf_init();
            
            while (odr->read(v))
            {
                bcf_unpack(v, BCF_UN_STR);
                
                if (bcf_get_end_pos1(v)<start1) 
                {
                    continue;
                }
                
                if (bcf_get_pos1(v)>end1)
                {
                    buffer.push_back(v);
                    v = NULL;
                    break;
                }
                
                overlaps = true;
                buffer.push_back(v);
                overlap_vars.push_back(v);
                v = bcf_init();
            }
            
            if (v)
            {
                bcf_destroy(v);
            }
        }
    }

    return overlaps;
};
