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

#include "bam_ordered_reader.h"

BAMOrderedReader::BAMOrderedReader(std::string bam_file, std::vector<GenomeInterval>& intervals)
:bam_file(bam_file), intervals(intervals), sam(0), hdr(0), idx(0), itr(0)
{
    const char* fname = bam_file.c_str();
    int len = strlen(fname);
    if ( !strcasecmp(".bam", fname+len-4) && !strcasecmp("-", fname) && !strcasecmp(".cram", fname+len-5))
    {
        fprintf(stderr, "[%s:%d %s] Not a BAM/CRAM file: %s\n", __FILE__, __LINE__, __FUNCTION__, bam_file.c_str());
        exit(1);
    }
        
    sam = sam_open(bam_file.c_str(), "r");
    hdr = sam_hdr_read(sam);
    s = bam_init1();

    idx = bam_index_load(bam_file.c_str());
    if (idx==0)
    {
        //fprintf(stderr, "[%s:%d %s] fail to load index for %s\n", __FILE__, __LINE__, __FUNCTION__, bam_file.c_str());
        index_loaded = false;
    }
    else
    {
        index_loaded = true;
    }

    str = {0,0,0};

    intervals_present =  intervals.size()!=0;
    interval_index = 0;

    random_access_enabled = intervals_present && index_loaded;
    
};

/**
 * Jump to interval. Returns false if not successful.
 *
 * @interval - string representation of interval.
 */
bool BAMOrderedReader::jump_to_interval(GenomeInterval& interval)
{
    if (index_loaded)
    {
        intervals_present = true;
        random_access_enabled = true;
        intervals.clear();
        intervals.push_back(interval);
        interval_index = 0;
        
        intervals[interval_index++].to_string(&str);
        itr = bam_itr_querys(idx, hdr, str.s);
        if (itr)
        {
            return true;
        }
    }

    return false;
};

/**
 * Initialize next interval.
 * Returns false only if all intervals are accessed.
 */
bool BAMOrderedReader::initialize_next_interval()
{
    while (interval_index!=intervals.size())
    {
        intervals[interval_index++].to_string(&str);
        itr = bam_itr_querys(idx, hdr, str.s);

        if (itr)
        {
            return true;
        }
    }

    return false;
};

/**
 * Reads next record, hides the random access of different regions from the user.
 */
bool BAMOrderedReader::read(bam1_t *s)
{
    if (random_access_enabled)
    {
        while(true)
        {
            if (itr && bam_itr_next(sam, itr, s)>=0)
            {
                return true;
            }
            else if (!initialize_next_interval())
            {
                return false;
            }
        }
    }
    else
    {
        if (sam_read1(sam, hdr, s)>=0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    return false;
};

/**
 * Closes the file.
 */
void BAMOrderedReader::close()
{
    sam_close(sam);
}
