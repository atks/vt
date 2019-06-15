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
    bcf_hdr_append_info_with_backup_naming(odr->hdr, "EXACT_OVERLAPS", "1", "Integer", "Number of exact overlapping variants with this variant.", true);
    bcf_hdr_append_info_with_backup_naming(odr->hdr, "FUZZY_OVERLAPS", "1", "Integer", "Number of fuzzy overlapping variants with this variant.", true);
    if (bcf_hdr_sync(odr->hdr)<0) 
    {
        fprintf(stderr, "[%s:%d %s] Cannot update header\n", __FILE__, __LINE__, __FUNCTION__);
        exit(1);
    }
    
    no_regions = 0;
    current_interval.seq = "";
    no_exact_overlaps = 0;
    no_fuzzy_overlaps = 0;
    no_nonoverlaps = 0;
    no_variants = 0;
};

/**
 * Constructor.
 */
OrderedBCFOverlapMatcher::OrderedBCFOverlapMatcher(std::string& file, std::vector<GenomeInterval>& intervals, std::string fexp)
{
    odr = new BCFOrderedReader(file, intervals);
    bcf_hdr_append_info_with_backup_naming(odr->hdr, "EXACT_OVERLAPS", "1", "Integer", "Number of exact overlapping variants with this variant.", true);
    bcf_hdr_append_info_with_backup_naming(odr->hdr, "FUZZY_OVERLAPS", "1", "Integer", "Number of fuzzy overlapping variants with this variant.", true);
    if (bcf_hdr_sync(odr->hdr)<0) 
    {
        fprintf(stderr, "[%s:%d %s] Cannot update header\n", __FILE__, __LINE__, __FUNCTION__);
        exit(1);
    }

    filter.parse(fexp.c_str());
    filter_exists = true;

    no_regions = 0;
    current_interval.seq = "";
    no_exact_overlaps = 0;
    no_fuzzy_overlaps = 0;
    no_nonoverlaps = 0;
    no_variants = 0;
};

/**
 * Destructor.
 */
OrderedBCFOverlapMatcher::~OrderedBCFOverlapMatcher()
{
    odr->close();
    delete odr;
};

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
            fprintf(stderr, "[%s:%d %s] cannot jump to %s\n", __FILE__, __LINE__, __FUNCTION__, current_interval.to_string().c_str());
            exit(1);
        }

        v = bcf_init();

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);
            
            if (filter_exists)
            {
                variant.classify(odr->hdr, v);
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }
            
            if (bcf_get_end1(v)<start1) continue;
            overlaps = overlaps || (bcf_get_pos1(v)<=end1);
            buffer.push_back(v);
            if (bcf_get_pos1(v)>end1) break;

            v = bcf_init();
            ++no_variants;
        }
    }
    else
    {
        std::list<bcf1_t*>::iterator i = buffer.begin();
        while (i!=buffer.end())
        {
            if (bcf_get_end1(*i)<start1)
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
            bcf_unpack(v, BCF_UN_INFO);
            
            if (filter_exists)
            {
                variant.classify(odr->hdr, v);
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }
            
            ++no_variants;
            if (bcf_get_end1(v)<start1) continue;
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
 * This ensures that all records in the reference VCF is processed to compute accurate overlap statistics.
 */
bool OrderedBCFOverlapMatcher::overlaps_with(int32_t rid, int32_t beg1, int32_t end1, std::vector<bcf1_t*>& overlap_vars)
{
    overlap_vars.clear();
    bool overlaps = false;

    bool need_to_read = true;

    //scythe records that occur prior to chrom:start1-end1
    std::list<bcf1_t*>::iterator i = buffer.begin();
    while (i!=buffer.end())
    {
        int32_t crid =  bcf_get_rid(*i);
        int32_t cbeg1 = bcf_get_pos1(*i);
        int32_t cend1 = bcf_get_end1(*i);

        //drop variant
        if (crid<rid || (crid==rid && cend1<beg1))
        {
            update_overlap_statistics(*i);
            bcf_destroy(*i);
            i = buffer.erase(i);
            continue;
        }
        
        //buffer is ahead of current record, no need to read new records
        if ((crid==rid && cbeg1>end1) || crid>rid)
        {
            need_to_read = false;
            break;
        }

        if (beg1==cbeg1 && end1==cend1)
        {
            increment_exact_overlap(*i);
        }
        else
        {
            increment_fuzzy_overlap(*i);
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
            bcf_unpack(v, BCF_UN_INFO);
            
            if (filter_exists)
            {
                variant.classify(odr->hdr, v);
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }
            
            ++no_variants;

            int32_t crid =  bcf_get_rid(v);
            int32_t cbeg1 = bcf_get_pos1(v);
            int32_t cend1 = bcf_get_end1(v);

            //drop variant
            if (crid<rid || (crid==rid && cend1<beg1))
            {
                update_overlap_statistics(v);
                continue;
            }

            //add and stop reading from file
            if ((crid==rid && cbeg1>end1) || crid>rid)
            {
                buffer.push_back(v);
                v = NULL;
                break;
            }

            //exact overlap
            if (beg1==cbeg1 && end1==cend1)
            {
                increment_exact_overlap(v);
            }
            else
            {
                increment_fuzzy_overlap(v);
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

    return overlaps;
};

/**
 * Returns true if chrom:start1-end1 overlaps with a region in the file and populates the overlapping variants.
 * This ensures that all records in the reference VCF is processed to compute accurate overlap statistics.
 * Flushed variants are written to odw.
 */
bool OrderedBCFOverlapMatcher::overlaps_with(int32_t rid, int32_t beg1, int32_t end1, std::vector<bcf1_t*>& overlap_vars, BCFOrderedWriter* odw)
{
    overlap_vars.clear();
    bool overlaps = false;

    bool need_to_read = true;

    //scythe records that occur prior to chrom:start1-end1
    std::list<bcf1_t*>::iterator i = buffer.begin();
    while (i!=buffer.end())
    {
        int32_t crid =  bcf_get_rid(*i);
        int32_t cbeg1 = bcf_get_pos1(*i);
        int32_t cend1 = bcf_get_end1(*i);

        //drop variant
        if (crid<rid || (crid==rid && cend1<beg1))
        {
            update_overlap_statistics(*i, odw);
            bcf_destroy(*i);
            i = buffer.erase(i);
            continue;
        }
        
        //buffer is ahead of current record, no need to read new records
        if ((crid==rid && cbeg1>end1) || crid>rid)
        {
            need_to_read = false;
            break;
        }

        if (beg1==cbeg1 && end1==cend1)
        {
            increment_exact_overlap(*i);
        }
        else
        {
            increment_fuzzy_overlap(*i);
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
            bcf_unpack(v, BCF_UN_INFO);
                        
            if (filter_exists)
            {
                variant.classify(odr->hdr, v);
                if (!filter.apply(odr->hdr, v, &variant))
                {
                    continue;
                }
            }
            ++no_variants;

            int32_t crid =  bcf_get_rid(v);
            int32_t cbeg1 = bcf_get_pos1(v);
            int32_t cend1 = bcf_get_end1(v);

            //drop variant
            if (crid<rid || (crid==rid && cend1<beg1))
            {
                update_overlap_statistics(v, odw);
                continue;
            }

            //add and stop reading from file
            if ((crid==rid && cbeg1>end1) || crid>rid)
            {
                buffer.push_back(v);
                v = NULL;
                break;
            }

            //exact overlap
            if (beg1==cbeg1 && end1==cend1)
            {
                increment_exact_overlap(v);
            }
            else
            {
                increment_fuzzy_overlap(v);
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

    return overlaps;
};

/**
 * Flushes remaining variants.
 */
void OrderedBCFOverlapMatcher::flush()
{
    //clear records from previous chromosome
    std::list<bcf1_t*>::iterator i = buffer.begin();
    while (i!=buffer.end())
    {
        update_overlap_statistics(*i);
        bcf_destroy(*i);
        i = buffer.erase(i);
    }

    v = bcf_init();

    while (odr->read(v))
    {        
        if (filter_exists)
        {
            variant.classify(odr->hdr, v);
            if (!filter.apply(odr->hdr, v, &variant))
            {
                continue;
            }
        }
        ++no_variants;
        update_overlap_statistics(v);
    }

    if (v)
    {
        bcf_destroy(v);
    }
}

/**
 * Flushes remaining variants.
 * Flushed variants are written to odw.
 */
void OrderedBCFOverlapMatcher::flush(BCFOrderedWriter* odw)
{
    //clear records from previous chromosome
    std::list<bcf1_t*>::iterator i = buffer.begin();
    while (i!=buffer.end())
    {
        update_overlap_statistics(*i, odw);
        bcf_destroy(*i);
        i = buffer.erase(i);
    }

    v = bcf_init();

    while (odr->read(v))
    {        
        if (filter_exists)
        {
            variant.classify(odr->hdr, v);
            if (!filter.apply(odr->hdr, v, &variant))
            {
                continue;
            }
        }
        
        ++no_variants;
        update_overlap_statistics(v, odw);
    }

    if (v)
    {
        bcf_destroy(v);
    }
}

/**
 * Closes the file.
 */
void OrderedBCFOverlapMatcher::close()
{
    odr->close();
}

/**
 * Increments the EXACT_OVERLAPS count of a variant record.
 */
void OrderedBCFOverlapMatcher::increment_exact_overlap(bcf1_t* v)
{
    int32_t n = 0;
    int32_t *count = NULL;
    bcf_unpack(v, BCF_UN_INFO);
    if (bcf_get_info_int32(odr->hdr, v, "EXACT_OVERLAPS", &count, &n)>0)
    {
//        std::cerr << "exact " << count[0] << " "  << n << "\n";
        count[0] = std::min(127, count[0]+1);
        bcf_update_info_int32(odr->hdr, v, "EXACT_OVERLAPS", count, n);
        free(count);
    }
    else
    {
        int32_t c = 1;
        bcf_update_info_int32(odr->hdr, v, "EXACT_OVERLAPS", &c, 1);
    }
}

/**
 * Increments the FUZZY_OVERLAPS count of a variant record.
 */
void OrderedBCFOverlapMatcher::increment_fuzzy_overlap(bcf1_t* v)
{
    int32_t n = 0;
    int32_t *count = NULL;
    bcf_unpack(v, BCF_UN_INFO);
    if (bcf_get_info_int32(odr->hdr, v, "FUZZY_OVERLAPS", &count, &n)>0)
    {
        count[0] = std::min(127, count[0]+1);
        bcf_update_info_int32(odr->hdr, v, "FUZZY_OVERLAPS", count, n);
        free(count);
    }
    else
    {
        int32_t c = 1;
        bcf_update_info_int32(odr->hdr, v, "FUZZY_OVERLAPS", &c, 1);
    }
}

/**
 * Updates the number of non overlapping and non overlapping variants.
 * This is always invoked when a variant is flushed.
 */
void OrderedBCFOverlapMatcher::update_overlap_statistics(bcf1_t* v)
{
    int32_t exact_n = 0;
    int32_t fuzzy_n = 0;
    int32_t *exact_count = NULL;
    int32_t *fuzzy_count = NULL;
    int32_t ret_exact_overlaps = bcf_get_info_int32(odr->hdr, v, "EXACT_OVERLAPS", &exact_count, &exact_n);
    int32_t ret_fuzzy_overlaps = bcf_get_info_int32(odr->hdr, v, "FUZZY_OVERLAPS", &fuzzy_count, &fuzzy_n);
    if (ret_exact_overlaps==-3 && ret_fuzzy_overlaps==-3)
    {
        ++no_nonoverlaps;
    }
    else if (ret_exact_overlaps==-3 && ret_fuzzy_overlaps>0)
    {
        ++no_fuzzy_overlaps;
    }
    //priority to exact matches
    else if (ret_exact_overlaps>0)
    {
        ++no_exact_overlaps;
    }
    else
    {
        std::cerr << "UNACCOUNTED\n";
        bcf_print(odr->hdr,v);
    }

    if (exact_count) free(exact_count);
    if (fuzzy_count) free(fuzzy_count);
}

/**
 * Updates the number of non overlapping and non overlapping variants.
 * This is always invoked when a variant is flushed.
 * Record is written to odw.
 */
void OrderedBCFOverlapMatcher::update_overlap_statistics(bcf1_t* v, BCFOrderedWriter* odw)
{
    int32_t exact_n = 0;
    int32_t fuzzy_n = 0;
    int32_t *exact_count = NULL;
    int32_t *fuzzy_count = NULL;
    int32_t ret_exact_overlaps = bcf_get_info_int32(odr->hdr, v, "EXACT_OVERLAPS", &exact_count, &exact_n);
    int32_t ret_fuzzy_overlaps = bcf_get_info_int32(odr->hdr, v, "FUZZY_OVERLAPS", &fuzzy_count, &fuzzy_n);
    if (ret_exact_overlaps==-3 && ret_fuzzy_overlaps==-3)
    {
        ++no_nonoverlaps;
    }
    else if (ret_exact_overlaps==-3 && ret_fuzzy_overlaps>0)
    {
        ++no_fuzzy_overlaps;
    }
    //priority to exact matches
    else if (ret_exact_overlaps>0)
    {
        ++no_exact_overlaps;
    }
    else
    {
        std::cerr << "UNACCOUNTED\n";
        bcf_print(odr->hdr,v);
    }
   
    if (exact_count) free(exact_count);
    if (fuzzy_count) free(fuzzy_count);

    if (odw) odw->write(v);
}

/**
 * Get number of exact overlap variants that has been printed and reset no_exact_overlaps.
 */
int32_t OrderedBCFOverlapMatcher::get_no_exact_overlaps()
{
    int32_t val = no_exact_overlaps;
    no_exact_overlaps = 0;

    return val;
}

/**
 * Get number of fuzzy overlap variants that has been printed and reset no_fuzzy_overlaps.
 */
int32_t OrderedBCFOverlapMatcher::get_no_fuzzy_overlaps()
{
    int32_t val = no_fuzzy_overlaps;
    no_fuzzy_overlaps = 0;

    return val;
}

/**
 * Get number of non-overlapping variants that has been printed and reset no_nonoverlaps.
 */
int32_t OrderedBCFOverlapMatcher::get_no_nonoverlaps()
{
    int32_t val = no_nonoverlaps;
    no_nonoverlaps = 0;

    return val;
}

/**
 * Is this record and exact match?.
 */
bool OrderedBCFOverlapMatcher::is_exact_match(int32_t rid, int32_t beg1, int32_t end1, bcf1_t* v)
{
    return (rid==bcf_get_rid(v) && beg1==bcf_get_pos1(v) && end1==bcf_get_end1(v));        
}
