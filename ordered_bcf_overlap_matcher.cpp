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

    bcf_hdr_sync(odr->hdr);
    no_regions = 0;
    current_interval.seq = "";
    no_exact_overlaps = 0;
    no_fuzzy_overlaps = 0;
    no_nonoverlaps = 0;
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
            fprintf(stderr, "[E:%s] cannot jump to %s\n", __FUNCTION__, current_interval.to_string().c_str());
            exit(1);
        }
//        std::cerr << "Jumped to chromosome " << chrom << "\n";

        v = bcf_init();

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);
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
            bcf_unpack(v, BCF_UN_INFO);
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
bool OrderedBCFOverlapMatcher::overlaps_with(std::string& chrom, int32_t beg1, int32_t end1, std::vector<bcf1_t*>& overlap_vars)
{
    overlap_vars.clear();
    bool overlaps = false;

    if (current_interval.seq!=chrom)
    {
        //clear records from previous chromosome
        std::list<bcf1_t*>::iterator i = buffer.begin();
        while (i!=buffer.end())
        {
            update_overlap_statistics(*i);
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
            bcf_unpack(v, BCF_UN_INFO);
            int32_t cbeg1 = bcf_get_pos1(v);
            int32_t cend1 = bcf_get_end_pos1(v);

            if (cend1<beg1)
            {
                continue;
            }

            if (cbeg1>end1)
            {
                buffer.push_back(v);
                v = NULL;
                break;
            }

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
    else
    {
        bool need_to_read = true;

        //scythe records that occur prior to chrom:start1-end1
        std::list<bcf1_t*>::iterator i = buffer.begin();
        while (i!=buffer.end())
        {
            int32_t cbeg1 = bcf_get_pos1(*i);
            int32_t cend1 = bcf_get_end_pos1(*i);

            if (cend1<beg1)
            {
                update_overlap_statistics(*i);
                bcf_destroy(*i);
                i = buffer.erase(i);
                continue;
            }

            if (cbeg1>end1)
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

                int32_t cbeg1 = bcf_get_pos1(v);
                int32_t cend1 = bcf_get_end_pos1(v);

                if (cend1<beg1)
                {
                    continue;
                }

                if (cbeg1>end1)
                {
                    buffer.push_back(v);
                    v = NULL;
                    break;
                }

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
//        ++no_nonoverlaps;
    }

    if (v)
    {
        bcf_destroy(v);
    }
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
//        std::cerr << "fuzzy " << count[0] << " "  << n << "\n";
        count[0] = std::min(127, count[0]+1);
//        std::cerr << "\tto " << count[0] << "\n";
//        bcf_print(odr->hdr, v);
        bcf_update_info_int32(odr->hdr, v, "FUZZY_OVERLAPS", count, n);
        free(count);
    }
    else
    {
        int32_t c = 1;
//        bcf_print(odr->hdr, v);
        bcf_update_info_int32(odr->hdr, v, "FUZZY_OVERLAPS", &c, 1);
    }
}

/**
 * Updates the number of non overlapping and non overlapping variants.
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
        bcf_print(odr->hdr, v);
        ++no_nonoverlaps;
    }
    if (ret_exact_overlaps==-3 && ret_fuzzy_overlaps>0)
    {
        ++no_fuzzy_overlaps;
    }
    if (ret_exact_overlaps>0)
    {
        ++no_exact_overlaps;
    }

    if (exact_count) free(exact_count);
    if (fuzzy_count) free(fuzzy_count);
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
bool OrderedBCFOverlapMatcher::is_exact_match(bcf1_t* v)
{
    int32_t exact_n = 0;
    int32_t *exact_count = NULL;
    int32_t ret_exact_overlaps = bcf_get_info_int32(odr->hdr, v, "EXACT_OVERLAPS", &exact_count, &exact_n);
    if (exact_count) free(exact_count);
    return (ret_exact_overlaps>0);
}
