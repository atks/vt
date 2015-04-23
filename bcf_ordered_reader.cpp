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

#include "bcf_ordered_reader.h"

BCFOrderedReader::BCFOrderedReader(std::string file_name, std::vector<GenomeInterval>& intervals)
{
    this->file_name = (file_name=="+")? "-" : file_name;
    file = NULL;
    hdr = NULL;
    idx = NULL;
    tbx = NULL;
    itr = NULL;

    this->intervals = intervals;
    interval_index = 0;
    index_loaded = false;

    file = hts_open(this->file_name.c_str(), "r");
    if (!file)
    {
        fprintf(stderr, "[%s:%d %s] Cannot open %s\n", __FILE__, __LINE__, __FUNCTION__, file_name.c_str());
        exit(1);    
    }    
    ftype = file->format;

    if (ftype.format!=vcf && ftype.format!=bcf)
    {
        fprintf(stderr, "[%s:%d %s] Not a VCF/BCF file: %s\n", __FILE__, __LINE__, __FUNCTION__, file_name.c_str());
        exit(1);
    }

    s = {0, 0, 0};
    if (file==NULL) exit(1);
    hdr = bcf_alt_hdr_read(file);
    if (!hdr) exit(1);

    intervals_present =  intervals.size()!=0;

    if (ftype.format==bcf)
    {   
        if ((idx = bcf_index_load(file_name.c_str())))
        {
            index_loaded = true;
        }
        else
        {    
            if (intervals_present)
            {
                fprintf(stderr, "[E:%s] index cannot be loaded for %s for random access\n", __FUNCTION__, file_name.c_str());
                exit(1);
            }
        }
    }
    else if (ftype.format==vcf)
    {
        if (ftype.compression==bgzf)
        {    
            if ((tbx = tbx_index_load(file_name.c_str())))
            {
                index_loaded = true;
            }
            else
            {
                if (intervals_present)
                {
                    fprintf(stderr, "[E:%s] index cannot be loaded for %s for random access\n", __FUNCTION__, file_name.c_str());
                    exit(1);
                }
            }
        }
        else
        {
            if (intervals_present)
            {
                fprintf(stderr, "[E:%s] no random access support for VCF file: %s\n", __FUNCTION__, file_name.c_str());
                exit(1);
            }
        }
    }

    random_access_enabled = intervals_present && index_loaded;
};

/**
 * Jump to interval. Returns false if not successful.
 *
 * @interval - string representation of interval.
 */
bool BCFOrderedReader::jump_to_interval(GenomeInterval& interval)
{
    if (index_loaded)
    {
        intervals_present = true;
        random_access_enabled = true;
        intervals.clear();
        intervals.push_back(interval);
        interval_index = 0;
        if (ftype.format==bcf)
        {
            intervals[interval_index++].to_string(&s);
            itr = bcf_itr_querys(idx, hdr, s.s);
            if (itr)
            {
                return true;
            }
        }
        else if (ftype.format==vcf && ftype.compression==bgzf)
        {
            intervals[interval_index++].to_string(&s);
            itr = tbx_itr_querys(tbx, s.s);
            if (itr)
            {
                return true;
            }
        }
    }

    return false;
};

/**
 * Gets sequence name of a record.
 */
const char* BCFOrderedReader::get_seqname(bcf1_t *v)
{
    return bcf_get_chrom(hdr, v);
};

/**
 * Checks if index is loaded.
 */
bool BCFOrderedReader::is_index_loaded()
{
    return index_loaded;
}

/**
 * Gets bcf header.
 */
bcf_hdr_t* BCFOrderedReader::get_hdr()
{
    return hdr;
};

/**
 * Initialize next interval.
 * Returns false only if all intervals are accessed.
 */
bool BCFOrderedReader::initialize_next_interval()
{
    while (interval_index!=intervals.size())
    {
        if (ftype.format==bcf)
        {
            intervals[interval_index++].to_string(&s);
            itr = bcf_itr_querys(idx, hdr, s.s);
            if (itr)
            {
                return true;
            }
        }
        else if (ftype.format==vcf && ftype.compression==bgzf)
        {
            intervals[interval_index++].to_string(&s);
            itr = tbx_itr_querys(tbx, s.s);
            if (itr)
            {
                return true;
            }
        }
    }

    return false;
};

/**
 * Reads next record, hides the random access of different regions from the user.
 */
bool BCFOrderedReader::read(bcf1_t *v)
{
    if (random_access_enabled)
    {
        if (ftype.format==bcf)
        {
            while(true)
            {
                if (itr && bcf_itr_next(file, itr, v)>=0)
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
            while(true)
            {
                if (itr && tbx_itr_next(file, tbx, itr, &s)>=0)
                {
                    vcf_parse1(&s, hdr, v);
                    return true;
                }
                else if (!initialize_next_interval())
                {
                    return false;
                }
            }
        }
    }
    else
    {
        if (bcf_read(file, hdr, v)==0)
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
 * Returns record to pool
 */
void BCFOrderedReader::store_bcf1_into_pool(bcf1_t* v)
{
    bcf_clear(v);
    pool.push_back(v);
}

/**
 * Gets record from pool, creates a new record if necessary
 */
bcf1_t* BCFOrderedReader::get_bcf1_from_pool()
{
    if(!pool.empty())
    {
        bcf1_t* v = pool.front();
        pool.pop_front();
        return v;
    }
    else
    {
        return bcf_init1();
    }
};

/**
 * Closes the file.
 */
void BCFOrderedReader::close()
{
    bcf_close(file);
    if (hdr) bcf_hdr_destroy(hdr);
    hdr = NULL;
}
