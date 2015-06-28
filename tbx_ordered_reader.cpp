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

#include "tbx_ordered_reader.h"

/**
 * Initialize files and intervals.
 *
 * @hts_file   name of the input file
 */
TBXOrderedReader::TBXOrderedReader(std::string& hts_file)
{
    this->hts_file = hts_file;
    interval_index = 0;

    hts = NULL;
    tbx = NULL;
    itr = NULL;

    s = {0, 0, 0};

    hts = hts_open(hts_file.c_str(), "r");

    index_loaded = false;
    if ((tbx = tbx_index_load(hts_file.c_str())))
    {
        index_loaded = true;
    }

    random_access_enabled = index_loaded;
};

/**
 * Initialize files and intervals.
 *
 * @hts_file   -  name of the input file
 * @intervals  -  list of intervals, if empty, all records are selected.
 */
TBXOrderedReader::TBXOrderedReader(std::string& hts_file, std::vector<GenomeInterval>& intervals)
{
    this->hts_file = hts_file;
    this->intervals = intervals;
    interval_index = 0;

    hts = NULL;
    tbx = NULL;
    itr = NULL;

    s = {0, 0, 0};

    hts = hts_open(hts_file.c_str(), "r");

    intervals_present =  intervals.size()!=0;

    if ((tbx = tbx_index_load(hts_file.c_str())))
    {
        index_loaded = true;
    }
    else
    {
        if (intervals_present)
        {
            fprintf(stderr, "[E:%s] index cannot be loaded for %s\n", __FUNCTION__, hts_file.c_str());
            exit(1);
        }
    }

    random_access_enabled = intervals_present && index_loaded;
};

/**
 * Jump to interval. Returns false if not successful.
 *
 * @interval - string representation of interval.
 */
bool TBXOrderedReader::jump_to_interval(GenomeInterval& interval)
{
    if (index_loaded)
    {
        intervals_present = true;
        random_access_enabled = true;
        intervals.clear();
        intervals.push_back(interval);
        interval_index = 0;

        intervals[interval_index++].to_string(&s);
        itr = tbx_itr_querys(tbx, s.s);
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
bool TBXOrderedReader::initialize_next_interval()
{
    while (interval_index!=intervals.size())
    {
        intervals[interval_index++].to_string(&s);
        itr = tbx_itr_querys(tbx, s.s);
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
bool TBXOrderedReader::read(kstring_t *s)
{
    if (random_access_enabled)
    {
        while(true)
        {
            if (itr && tbx_itr_next(hts, tbx, itr, s)>=0)
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
        if (hts_getline(hts, '\n', s) >= 0)
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
void TBXOrderedReader::close()
{
    hts_close(hts);
}


