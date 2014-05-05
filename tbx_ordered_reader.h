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

#ifndef TBX_ORDERED_READER_H
#define TBX_ORDERED_READER_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "hts_utils.h"
#include "genome_interval.h"

/**
 * A class for reading ordered tabixed files.
 *
 * Basically a record iterator that hides the
 * htslib interface from the programs in vt.
 *
 * The main impetus for this is that htslib
 * is currently incorporated as a very early
 * stage and is thus lacking many feature that
 * is useful to us at this point in time. This
 * allows us to isolate the changes required in
 * future to simply methods in hts_utils.
 *
 * The following cases are supported.
 *
 * 1) Input is an unindexed file which is not necessarily ordered.
 * 2) Input is an indexed file
 *
 * This class hides the handling of indices from
 * the user and also allows for the selection of
 * records in intervals in both cases 1 and 2.
 */

class TBXOrderedReader
{
    public:

    ///////
    //i/o//
    ///////
    std::string hts_file;
    htsFile *hts;
    tbx_t *tbx;
    hts_itr_t *itr;

    //for control
    bool intervals_present;
    bool index_loaded;
    bool random_access_enabled;

    //list of intervals
    std::vector<GenomeInterval> intervals;
    uint32_t interval_index;

    //shared objects for string manipulation
    kstring_t s;
    
    /**
     * Initialize files and intervals.
     *
     * @hts              name of the input file
     */
    TBXOrderedReader(std::string& hts_file);
   
    /**
     * Initialize files and intervals.
     *
     * @hts              name of the input file
     * @intervals        list of intervals, if empty, all records are selected.
     */
    TBXOrderedReader(std::string& hts_file, std::vector<GenomeInterval>& intervals);

    /**
     * Jump to interval. Returns false if not successful.
     *
     * @interval - string representation of interval.
     */
    bool jump_to_interval(GenomeInterval& interval);

    /**
     * Returns next vcf record.
     */
    bool read(kstring_t *s);

    /**
    * Initialize next interval.
    * Returns false only if all intervals are accessed.
    */
    bool initialize_next_interval();

    /**
     * Returns next set of vcf records at a start position.
     * Note that this function should never be used in conjunction with read(bcf1_t *v)
     */
    bool read_next_position(std::vector<bcf1_t *>& vs);

    /**
     * Closes the file.
     */
    void close();

    private:
};

#endif
