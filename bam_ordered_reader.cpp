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

char *samfaipath(const char *fn_ref)
{
    char *fn_list = 0;
    if (fn_ref == 0) return 0;
    fn_list = (char*) calloc(strlen(fn_ref) + 5, 1);
    strcat(strcpy(fn_list, fn_ref), ".fai");
    if (access(fn_list, R_OK) == -1) { // fn_list is unreadable
        if (access(fn_ref, R_OK) == -1) {
            fprintf(stderr, "[samfaipath] fail to read file %s.\n", fn_ref);
        } else {
            if (fai_build(fn_ref) == -1) {
                fprintf(stderr, "[samfaipath] fail to build FASTA index.\n");
                free(fn_list); fn_list = 0;
            }
        }
    }
    return fn_list;
};

/**
 * Initialize files, intervals and reference file.
 *
 * @file_name        name of the input VCF file
 * @intervals             list of intervals, if empty, all records are selected.
 * @ref_fasta_file  reference FASTA file for CRAM
 */
BAMOrderedReader::BAMOrderedReader(std::string file_name, std::vector<GenomeInterval>& intervals, std::string ref_fasta_file)
{
    this->file_name = (file_name=="+")? "-" : file_name;
    file = NULL;
    hdr = NULL;
    idx = NULL;
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

    if (ftype.format!=sam && ftype.format!=bam && ftype.format!=cram)
    {
        fprintf(stderr, "[%s:%d %s] Not a SAM/BAM/CRAM file: %s\n", __FILE__, __LINE__, __FUNCTION__, file_name.c_str());
        exit(1);
    }

    if (ref_fasta_file!="")
    {
        char* fai = samfaipath(ref_fasta_file.c_str());
        hts_set_fai_filename(file, fai);
        free(fai);
    }

    hdr = sam_hdr_read(file);
    
    s = bam_init1();

    idx = sam_index_load(file, file_name.c_str());
    if (idx)
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
        itr = sam_itr_querys(idx, hdr, str.s);

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
        itr = sam_itr_querys(idx, hdr, str.s);

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
            if (itr && sam_itr_next(file, itr, s)>=0)
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
        if (sam_read1(file, hdr, s)>=0)
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
    sam_close(file);
}
