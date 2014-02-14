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

#include "bcf_synced_reader.h"

/**
 * Constructor.
 *
 * @intervals - if empty, will add the contigs found in the header files
 */
BCFSyncedReader::BCFSyncedReader(std::vector<std::string>& vcf_files, std::vector<GenomeInterval>& intervals, bool sync_by_pos)
:vcf_files(vcf_files), intervals(intervals), sync_by_pos(sync_by_pos)
{
    nfiles = vcf_files.size();
    vcfs.resize(nfiles, 0);
    hdrs.resize(nfiles, 0);
    idxs.resize(nfiles, 0);
    tbxs.resize(nfiles, 0);
    itrs.resize(nfiles, 0);
    ftypes.resize(nfiles, -1);

    current_interval = "";
    current_pos1 = 0;

    buffer.resize(nfiles);
    s = {0, 0, 0};

    exists_selected_intervals = (intervals.size()!=0);
    for (uint32_t i=0; i<intervals.size(); ++i)
    {
        intervals_map[intervals[i].to_string()] = i;
    }
    intervals_index = 0;

    //1. check file type validity
    //2. loads indices
    //3. adds sequences found in all indexed files, this allows us to iterate through all sequences.
    for (int32_t i = 0; i<nfiles; ++i)
    {
        ftypes[i] = hts_file_type(vcf_files[i].c_str());
        vcfs[i] = bcf_open(vcf_files[i].c_str(), "r");
        hdrs[i] = bcf_alt_hdr_read(vcfs[i]);

        if (i==0)
        {
            if (!(ftypes[i] & (FT_VCF|FT_BCF|FT_STDIN)))
            {
                fprintf(stderr, "[E:%s:%d %s] %s not a VCF or BCF file\n", __FILE__, __LINE__, __FUNCTION__, vcf_files[i].c_str());
                exit(1);
            }

            if (!load_index(i))
            {
                fprintf(stderr, "[I:%s:%d %s] index cannot be loaded for %s\n", __FILE__, __LINE__, __FUNCTION__, vcf_files[i].c_str());
                exit(1);
            }

            if (!exists_selected_intervals)
            {
                //add sequences from file i
                add_interval(i);
            }
        }
        else
        {
            if (!(ftypes[i] & (FT_VCF_GZ|FT_BCF_GZ)))
            {
                fprintf(stderr, "[E:%s:%d %s] %s not a VCF_GZ or BCF file\n", __FILE__, __LINE__, __FUNCTION__, vcf_files[i].c_str());
                exit(1);
            }

            if (!load_index(i))
            {
                fprintf(stderr, "[E:%s:%d %s] index cannot be loaded for %s\n", __FILE__, __LINE__, __FUNCTION__, vcf_files[i].c_str());
                exit(1);
            }

            if (!exists_selected_intervals)
            {
                //add sequences from file i
                add_interval(i);
            }
        }
    }
}

/**
 * Populate sequence names from files.
 * Searches headers first folllowed by tabix.
 */
void BCFSyncedReader::add_interval(int32_t i)
{
    int32_t nseqs = 0;
    const char **seq_names = NULL;

    if (hdrs[i])
    {
        seq_names = bcf_hdr_seqnames(hdrs[i], &nseqs);
        for (uint32_t j=0; j<nseqs; ++j)
        {
            std::string seq(seq_names[j]);
            if (intervals_map.find(seq)==intervals_map.end())
            {
                intervals_map[seq] = intervals_map.size();
                intervals.push_back(GenomeInterval(seq));
            }
        }
        if(seq_names) free(seq_names);
    }

    if (tbxs[i])
    {
        seq_names = tbx_seqnames(tbxs[i], &nseqs);
        for (uint32_t j=0; j<nseqs; ++j)
        {
            std::string seq(seq_names[j]);
            if (intervals_map.find(seq)==intervals_map.end())
            {
                intervals_map[seq] = intervals_map.size();
                intervals.push_back(GenomeInterval(seq));
            }
        }
        if(seq_names) free(seq_names);
    }
}

/**
 * Populate sequence names from files.
 */
void BCFSyncedReader::remove_interval(std::string& interval)
{
    if (intervals_map.find(interval)!=intervals_map.end())
    {
        intervals_map[interval] = intervals_map.size();
        intervals.push_back(GenomeInterval(interval));
    }
}

/**
 * Load index for the ith file, returns true if successful
 */
bool BCFSyncedReader::load_index(int32_t i)
{
    if (ftypes[i]==FT_BCF_GZ)
    {
        if (!(idxs[i] = bcf_index_load(vcf_files[i].c_str())))
        {
            return false;
        }
    }
    else if (ftypes[i]==FT_VCF_GZ)
    {
        if (!(tbxs[i] = tbx_index_load(vcf_files[i].c_str())))
        {
            return false;
        }
    }

    return true;
}

/**
 * Gets sequence name of a record.
 */
const char* BCFSyncedReader::get_seqname(int32_t i, bcf1_t *v)
{
    return bcf_get_chrom(hdrs[i], v);
}

/**
 * Gets number of files read.
 */
int32_t BCFSyncedReader::get_nfiles()
{
    return nfiles;
}

/**
 * Gets current sequence being accessed.
 */
std::string BCFSyncedReader::get_current_sequence()
{
    std::vector<std::string> s;
    split(s, ":", current_interval);
    if (s.size()==0)
    {
        return "";
    }
    else
    {
        return s[0];
    }
}

/**
 * Gets current 1 based position being accessed.
 */
int32_t BCFSyncedReader::get_current_pos1()
{
    return current_pos1;
}

/**
 * Prints buffer.
 */
void BCFSyncedReader::print_buffer()
{
    for (int32_t i = 0; i<nfiles; ++i)
    {
        std::cerr << "#" << i << " ";
        for (std::list<bcf1_t *>::iterator j = buffer[i].begin(); j!=buffer[i].end(); ++j)
        {
            bcf1_t *v = *j;
            std::cerr << " " << v->rid << ":" << (v->pos+1) << ":" << v->d.allele[0] << ":" << v->d.allele[1];
        }
        std::cerr << "\n";
    }
}

/**
 * Closes files.
 */
void BCFSyncedReader::close()
{
    for (int32_t i = 0; i<nfiles; ++i)
    {
        bcf_close(vcfs[i]);
        bcf_hdr_destroy(hdrs[i]);
        bcf_itr_destroy(itrs[i]);
    }
}

/**
 * Inserts a record into pq.
 */
void BCFSyncedReader::insert_into_pq(int32_t i, bcf1_t *v)
{
    pq.push(new bcfptr(i, bcf_get_pos1(v), hdrs[i], v, sync_by_pos));
}

/**
 * Gets record from pool, creates a new record if necessary
 */
bcf1_t* BCFSyncedReader::get_bcf1_from_pool()
{
    if(!pool.empty())
    {
        bcf1_t* v = pool.front();
        pool.pop_front();
        bcf_clear(v);
        return v;
    }
    else
    {
        bcf1_t* v = bcf_init1();
        bcf_clear(v);
        return v;
    }
}

/**
 * Returns record to pool
 */
void BCFSyncedReader::store_bcf1_into_pool(bcf1_t* v)
{
    pool.push_back(v);
    v = 0;
}

int32_t bcfptr_cmp(bcfptr *a, bcfptr *b)
{
    if (a->pos1 == b->pos1)
    {
        if (a->alleles.l!=0 && b->alleles.l!=0)
        {
            int32_t d = strcmp(a->alleles.s, b->alleles.s);
            return d;
        }
        else
        {
            return 0;
        }
    }

    return a->pos1 >= b->pos1 ? 1 : -1;
}

/**
 * Ensures that buffer for each file contains at least records of 2 different positions
 * Updates the latest position.  [store latest and second latest]
 * returns false when all files are read through.
 * Note that these bcf1_t memory allocation are handled by BCFSyncedReader.
 */
bool BCFSyncedReader::read_next_position(std::vector<bcfptr*>& current_recs)
{
    //put records in pool
    for (uint32_t i=0; i<current_recs.size(); ++i)
    {
        store_bcf1_into_pool(current_recs[i]->v);
        delete(current_recs[i]);
    }
    current_recs.clear();

    //process records in priority queue or initialize next interval if pq is empty
    //initialize_next_interval tops up the pq
    if (pq.size()!=0 || initialize_next_interval())
    {
        //dequeue pqueue most recent position and return it
        bcfptr* variant = pq.top();
        bcfptr* cvariant = variant;

        while (bcfptr_cmp(cvariant, variant)==0)
        {
            bcfptr *b = pq.top();
            current_recs.push_back(b);

            buffer[b->file_index].remove(b->v);
            fill_buffer(b->file_index);
            pq.pop();

            if (pq.size()==0)
            {
                break;
            }
            else
            {
                cvariant = pq.top();
            }
        }

        //current_pos1 = current_recs.front().pos1;

        return true;
    }
    else //end of contig or eof for all files
    {
        return false;
    }
}

/**
 * Initialize buffer for next interval.
 * This should only be invoked if the buffer is empty.
 * Returns true if successful.
 */
bool BCFSyncedReader::initialize_next_interval()
{
    while (intervals_index < intervals.size())
    {
        GenomeInterval interval = intervals[intervals_index++];

        for (int32_t i = 0; i<nfiles; ++i)
        {
            int32_t ftype = hts_file_type(vcf_files[i].c_str());
            hts_itr_destroy(itrs[i]);
            itrs[i] = 0;
            interval.to_string(&s);

            if (ftype==FT_BCF_GZ)
            {
                itrs[i] = bcf_itr_querys(idxs[i], hdrs[i], s.s);
            }
            else if (ftype==FT_VCF_GZ)
            {
                itrs[i] = tbx_itr_querys(tbxs[i], s.s);
            }

            fill_buffer(i);
        }

        //make sure pq is not empty
        //it is possible for the pq to be empty as iterators may be returned
        //as the sequence might be a valid sequence stated in the header
        if (pq.size()!=0)
        {
            return true;
        }
    }

    return false;
}

/**
 * Gets records for the most recent position and fills up the buffer from file i.
 * returns true if buffer is filled or it is not necessary to fill buffer.
 * returns false if no more records are found to fill buffer
 */
void BCFSyncedReader::fill_buffer(int32_t i)
{
    //not necessary to fill buffer
    if (buffer[i].size()>=2)
        return;

    int32_t pos1 = buffer[i].size()==0 ? 0 : bcf_get_pos1(buffer[i].front());

    if (ftypes[i]==FT_BCF_GZ)
    {
        bcf1_t *v = get_bcf1_from_pool();
        bool populated = false;
        while (itrs[i] && bcf_itr_next(vcfs[i], itrs[i], v) >= 0)
        {
            populated = true;
            bcf_unpack(v, BCF_UN_STR);
            buffer[i].push_back(v);
            insert_into_pq(i, v);

            if (pos1==0)
            {
                pos1 = bcf_get_pos1(v);
            }

            if (bcf_get_pos1(v)!=pos1)
            {
                break;
            }
            v = get_bcf1_from_pool();
            populated = false;
        }
        if (!populated)
            store_bcf1_into_pool(v);
    }
    else if (ftypes[i]==FT_VCF_GZ)
    {
        while (itrs[i] && tbx_itr_next(vcfs[i], tbxs[i], itrs[i], &s) >= 0)
        {
            bcf1_t *v = get_bcf1_from_pool();
            vcf_parse(&s, hdrs[i], v);

            bcf_unpack(v, BCF_UN_STR);
            buffer[i].push_back(v);
            insert_into_pq(i, v);

            if (pos1==0)
            {
                pos1 = bcf_get_pos1(v);
            }

            if (bcf_get_pos1(v)!=pos1)
            {
                break;
            }
        }
    }
}
