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
BCFSyncedReader::BCFSyncedReader(std::vector<std::string>& file_names, std::vector<GenomeInterval>& intervals, bool sync_by_pos)
:file_names(file_names), intervals(intervals), sync_by_pos(sync_by_pos)
{
    nfiles = file_names.size();
    files.resize(nfiles, 0);
    hdrs.resize(nfiles, 0);
    idxs.resize(nfiles, 0);
    tbxs.resize(nfiles, 0);
    itrs.resize(nfiles, 0);
    ftypes.resize(nfiles);

    current_interval = "";
    current_pos1 = 0;

    buffer.resize(nfiles);
    s = {0, 0, 0};

    random_access = (intervals.size()!=0);
    for (size_t i=0; i<intervals.size(); ++i)
    {
        intervals_map[intervals[i].to_string()] = i;
    }
    intervals_index = 0;

    uint32_t no_stdins = 0;

//    bool toexit = false;
    for (size_t i = 0; i<nfiles; ++i)
    {
        if (file_names[0]=="+")
        {
            file_names[0]="-";
            ++no_stdins;
        }
        else if (file_names[0]=="-")
        {
            file_names[0]="-";
            ++no_stdins;
        }

        if (no_stdins>1)
        {
            fprintf(stderr, "[E:%s:%d %s] BCFSyncedReader does not support reading from more than one STDIN stream\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }

        files[i] = hts_open(file_names[i].c_str(), "r");
        if (files[i]==NULL)
        {
            fprintf(stderr, "[%s:%d %s] Cannot open %s\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
            exit(1);
//            toexit = true;        
        }
        ftypes[i] = files[i]->format;

        //check format
        if (ftypes[i].format!=vcf && ftypes[i].format!=bcf)
        {
            fprintf(stderr, "[E:%s:%d %s] %s not a VCF or BCF file\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
            exit(1);
        }

        //read header
        hdrs[i] = bcf_alt_hdr_read(files[i]);
        if (!hdrs[i])
        {
            fprintf(stderr, "[E:%s:%d %s] header cannot be read for %s\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
            exit(1);
        }

        //load index if intervals are specified
        if (random_access && !load_index(i))
        {
            fprintf(stderr, "[E:%s:%d %s] index cannot be loaded for %s for random access\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
            exit(1);
        }

        //check contigs consistency
        if (i)
        {
            int32_t nseqs0;
            const char ** seqnames0 = bcf_hdr_seqnames(hdrs[i], &nseqs0);

            int32_t nseqs;
            const char ** seqnames = bcf_hdr_seqnames(hdrs[i], &nseqs);

            if (nseqs0==0 || nseqs==0 || nseqs0!=nseqs)
            {
                fprintf(stderr, "[E:%s:%d %s] contigs in header not consistent with first file for %s\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
                exit(1);
            }

            for (size_t j=0; j<nseqs; ++j)
            {
                if (strcmp(seqnames0[j], seqnames[j]))
                {
                    fprintf(stderr, "[E:%s:%d %s] contigs in header not consistent with first file for %s\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
                    exit(1);
                }
            }

            free(seqnames0);
            free(seqnames);
        }
    }
    
//    if (toexit)
//    {
//        //fprintf(stderr, "[E:%s:%d %s] BCFSyncedReader does not support reading from more than one STDIN stream\n", __FILE__, __LINE__, __FUNCTION__);
//        exit(1);
//        //toexit = true;
//    }
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
        for (size_t j=0; j<nseqs; ++j)
        {
            std::string seq(seq_names[j]);
            if (intervals_map.find(seq)==intervals_map.end())
            {
                intervals_map[seq] = intervals_map.size();
                intervals.push_back(GenomeInterval(seq));
            }
        }
        if (seq_names) free(seq_names);
    }

    if (tbxs[i])
    {
        seq_names = tbx_seqnames(tbxs[i], &nseqs);
        for (size_t j=0; j<nseqs; ++j)
        {
            std::string seq(seq_names[j]);
            if (intervals_map.find(seq)==intervals_map.end())
            {
                intervals_map[seq] = intervals_map.size();
                intervals.push_back(GenomeInterval(seq));
            }
        }
        if (seq_names) free(seq_names);
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
    if (ftypes[i].format==bcf && ftypes[i].compression==bgzf)
    {
        if (!(idxs[i] = bcf_index_load(file_names[i].c_str())))
        {
            return false;
        }
    }
    else if (ftypes[i].format==vcf && ftypes[i].compression==bgzf)
    {
        if (!(tbxs[i] = tbx_index_load(file_names[i].c_str())))
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
    for (size_t i=0; i<nfiles; ++i)
    {
        bcf_close(files[i]);
        bcf_hdr_destroy(hdrs[i]);
        if (idxs[i]) hts_idx_destroy(idxs[i]);
        if (tbxs[i]) tbx_destroy(tbxs[i]);
        bcf_itr_destroy(itrs[i]);
    }

    while (pool.size()!=0)
    {
        bcf_destroy(pool.front());
        pool.pop_front();
    }
}

/**
 * Inserts a record into pq.
 */
void BCFSyncedReader::insert_into_pq(int32_t i, bcf1_t *v)
{
    pq.push(new bcfptr(i, bcf_get_rid(v), bcf_get_pos1(v), hdrs[i], v, sync_by_pos));
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

/**
 * Compares records based on type of comparison.
 */
int32_t BCFSyncedReader::bcfptr_cmp(bcfptr *a, bcfptr *b)
{
    if (random_access)
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
    else
    {
        if (a->rid == b->rid)
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

        return a->rid >= b->rid ? 1 : -1;
    }
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
    for (size_t i=0; i<current_recs.size(); ++i)
    {
        store_bcf1_into_pool(current_recs[i]->v);
        delete current_recs[i];
    }
    current_recs.clear();

    //process records in priority queue or initialize next interval if pq is empty
    //initialize_next_interval tops up the pq
    //initialize_next_interval will never be invoked until the end for non indexed reading
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

        return true;
    }
    else
    {
        //end of contig or eof for all files
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
    if (random_access)
    {
        while (intervals_index < intervals.size())
        {
            GenomeInterval interval = intervals[intervals_index++];

            for (size_t i=0; i<nfiles; ++i)
            {
                hts_itr_destroy(itrs[i]);
                itrs[i] = 0;
                interval.to_string(&s);

                if (ftypes[i].format==bcf)
                {
                    itrs[i] = bcf_itr_querys(idxs[i], hdrs[i], s.s);
                }
                else if (ftypes[i].format==vcf)
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
    else
    {
        for (size_t i=0; i<nfiles; ++i)
        {
            fill_buffer(i);
        }

        if (pq.size()!=0)
        {
            return true;
        }

        return false;
    }
}

/**
 * Gets records for the most recent position and fills up the buffer from file i.
 * returns true if buffer is filled or it is not necessary to fill buffer.
 * returns false if no more records are found to fill buffer
 */
void BCFSyncedReader::fill_buffer(int32_t i)
{
    if (buffer[i].size()>=2)
        return;

    if (random_access)
    {
        int32_t pos1 = buffer[i].size()==0 ? 0 : bcf_get_pos1(buffer[i].front());

        if (ftypes[i].format==bcf)
        {
            bcf1_t *v = get_bcf1_from_pool();
            bool populated = false;

            while (itrs[i] && bcf_itr_next(files[i], itrs[i], v)>=0)
            {
                populated = true;
                bcf_unpack(v, BCF_UN_STR);
                
                //check to ensure order
                if (!buffer[i].empty())
                {
                    if (!bcf_is_in_order(buffer[i].back(), v))
                    {
                        fprintf(stderr, "[E:%s:%d %s] VCF file not in order: %s\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
                        exit(1);
                    }
                }
                
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
        else if (ftypes[i].format==vcf)
        {
            while (itrs[i] && tbx_itr_next(files[i], tbxs[i], itrs[i], &s)>=0)
            {
                bcf1_t *v = get_bcf1_from_pool();
                vcf_parse(&s, hdrs[i], v);

                bcf_unpack(v, BCF_UN_STR);
                
                //check to ensure order
                if (!buffer[i].empty())
                {
                    if (!bcf_is_in_order(buffer[i].back(), v))
                    {
                        fprintf(stderr, "[E:%s:%d %s] VCF file not in order: %s\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
                        exit(1);
                    }
                }
                
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
    else
    {
        int32_t rid = buffer[i].size()==0 ? -1 : bcf_get_rid(buffer[i].front());
        int32_t pos1 = buffer[i].size()==0 ? 0 : bcf_get_pos1(buffer[i].front());

        bcf1_t *v = get_bcf1_from_pool();
        bool populated = false;

        while (bcf_read(files[i], hdrs[i], v)>=0)
        {
            populated = true;
            bcf_unpack(v, BCF_UN_STR);
            
            //check to ensure order
            if (!buffer[i].empty())
            {
                if (!bcf_is_in_order(buffer[i].back(), v))
                {
                    fprintf(stderr, "[E:%s:%d %s] VCF file not in order: %s\n", __FILE__, __LINE__, __FUNCTION__, file_names[i].c_str());
                    exit(1);
                }
            }
            
            buffer[i].push_back(v);
            insert_into_pq(i, v);

            if (rid==-1)
            {
                rid = bcf_get_rid(v);
                pos1 = bcf_get_pos1(v);
            }

            if (bcf_get_rid(v)!=rid || bcf_get_pos1(v)!=pos1)
            {
                break;
            }

            v = get_bcf1_from_pool();
            populated = false;
        }

        if (!populated)
            store_bcf1_into_pool(v);
    }
}
