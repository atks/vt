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

#include "synced_reader.h"

/**
 * Constructor.
 */
SyncedReader::SyncedReader(std::vector<std::string>& _vcf_files, std::vector<GenomeInterval>& _intervals)
{
    vcf_files = _vcf_files;
    nfiles = vcf_files.size();
    vcfs.resize(nfiles, NULL);    
    hdrs.resize(nfiles, NULL);
    idxs.resize(nfiles, NULL);
    tbxs.resize(nfiles, NULL);
    itrs.resize(nfiles);
    ftypes.resize(nfiles, -1);
    
    intervals = _intervals;
    interval_index = 0;
    
    indexed_first_file = false;
    current_interval = "";
    current_pos1 = 0;
    
    diff_seq_v = NULL;
    diff_seq_name = "";
    
    buffer.resize(nfiles);
    s = {0, 0, 0};
    
    std::map<char*, uint32_t> sequences;
    bool no_selected_intervals = (intervals.size()==0);
            
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
            
            if (load_index(i))
            {
                indexed_first_file = true;  
            }
            else
            {
                fprintf(stderr, "[I:%s:%d %s] index cannot be loaded for %s\n", __FILE__, __LINE__, __FUNCTION__, vcf_files[i].c_str());
            }
        
            if (!indexed_first_file && no_selected_intervals)
            {
                add_sequence(i);
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
            
            if (!indexed_first_file && no_selected_intervals)
            {
                add_sequence(i);
            }
        }
    }
}

/**
 * Populate sequence names from files.
 */
void SyncedReader::add_sequence(int32_t i)
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
 * Load index for the ith file, returns true if successful
 */
bool SyncedReader::load_index(int32_t i)
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
const char* SyncedReader::get_seqname(int32_t i, bcf1_t *v)
{
    return bcf_get_chrom(hdrs[i], v);
}                   

/**
 * Gets current sequence being accessed.
 */
std::string SyncedReader::get_current_sequence()
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
int32_t SyncedReader::get_current_pos1()
{
    return current_pos1;
}

/**
 * Returns true if there are more contigs to process.
 */
bool SyncedReader::more_intervals()
{
    return interval_index==intervals.size();
}

/**
 * Initialize buffer for next interval.  
 * This should only be invoked if the buffer is empty.
 * Returns true if successful.
 */
bool SyncedReader::initialize_next_interval()
{
    if (indexed_first_file)
    {
        neofs = 0;
    
    	//update iterators to point at the next region
    	if (diff_seq_v==NULL)
    	{
    	    
    	}    
    	
    	GenomeInterval interval(diff_seq_name);
            
	    for (int32_t i = 0; i<nfiles; ++i)
    	{
        	int32_t ftype = hts_file_type(vcf_files[i].c_str());
			hts_itr_destroy(itrs[i]); 
			itrs[i] = 0;
        	
        	if (ftype==FT_BCF_GZ)
        	{
        	    int tid = bcf_hdr_name2id(hdrs[i], interval.seq.c_str());
            	itrs[i] = bcf_itr_queryi(idxs[i], tid, interval.start1, interval.end1+1);
            	if (itrs[i]) 
            	{
            	    ++neofs;
    		    }
			}
			else if (ftype==FT_VCF_GZ)
	    	{
	    	    int tid = tbx_name2id(tbxs[i], interval.seq.c_str());
            	itrs[i] = tbx_itr_queryi(tbxs[i], tid, interval.start1, interval.end1+1);
            	if (itrs[i])
                {
                    ++neofs;
    		    }
	    	}
	    }

    	if (neofs!=nfiles)
    	{
	    	//fill buffer
    		for (int32_t i = 0; i<nfiles; ++i)
    		{
        		fill_buffer(i);
			}
			
			if (pq.size()!=0)
			{
				return true;
			}
		}
        
        
        return true;
    }
    else
    {   
        while (interval_index<intervals.size())
        {
        	neofs = 0;
    
        	//update iterators to point at the next region
        	GenomeInterval interval = intervals[interval_index];
    		interval_index++;
                
    	    for (int32_t i = 0; i<nfiles; ++i)
        	{
            	int32_t ftype = hts_file_type(vcf_files[i].c_str());
    			hts_itr_destroy(itrs[i]); 
    			itrs[i] = 0;
            	
            	if (ftype==FT_BCF_GZ)
            	{
            	    int tid = bcf_hdr_name2id(hdrs[i], interval.seq.c_str());
                	itrs[i] = bcf_itr_queryi(idxs[i], tid, interval.start1, interval.end1+1);
                	if (itrs[i]) 
                	{
                	    ++neofs;
        		    }
    			}
    			else if (ftype==FT_VCF_GZ)
    	    	{
    	    	    int tid = tbx_name2id(tbxs[i], interval.seq.c_str());
                	itrs[i] = tbx_itr_queryi(tbxs[i], tid, interval.start1, interval.end1+1);
                	if (itrs[i])
                    {
                        ++neofs;
        		    }
    	    	}
    	    }
    
        	if (neofs!=nfiles)
        	{
    	    	//fill buffer
        		for (int32_t i = 0; i<nfiles; ++i)
        		{
            		fill_buffer(i);
    			}
    			
    			if (pq.size()!=0)
    			{
    				return true;
    			}
    		}
    	}
    }
	return false;
}

/**
 * Prints buffer.
 */
void SyncedReader::print_buffer()
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
 * Inserts a record into pq.
 */  
void SyncedReader::insert_into_pq(int32_t i, bcf1_t *v)
{
    bcfptr b;
    b.file_index = i;
    b.pos1 = bcf_get_pos1(v);
    b.v = v;
    pq.push(b);
}

/**
 * Gets record from pool, creates a new record if necessary
 */
bcf1_t* SyncedReader::get_bcf1_from_pool()
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
}

/**
 * Returns record to pool 
 */ 
void SyncedReader::store_bcf1_into_pool(bcf1_t* v)
{
    pool.push_back(v);
}

/**
 * Ensures that buffer for each file contains at least records of 2 different positions
 * Updates the latest position.  [store latest and second latest]
 * returns false when all files are read through.  
 * Note that these bcf1_t memory allocation is handled by SyncedReader.
 */
bool SyncedReader::read_next_position(std::vector<bcfptr>& current_recs)
{
    //put records in pool
    for (uint32_t i=0; i<current_recs.size(); ++i)
    {
       store_bcf1_into_pool(current_recs[i].v);
    }
    current_recs.clear();

	//process records in priority queue or invoke next region
    if (pq.size()!=0 || initialize_next_interval())
    {
        //dequeue pqueue most recent position and return it
        int32_t pos1 = pq.top().pos1;
        int32_t cpos1 = pos1;
        
        while (cpos1==pos1)
        {
           bcfptr b;
           b.pos1 = pq.top().pos1;
           b.v = pq.top().v;
           b.file_index = pq.top().file_index;
           
           current_recs.push_back(b);
           pop_and_push_rec(b);
           cpos1 = pq.size()!=0 ? pq.top().pos1 : -1;
        } 
        
        current_pos1 = current_recs.front().pos1;
    
	    return true;
    }
    else //endof contig or eof for all files
    {
        return false;
    }
}

/**
 * Updates pq, buffer simultaneously.
 */
void SyncedReader::pop_and_push_rec(bcfptr b)
{
    buffer[b.file_index].remove(b.v);
    fill_buffer(b.file_index);
    pq.pop();
}

/** 
 * Gets records for the most recent position and fills up the buffer from file i.
 * returns true if buffer is filled or it is not necessary to fill buffer.
 * returns false if no more records are found to fill buffer
 */
void SyncedReader::fill_buffer(int32_t i)
{
    //not necessary to fill buffer
    if (buffer[i].size()>=2)
    {
        return;
    }    
    
    int32_t pos1 = buffer[i].size()==0 ? 0 : bcf_get_pos1(buffer[i].front());
    
    if (!indexed_first_file || i>0)
    {
        if (ftypes[i]==FT_BCF_GZ)
        {   
            bcf1_t *v = get_bcf1_from_pool();
            while (itrs[i] && bcf_itr_next(vcfs[i], itrs[i], v) >= 0)
            {
                bcf_get_pos1(v);
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
            }            
            store_bcf1_into_pool(v); 
        }
        else if (ftypes[i]==FT_VCF_GZ)
        {    
            while (itrs[i] && tbx_itr_next(vcfs[i], tbxs[i], itrs[i], &s) >= 0)
            {
                bcf1_t *v = get_bcf1_from_pool();
                vcf_parse1(&s, hdrs[i], v);
                bcf_unpack(v, BCF_UN_FLT);
                bcf_get_pos1(v);
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
        bcf1_t *v = get_bcf1_from_pool();
        while (bcf_read(vcfs[0], hdrs[0], v) == 0)
        {
            int32_t rid = bcf_get_rid(v);
        
            if (rid==current_rid)
            {
                int32_t pos1 = bcf_get_pos1(v);
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
            else
            {      
                diff_seq_name = std::string(bcf_get_chrom(hdrs[0], v));
                diff_seq_v = v;
            }
            
            v = get_bcf1_from_pool();
        }            
        store_bcf1_into_pool(v); 
    }
}
