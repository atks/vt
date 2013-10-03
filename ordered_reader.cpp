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

#include "ordered_reader.h"

OrderedReader::OrderedReader(std::string _vcf_file, std::vector<std::string>& _intervals)
{
    vcf_file = _vcf_file;
    intervals = _intervals;
    interval_index = 0;

    vcf = NULL;
    vcfgz = NULL;
    hdr = NULL;
    idx = NULL;
    tbx = NULL;
    itr = NULL;

    ss.str("");
    s.s = 0; s.l = s.m = 0;
    
    //initialize indices
    vcf_ftype = zfile_type(vcf_file.c_str());
    need_random_access = intervals.size()!=0;    

    if (need_random_access)
    {
        if (vcf_ftype==IS_STDIN)
        {
            fprintf(stderr, "[E::%s] Random access is not supported for STDIN\n", __func__);
            exit(1);
        }
        else if (vcf_ftype==IS_VCF)
        {
            fprintf(stderr, "[E::%s] Random access is not supported for non indexed VCF file %s\n", __func__, vcf_file.c_str());
            exit(1);
        }
        else if (vcf_ftype==IS_BCF)
        {
            vcf = vcf_open(vcf_file.c_str(), modify_mode(vcf_file.c_str(), 'r'), 0);
            hdr = vcf_hdr_read(vcf);
            if (!(idx = bcf_index_load(vcf_file.c_str())))
            {
                fprintf(stderr, "[E::%s] fail to load index for %s\n", __func__, vcf_file.c_str());
                exit(1);
            }
        }
        else if (vcf_ftype==IS_VCF_GZ)
        {
            vcf = vcf_open(vcf_file.c_str(), modify_mode(vcf_file.c_str(), 'r'), 0);
            hdr = vcf_hdr_read(vcf);
            vcf_close(vcf);
            vcf = NULL;
            vcfgz = xbgzf_open(vcf_file.c_str(), "r");
            if (!(tbx=tbx_index_load(vcf_file.c_str())))
            {
                fprintf(stderr, "[E::%s] fail to load index for %s\n", __func__, vcf_file.c_str());
                exit(1);
            }
        }
        else
        {
            fprintf(stderr, "[E::%s] %s is not a VCF or BCF file\n", __func__, vcf_file.c_str());
            exit(1);
        }        
    }
    else
    {
        if (vcf_ftype==IS_STDIN || vcf_ftype==IS_VCF || vcf_ftype==IS_BCF || vcf_ftype==IS_VCF_GZ)
        {
            vcf = vcf_open(vcf_file.c_str(), modify_mode(vcf_file.c_str(), 'r'), 0);
            hdr = vcf_hdr_read(vcf);
        }
        else
        {
            fprintf(stderr, "[E::%s] %s is not a VCF or BCF file\n", __func__, vcf_file.c_str());
            exit(1);
        } 
    } 
    
    initialize_next_interval();
};

/**
 * Gets sequence name of a record
 */
const char* OrderedReader::get_seqname(bcf1_t *v)
{
    return bcf_get_chrom(hdr, v);
};                   

/**
 * Initialize next interval.
 * Returns false only if all intervals are accessed.
 */
bool OrderedReader::initialize_next_interval()
{
    if (interval_index==intervals.size())
    {
        return false;
    }    
    
    //update iterators to point at the next region
    ss.str("");
    ss << intervals[interval_index];
    std::cerr << "accessing " << intervals[interval_index] << "\n";
    interval_index++;
    
    //go to next region
    if (vcf_ftype==IS_BCF)
    {
		if (!(itr = bcf_itr_querys(idx, hdr, ss.str().c_str())))
		{
			return initialize_next_interval();
		}
	}
	else //vcf gz
    {
        if (!(itr = tbx_itr_querys(tbx, ss.str().c_str())))
        {
			return initialize_next_interval();
		}
    }
    
    return true;
};

/**
 * Reads next record, hides the random access of different regions from the user.
 */
bool OrderedReader::read1(bcf1_t *v)
{    
    if (need_random_access)
    {
        //go to next region
        if (vcf_ftype==IS_BCF)
        {
            if (bcf_itr_next((BGZF*)vcf->fp, itr, v)>=0)
            {
                return true;
            }    
            else
            {
                if (initialize_next_interval())
                {    
                    return read1(v);
                }
                else
                {
                    return false;
                }
            }
    	}
    	else //vcf gz
        { 
            if (tbx_itr_next(vcfgz, tbx, itr, &s) >= 0)
            {   
                vcf_parse1(&s, hdr, v);
                return true;
            }
            else
            {
                if (initialize_next_interval())
                {    
                    return read1(v);
                }
                else
                {
                    return false;
                }
            }
        }
    }
    else
    {   
        if (vcf_read1(vcf, hdr, v)==0)
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
