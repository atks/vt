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

#include "bcf_ordered_writer.h"

BCFOrderedWriter::BCFOrderedWriter(std::string output_vcf_file_name, int32_t window, int32_t compression)
{
    this->file_name = output_vcf_file_name;
    this->window = window;
    file = NULL;

    kstring_t mode = {0,0,0};
    kputc('w', &mode);

    if (file_name=="+")
    {
        kputs("bu", &mode);
        file_name = "-";
    }
    else if (file_name=="-")
    {
        //do nothing
    }
    else
    {
        if (str_ends_with(file_name, ".vcf"))
        {
            //do nothing
        }
        else if (str_ends_with(file_name, ".vcf.gz"))
        {
            kputc('z', &mode);
            if (compression!=6 && compression!=-1)
            {
                kputw(compression, &mode);
            }
            else if (compression==-1)
            {
                kputw(0, &mode);
            }
        }
        else if (str_ends_with(file_name, ".bcf"))
        {
            kputc('b', &mode);
            if (compression!=6 && compression!=-1)
            {
                kputw(compression, &mode);
            }
            else if (compression==-1)
            {
                kputc('u', &mode);
            }
        }
        else if (str_ends_with(file_name, ".ubcf"))
        {
            kputs("bu", &mode);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Not a VCF/BCF file: %s\n", __FILE__,__LINE__,__FUNCTION__, file_name.c_str());
            exit(1);
        }
    }

    file = bcf_open(file_name.c_str(), mode.s);
    if (file==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot open VCF/BCF file for writing: %s\n", __FILE__,__LINE__,__FUNCTION__, file_name.c_str());
        exit(1);
    }

    hdr = bcf_hdr_init("w");
    bcf_hdr_set_version(hdr, "VCFv4.2");
    linked_hdr = false;
}

/**
 * Duplicates a hdr and sets it.
 */
void BCFOrderedWriter::set_hdr(bcf_hdr_t *hdr)
{
    if (this->hdr)
    {
        bcf_hdr_destroy(this->hdr);
    }
    this->hdr = bcf_hdr_dup(hdr);
    linked_hdr = false;
}

/**
 * Links a header.  This is useful when the VCF file being read has an incomplete header.
 * As the VCF records are read, the incomplete header will be fixed with string type assumptions
 * and the VCF records can be written out without any failure.  The header in the VCF file being
 * written will be incomplete nonetheless and the user should use an alt header when reading the
 * file to bypass the problem.
 */
void BCFOrderedWriter::link_hdr(bcf_hdr_t *hdr)
{
    this->hdr = hdr;
    linked_hdr = true;
}

/**
 * Reads next record, hides the random access of different regions from the user.
 */
void BCFOrderedWriter::write_hdr()
{
    if (bcf_hdr_write(file, hdr))
    {
        fprintf(stderr, "[%s:%d %s] writing of header failed.\n",
                                          __FILE__,
                                          __LINE__,
                                          __FUNCTION__);
        exit(1);
    }
}

/**
 * Reads next record, hides the random access of different regions from the user.
 */
void BCFOrderedWriter::write(bcf1_t *v)
{
    //place into appropriate position in the buffer
    if (window)
    {
        if (!buffer.empty())
        {
            if (bcf_get_rid(v)==bcf_get_rid(buffer.back()))
            {
                std::list<bcf1_t*>::iterator i = buffer.begin();

                for (i=buffer.begin(); i!=buffer.end(); ++i)
                {
                    //equal sign ensures records are kept in original order
                    if (bcf_get_pos1(v)>=bcf_get_pos1(*i))
                    {
                        buffer.insert(i,v);
                        flush(false);
                        return;
                    }
                }

                if (i==buffer.end())
                {
                    int32_t cutoff_pos1 =  std::max((int32_t) bcf_get_pos1(buffer.front())-window,1);
                    if (bcf_get_pos1(v)<cutoff_pos1)
                    {
                        fprintf(stderr, "[%s:%d %s] Might not be sorted for window size %d at current record %s:%d < %d (%d [last record] - %d), please increase window size to at least %d.\n",
                                          __FILE__,
                                          __LINE__,
                                          __FUNCTION__,
                                          window,
                                          bcf_get_chrom(hdr, v),
                                          bcf_get_pos1(v),
                                          (int32_t) cutoff_pos1,
                                          (int32_t) bcf_get_pos1(buffer.front()),
                                          (int32_t) window,
                                          ((int32_t) bcf_get_pos1(buffer.front()))-bcf_get_pos1(v)+1);
                    }
                }

                buffer.insert(i,v);
                flush(false);
            }
            else
            {
                flush(true);
                buffer.push_front(v);
            }
        }
        else
        {
            buffer.push_front(v);
        }

        v = NULL;
    }
    else
    {
        //todo:  add a mechanism to populate header similar to vcf_parse in vcf_format which is called by bcf_write
        if (bcf_write(file, hdr, v))
        {
            fprintf(stderr, "[%s:%d %s] writing of VCF record failed.\n",
                                              __FILE__,
                                              __LINE__,
                                              __FUNCTION__);
            exit(1);
        }
    }
}

/**
 * Flush writable records from buffer.
 */
void BCFOrderedWriter::flush()
{
    flush(true);
}

/**
 * Returns record to pool.
 */
void BCFOrderedWriter::store_bcf1_into_pool(bcf1_t* v)
{
    bcf_clear(v);
    pool.push_back(v);
    v = NULL;
}

/**
 * Gets record from pool, creates a new record if necessary
 */
bcf1_t* BCFOrderedWriter::get_bcf1_from_pool()
{
    if(!pool.empty())
    {
        bcf1_t* v = pool.front();
        pool.pop_front();
        return v;
    }
    else
    {
        bcf1_t* v = bcf_init();
        bcf_clear(v);
        return v;
    }
};

/**
 * Flush writable records from buffer.
 */
void BCFOrderedWriter::flush(bool force)
{
    if (force)
    {
        while (!buffer.empty())
        {
            if (bcf_write(file, hdr, buffer.back()))
            {
                fprintf(stderr, "[%s:%d %s] writing of VCF record failed.\n",
                                                  __FILE__,
                                                  __LINE__,
                                                  __FUNCTION__);
                exit(1);
            }
            bcf_destroy(buffer.back());
            //store_bcf1_into_pool(buffer.back());
            buffer.pop_back();
        }
    }
    else
    {
        if (buffer.size()>1)
        {
            int32_t cutoff_pos1 =  std::max((int32_t) bcf_get_pos1(buffer.front())-window,1);

            while (buffer.size()>1)
            {
                if (bcf_get_pos1(buffer.back())<=cutoff_pos1)
                {
                    if (bcf_write(file, hdr, buffer.back()))
                    {
                        fprintf(stderr, "[%s:%d %s] writing of VCF record failed.\n",
                                                          __FILE__,
                                                          __LINE__,
                                                          __FUNCTION__);
                        exit(1);
                    }
                    bcf_destroy(buffer.back());
                    //store_bcf1_into_pool(buffer.back());
                    buffer.pop_back();
                }
                else
                {
                    return;
                }
            }
        }
    }
}

/**
 * Closes the file.
 */
void BCFOrderedWriter::close()
{
    flush(true);
    bcf_close(file);
    if (!linked_hdr && hdr) bcf_hdr_destroy(hdr);
//    while (buffer.size()!=0)
//    {
//        bcf_destroy(buffer.pop_back());
//        buffer.pop_back();
//        if (bcf_get_pos1(buffer.back())<=cutoff_pos1)
//        {
//            bcf_write(file, hdr, buffer.back());
//            store_bcf1_into_pool(buffer.back());
//            buffer.pop_back();
//        }
//        else
//        {
//            return;
//        }
//    }
}
