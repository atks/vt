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

BCFOrderedWriter::BCFOrderedWriter(std::string input_vcf_file, int32_t window)
{
    this->vcf_file = input_vcf_file;
    this->window = window;
    vcf = NULL;

    s = {0, 0, 0};

    int32_t ftype = hts_file_type(vcf_file.c_str());
    if (!strcmp("+", vcf_file.c_str())) ftype = FT_BCF;

    if (!(ftype & (FT_VCF|FT_BCF|FT_STDIN)))
    {
        fprintf(stderr, "[%s:%d %s] Not a VCF/BCF file: %s\n", __FILE__,__LINE__,__FUNCTION__, vcf_file.c_str());
        exit(1);
    }

    kstring_t *mode = &s;
    kputc('w', mode);
    if (!strcmp("+", vcf_file.c_str()))
    {
        kputs("bu", mode);
        vcf_file = "-";
    }
    if (ftype & FT_BCF) kputc('b', mode);
    if (ftype & FT_GZ) kputc('z', mode);
    vcf = bcf_open(vcf_file.c_str(), mode->s);
    if (vcf==NULL) exit(1);
    
    hdr = bcf_hdr_init("w");
    bcf_hdr_append(hdr, "##fileformat=VCFv4.1");
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
    bcf_hdr_write(vcf, hdr);
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
            //same chromosome?
            if (bcf_get_rid(v)==bcf_get_rid(buffer.back()))
            {
                std::list<bcf1_t*>::iterator i;
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

                //check order
                if (i==buffer.end())
                {
                    int32_t cutoff_pos1 =  std::max(bcf_get_pos1(buffer.front())-window,1);
                    if (bcf_get_pos1(buffer.back())<cutoff_pos1)
                    {
                         std::cerr << "Might not be sorted\n";
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
    }
    else
    {
         bcf_write(vcf, hdr, v);
         store_bcf1_into_pool(v);
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
 * Returns record to pool
 */
void BCFOrderedWriter::store_bcf1_into_pool(bcf1_t* v)
{
    pool.push_back(v);
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
        bcf_clear(v);
        return v;
    }
    else
    {
        bcf1_t *v = bcf_init();
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
            bcf_write(vcf, hdr, buffer.back());
            store_bcf1_into_pool(buffer.back());
            buffer.pop_back();
        }
    }
    else
    {
        if (buffer.size()>=2)
        {
            int32_t cutoff_pos1 =  std::max(bcf_get_pos1(buffer.front())-window,1);

            while (buffer.size()>1)
            {
                if (bcf_get_pos1(buffer.back())<=cutoff_pos1)
                {
                    bcf_write(vcf, hdr, buffer.back());
                    store_bcf1_into_pool(buffer.back());
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
    bcf_close(vcf);
    if (!linked_hdr && hdr) bcf_hdr_destroy(hdr);
}
