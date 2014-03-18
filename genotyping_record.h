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

#ifndef GENOTYPING_RECORD_H
#define GENOTYPING_RECORD_H

#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "bcf_ordered_writer.h"

/**
 * A generic record that holds information for genotyping a 
 * variant across multiple samples. 
 *
 * Maintains read information and allows for additional reads
 * till VCF record can be printed out.
 */
class GenotypingRecord
{
    public:
    bcf_hdr_t *h;
    bcf1_t *v;
    faidx_t *fai;
    std::string error_msg;
        
    /**
     * Constructor.
     * @fai - fai index.
     */
    GenotypingRecord(faidx_t *fai=NULL){};

    /**
     * Constructor.
     * @h - header of the candidate vcf record.
     * @v - candidate VCF record.
     * @fai - fai index.
     */
    GenotypingRecord(bcf_hdr_t *h, bcf1_t *v, faidx_t *fai=NULL){};

    /**
     * Destructor.
     */
    virtual ~GenotypingRecord(){};

    /**
     * Initializes a candidate variant for genotyping.
     */
    virtual bool initialize(bcf_hdr_t *h, bcf1_t *v, faidx_t *fai=NULL)=0;

    /**
     * Initializes a candidate VCF record. Returns false if failure.
     */
    virtual bool set(bcf1_t *v)=0;

    /**
     * Genotypes a read and add to body of evidence.
     */
    virtual void genotype(bam1_t *s)=0;

    /**
     * Prints record.
     */
    virtual void print(BCFOrderedWriter *odw)=0;

    /**
     * Clears this record.
     */
    virtual void clear()=0;
};

#endif