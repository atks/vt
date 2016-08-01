/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#include "bcf_single_genotyping_buffered_reader.h"

/**
 * Constructor.
 */
BCFSingleGenotypingBufferedReader::BCFSingleGenotypingBufferedReader(std::string input_vcf_file, std::vector<GenomeInterval>& intervals, std::string output_vcf_file)
{
    /////////////////////
    //io initialization//
    /////////////////////
    odr = new BCFOrderedReader(input_vcf_file, intervals);

    //////////////////////////
    //options initialization//
    //////////////////////////
    output_annotations = false;

    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_snps_genotyped = 0;
    no_indels_genotyped = 0;
    no_vntrs_genotyped = 0;

    ////////////////////////
    //tools initialization//
    ////////////////////////
    vm = new VariantManip();
//    fai = fai_load(ref_fasta_file.c_str());
//    if (fai==NULL)
//    {
//        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
//        exit(1);
//    }
}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 *
 * The VCF records in the buffer must never occur before
 */
void BCFSingleGenotypingBufferedReader::process_read(bam_hdr_t *h, bam1_t *s)
{
    //wrap bam1_t in AugmentBAMRecord
    as.initialize(h, s);

    uint32_t tid = bam_get_tid(s);
    uint32_t beg1 = as.beg1;
    uint32_t end1 = as.end1;

    //collect statistics for variant records that are in the buffer and overlap with the read
    GenotypingRecord* g;
    for (std::list<GenotypingRecord*>::iterator i=buffer.begin(); i!=buffer.end(); ++i)
    {
        g = *i;

//        std::cerr << g->pos1 << " " << g->beg1 << " " << g->end1 << " ";

        //same chromosome
        if (tid==g->rid)
        {
            if (end1 < g->beg1)
            {
                //can terminate
                return;
            }
            else if (beg1 > g->end1)
            {
                //this should not occur if the buffer was flushed before invoking process read
                continue;
            }
            //else if (beg1 <= g->beg1 && g->end1 <= end1)
            else if (beg1 <= g->pos1 && g->pos1 <= end1)
            {
//                collect_sufficient_statistics(*i, as);
            }
            else
            {
                //other types of overlap, just ignore
            }

//            std::cerr << "\n";
        }
        //prior chromosome
        else if (tid<g->rid)
        {
            //this should not occur if the buffer was flushed before invoking process read
            return;
        }
        //latter chromosome
        else if (tid>g->rid)
        {
            //in case if the buffer has a VCF record later in the list which matches it
            continue;
        }
    }

    //you will only reach here if a read occurs after or overlaps the last record in the buffer
    //adding new VCF records and collecting statistics if necessary
    bcf1_t *v = bcf_init();
    while (odr->read(v))
    {
        int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
        g = create_genotyping_record(odr->hdr, v, 2, variant);
        buffer.push_back(g);

        if (tid==g->rid)
        {
            //if (end1>=g->beg1 && pos1<=g->end1)
            if (beg1 <= g->pos1 && g->pos1 <= end1)
            {
//                collect_sufficient_statistics(g, as);
            }
        }

        //VCF record occurs after the read
        if (tid < g->rid || end1 < g->beg1)
        {
            return;
        }
        else
        {
            v = bcf_init();
        }
    }

    //this means end of file
    bcf_destroy(v);
}

/**
 * Flush records.
 */
void BCFSingleGenotypingBufferedReader::flush(bam_hdr_t *h, bam1_t *s, bool flush_all)
{
    if (flush_all)
    {
        //read all the remaining from the reference genotyping file
        bcf1_t *v = bcf_init();
        while (odr->read(v))
        {
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            GenotypingRecord* g = create_genotyping_record(odr->hdr, v, 2, variant);
            buffer.push_back(g);
            v = bcf_init();
        }
        bcf_destroy(v);

        GenotypingRecord* g;
        while (!buffer.empty())
        {
            g = buffer.front();
//            genotype_and_print(g);
            delete g;
            buffer.pop_front();
        }
    }
    else
    {
        //std::cerr << "partial flush\n";

        uint32_t tid = bam_get_tid(s);
        GenotypingRecord* g;

        while (!buffer.empty())
        {
            g = buffer.front();

            if (tid==g->rid)
            {
                if (bam_get_pos1(s) > g->end1)
                {
//                    genotype_and_print(g);
                    delete g;
                    buffer.pop_front();
                }
                else
                {
                    return;
                }
            }
            else if (tid>g->rid)
            {
//                genotype_and_print(g);
                delete g;
                buffer.pop_front();
            }
            else
            {
                return;
            }
        }
    }
}

/**
 * Create appropriate genotyping record.
 */
GenotypingRecord* BCFSingleGenotypingBufferedReader::create_genotyping_record(bcf_hdr_t* h, bcf1_t* v, uint32_t ploidy, Variant& variant)
{
    GenotypingRecord* g = NULL;
    if (variant.type==VT_SNP)
    {
        g = new SNPGenotypingRecord(h, v, 1, 2, NULL);
    }

    return g;
    
}





























