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

#include "bcf_genotyping_buffered_reader.h"

/**
 * Constructor.
 */
BCFGenotypingBufferedReader::BCFGenotypingBufferedReader(std::string filename, std::vector<GenomeInterval>& intervals)
{
    /////////////////////
    //io initialization//
    /////////////////////

    odr = new BCFOrderedReader(filename, intervals);
    
    ////////////////////////
    //tools initialization//
    ////////////////////////
    vm = new VariantManip();

}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 *
 * The VCF records in the buffer must never occur before
 */
void BCFGenotypingBufferedReader::process_read(bam_hdr_t *h, bam1_t *s)
{
    uint32_t tid = bam_get_tid(s);
    uint32_t pos1 = bam_get_pos1(s);
    uint32_t end1 = bam_get_end_pos1(s);
    
    GenotypingRecord* g;            
    for (std::list<GenotypingRecord*>::iterator i=buffer.begin(); i!=buffer.end(); ++i)
    {
        g = *i;
        if (tid==g->rid)
        {   
            if (end1<g->pos1)
            {
                //can terminate
                return;
            }
            else if  (pos1>g->end1)
            {
                //this should not occur if the buffer was flushed before invoking process read
                continue;
            }
            else
            {
                collect_sufficient_statistics(*i, s);
            }    
        }
        else 
        {
            return;
        }    
    }
    
    bcf1_t *v = bcf_init();
    while (odr->read(v))
    {    
        int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
        g = new GenotypingRecord(v, vtype);
        
        if (tid!=g->rid || end1<g->pos1)
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
 * Collects sufficient statistics from read for variants to be genotyped.
 */
void BCFGenotypingBufferedReader::collect_sufficient_statistics(GenotypingRecord *g, bam1_t *s)
{
    if (g->vtype==VT_SNP)
    {
        int32_t vpos1 = g->pos1;
        uint32_t pos1 = bam_get_pos1(s);


        uint8_t* seq = bam_get_seq(s);
        
        uint8_t* qual = bam_get_qual(s);
        int32_t l_qseq = bam_get_l_qseq(s);
        uint32_t* cigar = bam_get_cigar(s);
        char strand = bam_is_rev(s) ? '-' : '+';
        
        
        //get base
        /**
         * Gets the end position of the last mapped base in the read.
         */
            int32_t end_pos1 = bam_get_pos1(s);
            int32_t n_cigar_op = bam_get_n_cigar_op(s);
            if (n_cigar_op)
            {
                uint32_t *cigar = bam_get_cigar(s);
                for (int32_t i = 0; i < n_cigar_op; ++i)
                {
                    int32_t opchr = bam_cigar_opchr(cigar[i]);
                    
                    if (opchr=='M' || opchr=='D' || opchr=='N' || opchr=='=' || opchr=='X')
                    {
                        end_pos1 += bam_cigar_oplen(cigar[i]);
                    }
                }
            }
            
        //    return end_pos1-1;
        
    
        
        
            
        //nm    
        uint8_t *nm_aux;
        int32_t nm;
        if ((nm_aux=bam_aux_get(s, "NM")) &&  (nm = bam_aux2i(nm_aux)))
        {    
            g->no_mismatches.push_back(nm);
        }
            
    }   
    else if (g->vtype==VT_INDEL) 
    {
        
    }
    else if (g->vtype==VT_VNTR) 
    {
        
    }
}

/**
 * Flush records.
 */
void BCFGenotypingBufferedReader::flush(BCFOrderedWriter* odw, bam_hdr_t *h, bam1_t *s, bool flush_all)
{
    if (flush_all)
    {
        GenotypingRecord* g;            
        while (!buffer.empty())
        {
            g = buffer.front();
            genotype_and_print(odw, g);
            buffer.pop_front();
        }
    }    
    else
    {
        uint32_t tid = bam_get_tid(s);
        GenotypingRecord* g;            
        
        while (!buffer.empty())
        {
            g = buffer.front();
            
            std::cerr << " " << buffer.size() << "\n";
            if (tid==g->rid)
            {
                if (bam_get_end_pos1(s)<g->pos1)
                {
                    genotype_and_print(odw, g);
                    buffer.pop_front();
                }
                else 
                {
                    return;
                }    
            }
            else if (tid<g->rid)
            {
                genotype_and_print(odw, g);
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
 * Genotype variant and print to odw.
 */
void BCFGenotypingBufferedReader::genotype_and_print(BCFOrderedWriter* odw, GenotypingRecord* g)
{
    if (g->vtype==VT_SNP)
    {
        bcf1_t *v = bcf_init();
        
       
        bcf_set_rid(v, bcf_get_rid(g->v));
        
        kstring_t new_alleles = {0,0,0};    
        char** alleles = bcf_get_allele(g->v);
        for (size_t i=0; i<bcf_get_n_allele(g->v); ++i)
        {
            if (i) kputc(',', &new_alleles);
            kputs(alleles[i], &new_alleles);
        }
        bcf_update_alleles_str(odw->hdr, v, new_alleles.s);
        
        if (new_alleles.l) free(new_alleles.s);
        
        odw->write(g->v);
    }   
    else if (g->vtype==VT_INDEL) 
    {
        
    }
    else if (g->vtype==VT_VNTR) 
    {
        
    }
}