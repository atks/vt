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
    //stats initialization//
    ////////////////////////
    no_snps_genotyped = 0;
    no_indels_genotyped = 0;
    no_vntrs_genotyped = 0;


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
//    if (buffer.size()>1)
//    {
//        std::cerr << "size : " << buffer.size() << "\n";
//    }

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
            else if (pos1>g->end1)
            {
                //this should not occur if the buffer was flushed before invoking process read
                continue;
            }
            else
            {
//                std::cerr << "\tcollect statistics\n";
                collect_sufficient_statistics(*i, s);
            }
        }
        else if (tid<g->rid)
        {
            return;
        }
        else if (tid==g->rid)
        {
            continue;
        }
    }

    //adding new VCF records
    bcf1_t *v = bcf_init();
    while (odr->read(v))
    {
        //std::cerr << "\tgot read " << "\n";
        int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
        g = new GenotypingRecord(v, vtype);
        
        collect_sufficient_statistics(g, s);
        buffer.push_back(g);

        

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
        
        char strand = bam_is_rev(s) ? 'R' : 'F';
        char allele = 'R';
        uint32_t qual = 30; 
        uint32_t cycle = 10; 
        uint8_t *nm_aux;
        int32_t no_mismatches = 0;
        ((nm_aux=bam_aux_get(s, "NM")) &&  (no_mismatches = bam_aux2i(nm_aux)));
        
        int32_t n_cigar_op = bam_get_n_cigar_op(s);
        if (n_cigar_op)
        {
            int32_t vpos1 = g->pos1;
            int32_t cpos1 = bam_get_pos1(s);
            int32_t rpos1 = 0;
            
            uint32_t *cigar = bam_get_cigar(s);
            for (int32_t i = 0; i < n_cigar_op; ++i)
            {
                int32_t opchr = bam_cigar_opchr(cigar[i]);
                int32_t oplen = bam_cigar_oplen(cigar[i]);
                
                if (opchr=='M')
                {
                    if (vpos1>=cpos1 && vpos1<=cpos1+oplen)
                    {
                        uint8_t* bseq = bam_get_seq(s);
                        uint8_t* bqual = bam_get_qual(s);
                        int32_t l_qseq = bam_get_l_qseq(s);
        
                        rpos1 += vpos1-cpos1;
                        
                        allele = bam_base2char(bam_seqi(bseq, rpos1)) == bcf_get_allele(g->v)[0][0] ? 'R' : 'A';
                        qual = bqual[rpos1];
                       // std::cerr << "read length " << bam_get_l_qseq(s) << " " << rpos1 << "\n";
                        cycle = strand == 'F' ? rpos1 : (bam_get_l_qseq(s) - rpos1);
                    
                        break;   
                    }
                    
                    cpos1 += oplen;
                    rpos1 += oplen;
                }
                else if (opchr=='D' || opchr=='N')
                {
                    cpos1 += oplen;
                }
                else if (opchr=='I')
                {
                    rpos1 += oplen;
                }
            }
        }

        if (allele=='A') 
        {
            ++g->no_nonref;
        }
        
        if (strand=='F') 
        {
            ++g->depth_fwd;
        }
        else
        {
            ++g->depth_rev;
        }
        
        g->base_qualities_sum += qual;

        ++g->depth;
        g->quals.push_back(qual);
        g->cycles.push_back(cycle);
        g->alleles.append(1, allele);
        g->strands.append(1, strand);
        g->no_mismatches.push_back(no_mismatches);
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
                if (bam_get_pos1(s)>g->end1)
                {
                    genotype_and_print(odw, g);
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
                genotype_and_print(odw, g);
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
 * Genotype variant and print to odw.
 */
void BCFGenotypingBufferedReader::genotype_and_print(BCFOrderedWriter* odw, GenotypingRecord* g)
{
    if (g->vtype==VT_SNP)
    {
        bcf1_t *v = bcf_init();
        bcf_clear(v);
        bcf_set_n_sample(v, 1);


        bcf_set_rid(v, bcf_get_rid(g->v));
        bcf_set_pos1(v, bcf_get_pos1(g->v));

        kstring_t new_alleles = {0,0,0};
        char** alleles = bcf_get_allele(g->v);
        for (size_t i=0; i<bcf_get_n_allele(g->v); ++i)
        {
            if (i) kputc(',', &new_alleles);
            kputs(alleles[i], &new_alleles);
        }
        bcf_update_alleles_str(odw->hdr, v, new_alleles.s);

        if (new_alleles.l) free(new_alleles.s);
        
        if (g->no_nonref)
        {
            //bcf_print(odr->hdr, v);
    
            //depth
            bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);    
    
            //quals
            bcf_update_format_int32(odw->hdr, v, "BQ", &g->quals[0], g->quals.size());    
    
            //cycles
            bcf_update_format_int32(odw->hdr, v, "CY", &g->cycles[0], g->cycles.size());    
    
            //strand
            char* str = const_cast<char*>(g->strands.c_str());
            bcf_update_format_string(odw->hdr, v, "ST", const_cast<const char**>(&str), 1);
    
            //alleles
            str = const_cast<char*>(g->alleles.c_str());
            bcf_update_format_string(odw->hdr, v, "AL", const_cast<const char**>(&str),1);
    
            //no of mismatches
            bcf_update_format_int32(odw->hdr, v, "NM", &g->no_mismatches[0], g->no_mismatches.size());   
    
        }
        else
        {
            //depth
            bcf_update_format_int32(odw->hdr, v, "DPF", &g->depth_fwd, 1);    
            bcf_update_format_int32(odw->hdr, v, "DPR", &g->depth_rev, 1);    
    
            //qualsum
            bcf_update_format_int32(odw->hdr, v, "BQSUM", &g->base_qualities_sum, 1);    
        }
        
        odw->write(v);
        bcf_destroy(v);
        
        
        ++no_snps_genotyped;
    }
    else if (g->vtype==VT_INDEL)
    {
        ++no_indels_genotyped;
    }
    else if (g->vtype==VT_VNTR)
    {
        ++no_vntrs_genotyped;
    }
}