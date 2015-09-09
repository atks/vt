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
    uint32_t tid = bam_get_tid(s);
    uint32_t pos1 = bam_get_pos1(s);
    uint32_t end1 = bam_get_end_pos1(s);

    //wrap bam1_t in AugmentBAMRecord
    as.initialize(h, s);

    //collect statistics for variant records that are in the buffer and overlap with the read
    GenotypingRecord* g;
    for (std::list<GenotypingRecord*>::iterator i=buffer.begin(); i!=buffer.end(); ++i)
    {
        g = *i;
        //same chromosome
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
                collect_sufficient_statistics2(*i, as);
            }
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
        g = new GenotypingRecord(v, vtype);
        buffer.push_back(g);
        
        if (tid==g->rid)
        {    
            if (end1>=g->pos1 && pos1<=g->end1)
            {
                collect_sufficient_statistics2(g, as);
            }
        }
        
        //VCF record occurs after the read
        if (tid<g->rid || end1<g->pos1)
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
void BCFGenotypingBufferedReader::collect_sufficient_statistics2(GenotypingRecord *g, AugmentedBAMRecord& as)
{
    if (g->vtype==VT_SNP)
    {
        if (bcf_get_n_allele(g->v)==2)
        {
            bam1_t *s = as.s;

            char strand = bam_is_rev(s) ? 'R' : 'F';
            char allele = 'R';
            uint32_t pos1 = bam_get_pos1(s);
            uint8_t* seq = bam_get_seq(s);
            uint8_t* qual = bam_get_qual(s);
            uint32_t rlen = bam_get_l_qseq(s);
            uint8_t mapq = bam_get_mapq(s);
            uint32_t q = 30;
            uint32_t cycle = 10;

            std::vector<uint32_t>& aug_cigar = as.aug_cigar;
            std::vector<std::string>& aug_ref = as.aug_ref;
            std::vector<std::string>& aug_alt = as.aug_alt;

            int32_t vpos1 = g->pos1;
            int32_t cpos1 = bam_get_pos1(s);
            int32_t rpos0 = 0;


            for (uint32_t i=0; i<aug_cigar.size(); ++i)
            {
                uint32_t oplen = bam_cigar_oplen(aug_cigar[i]);
                char opchr = bam_cigar_opchr(aug_cigar[i]);

                if (opchr=='S')
                {
                    rpos0 += oplen;
                }
                else if (opchr=='=')
                {
                    if (vpos1>=cpos1 && vpos1<=(cpos1+oplen-1))
                    {
                        rpos0 += vpos1-cpos1;
                        allele = 'R';
                        q = qual[rpos0];
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }
                    else
                    {
                        cpos1 += oplen;
                        rpos0 += oplen;
                    }

                }
                else if (opchr=='X')
                {
                    if (vpos1==cpos1)
                    {
//                        std::cerr << i << ") " << vpos1 << "," << cpos1 << "," << rpos0 << " : " << aug_ref[i] << "/" << aug_alt[i] << " vs " <<bcf_get_allele(g->v)[1][0] << "\n";

                        allele = aug_alt[i].at(0) == bcf_get_allele(g->v)[1][0] ? 'A' : 'O';
                        q = qual[rpos0];

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);

                        break;
                    }

                    ++cpos1;
                    ++rpos0;

                }
                else if (opchr=='I')
                {
                    rpos0 += oplen;
                }
                else if (opchr=='D')
                {
                    cpos1 += oplen;
                }
                else
                {
                    std::cerr << "unrecognized cigar state " << opchr << "\n";
        //            exit(1);
                }
            }

            if (allele=='R')
            {
                if (strand=='F')
                {
                    ++g->depth_fwd;
                    ++g->allele_depth_fwd[0];
                }
                else
                {
                    ++g->depth_rev;
                    ++g->allele_depth_rev[0];
                }
            }
            else //non ref
            {
                ++g->no_nonref;

                if (strand=='F')
                {
                    ++g->depth_fwd;
                    ++g->allele_depth_fwd[1];
                }
                else
                {
                    ++g->depth_rev;
                    ++g->allele_depth_rev[1];
                }
            }

            //update no. of mismatches
            uint32_t no_mismatches = as.no_mismatches;
            if (allele != 'R' && q<20)
            {
                ++no_mismatches;
            }    

            if (allele != 'R' && no_mismatches==0)
            {
                std::cerr << "soemthing wrong\n";
            }    

            g->base_qualities_sum += q;

            ++g->depth;
            g->quals.push_back(q);
            g->map_quals.push_back(mapq);
            g->cycles.push_back(cycle);
            g->alleles.append(1, allele);
            g->strands.append(1, strand);
            g->no_mismatches.push_back(no_mismatches);
        }
        else //multiallelic
        {
        }
    }
    else if (g->vtype==VT_INDEL)
    {
        if (bcf_get_n_allele(g->v)==2)
        {
            bam1_t *s = as.s;

            char strand = bam_is_rev(s) ? 'R' : 'F';
            char allele = 'R';
            uint32_t pos1 = bam_get_pos1(s);
            uint8_t* seq = bam_get_seq(s);
            uint8_t* qual = bam_get_qual(s);
            uint32_t rlen = bam_get_l_qseq(s);
            uint8_t mapq = bam_get_mapq(s);


            int32_t dlen = g->dlen;
            uint32_t len = g->len;
            std::string& indel = g->indel;

            uint32_t q = len*30;
            uint32_t cycle = 10;

            std::vector<uint32_t>& aug_cigar = as.aug_cigar;
            std::vector<std::string>& aug_ref = as.aug_ref;
            std::vector<std::string>& aug_alt = as.aug_alt;

            int32_t vpos1 = g->pos1;

            int32_t cpos1 = bam_get_pos1(s);
            int32_t rpos0 = 0;

            for (uint32_t i=0; i<aug_cigar.size(); ++i)
            {
                uint32_t oplen = bam_cigar_oplen(aug_cigar[i]);
                char opchr = bam_cigar_opchr(aug_cigar[i]);

                if (opchr=='S')
                {
                    rpos0 += oplen;
                }
                else if (opchr=='=')
                {
                    if (vpos1>=cpos1 && vpos1<=(cpos1+oplen-1))
                    {
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                    }

                    cpos1 += oplen;
                    rpos0 += oplen;
                }
                else if (opchr=='X')
                {
                    if (cpos1-1==vpos1)
                    {
                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                    }

                    ++cpos1;
                    ++rpos0;
                }
                else if (opchr=='I')
                {
//                    std::cerr << "dlen  " << dlen << "\n";
//                    std::cerr << "cpos1 " << cpos1 << "\n";
//                    std::cerr << "vpos1 " << vpos1 << "\n";
//                    std::cerr << "indel " << aug_alt[i] << "\n";

                    if (dlen>0 && cpos1-1==vpos1)
                    {
                        if (indel==aug_alt[i])
                        {
                            q = len*30;
                            allele = 'A';
                        }
                        else
                        {
                            q = abs(len-aug_ref[i].size())*30;
                            allele = 'O';
                        }

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }
                    else if (dlen<0 && cpos1-1==vpos1)
                    {
                        q = 30;
                        allele = 'I';

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }

                    rpos0 += oplen;
                }
                else if (opchr=='D')
                {
//                    std::cerr << "dlen  " << dlen << "\n";
//                    std::cerr << "cpos1 " << cpos1 << "\n";
//                    std::cerr << "vpos1 " << vpos1 << "\n";
//                    std::cerr << "indel " << aug_ref[i] << "\n";

                    if (dlen<0 && cpos1-1==vpos1)
                    {
                        if (indel==aug_ref[i])
                        {
                            q = len*30;
                            allele = 'A';
                        }
                        else
                        {
                            q = abs(len-aug_ref[i].size())*30;
                            allele = 'O';
                        }

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }
                    else if (dlen>0 && cpos1-1==vpos1)
                    {
                        q = 30;
                        allele = 'D';

                        cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                        break;
                    }

                    cpos1 += oplen;
                }
                else
                {
                    std::cerr << "unrecognized cigar state " << opchr << "\n";
                }
            }

            if (allele=='R')
            {
                if (strand=='F')
                {
                    ++g->depth_fwd;
                    ++g->allele_depth_fwd[0];
                }
                else
                {
                    ++g->depth_rev;
                    ++g->allele_depth_rev[0];
                }
            }
            else //non ref
            {
                ++g->no_nonref;

                if (strand=='F')
                {
                    ++g->depth_fwd;
                    ++g->allele_depth_fwd[1];
                }
                else
                {
                    ++g->depth_rev;
                    ++g->allele_depth_rev[1];
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

            g->base_qualities_sum += q;

            ++g->depth;
            g->quals.push_back(q);
            g->cycles.push_back(cycle);
            g->alleles.append(1, allele);
            g->strands.append(1, strand);
            g->no_mismatches.push_back(as.no_mismatches);
        }
        else //multiallelic
        {

        }
    }
    else if (g->vtype==VT_VNTR)
    {

    }
}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 */
void BCFGenotypingBufferedReader::collect_sufficient_statistics(GenotypingRecord *g, bam1_t *s)
{
    if (g->vtype==VT_SNP)
    {
        if (bcf_get_n_allele(g->v)==2)
        {
            char strand = bam_is_rev(s) ? 'R' : 'F';
            char allele = 'R';
            uint32_t pos1 = bam_get_pos1(s);
            uint8_t* seq = bam_get_seq(s);
            uint8_t* qual = bam_get_qual(s);
            uint8_t mapq = bam_get_mapq(s);
            //uint32_t qual = 30;
            uint32_t q = 30;
            uint32_t cycle = 10;
            uint8_t *nm_aux;
            int32_t no_mismatches = 0;
            ((nm_aux=bam_aux_get(s, "NM")) &&  (no_mismatches = bam_aux2i(nm_aux)));

            int32_t n_cigar_op = bam_get_n_cigar_op(s);
            int32_t no_mismatches1 = 0;
            if (n_cigar_op)
            {
                int32_t vpos1 = g->pos1;
                int32_t cpos1 = bam_get_pos1(s);
                int32_t rpos0 = 0;

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

                            rpos0 += vpos1-cpos1;

                            std::cerr << bcf_get_allele(g->v)[0][0]  << " vs " << bam_base2char(bam_seqi(bseq, rpos0)) << "\n";

                            char obs_allele = bam_base2char(bam_seqi(bseq, rpos0));

//                            std::cerr << obs_allele << " " << bcf_get_allele(g->v)[0][0] << " " << bcf_get_allele(g->v)[1][0] << "\n";

                            if (obs_allele==bcf_get_allele(g->v)[0][0])
                            {
                                allele = 'R';
                            }
                            else if (obs_allele==bcf_get_allele(g->v)[1][0])
                            {
                                allele = 'A';
                            }
                            else
                            {
                                allele = 'O';
                            }
                            q = bqual[rpos0];
                           // std::cerr << "read length " << bam_get_l_qseq(s) << " " << rpos0 << "\n";
                            cycle = strand == 'F' ? rpos0 : (bam_get_l_qseq(s) - rpos0);

                         //   break;
                        }

                        cpos1 += oplen;
                        rpos0 += oplen;
                    }
                    else if (opchr=='D' || opchr=='N')
                    {
                        cpos1 += oplen;
                    }
                    else if (opchr=='I')
                    {
                        ++no_mismatches1;
                        rpos0 += oplen;
                    }
                    else if (opchr=='S')
                    {
                        rpos0 += oplen;
                    }
                }
            }

            uint8_t *md_aux;
            char* md = 0;
            ((md_aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(md_aux)));
            char* mdp = md;
            bool indel = false;
            uint8_t* bqual = bam_get_qual(s);
            uint32_t rpos0 = 0;
            std::string digit_string;
            uint32_t no_matches = 0;
            while (*mdp)
            {
                if (isdigit(*mdp))
                {
                    str2uint32(digit_string, no_matches);
                    rpos0 += no_matches;
                    indel = false;
                }
                else if (*mdp=='N')
                {
                    ++rpos0;
                    //ignore
                }
                else if (*mdp=='^')
                {
                    ++no_mismatches1;
                    indel = true;
                }
                else //alphabet
                {
                    if (!indel)
                    {
                        ++no_mismatches1;
                    }
                    else
                    {
                        ++rpos0;
                    }
                }

                ++mdp;
            }

            if (false && no_mismatches != no_mismatches1 )
            {
                if (n_cigar_op)
                {
                    uint32_t *cigar = bam_get_cigar(s);
                    for (int32_t i = 0; i < n_cigar_op; ++i)
                    {
                        int32_t opchr = bam_cigar_opchr(cigar[i]);
                        int32_t oplen = bam_cigar_oplen(cigar[i]);

                        std::cerr << oplen << ((char)opchr);
                    }
                    std::cerr << "\t";
                }

                std::cerr << md << " NM=" << no_mismatches << " NMfromMD=" << no_mismatches1 << "\n";
            }

            no_mismatches = no_mismatches1;

            if (allele=='R')
            {
                if (strand=='F')
                {
                    ++g->depth_fwd;
                    ++g->allele_depth_fwd[0];
                }
                else
                {
                    ++g->depth_rev;
                    ++g->allele_depth_rev[0];
                }
            }
            else //non ref
            {
                ++g->no_nonref;

                if (strand=='F')
                {
                    ++g->depth_fwd;
                    ++g->allele_depth_fwd[1];
                }
                else
                {
                    ++g->depth_rev;
                    ++g->allele_depth_rev[1];
                }
            }

            g->base_qualities_sum += q;

            ++g->depth;
            g->quals.push_back(q);
            g->map_quals.push_back(mapq);
            g->cycles.push_back(cycle);
            g->alleles.append(1, allele);
            g->strands.append(1, strand);
            g->no_mismatches.push_back(no_mismatches);
        }
        else //multiallelic
        {
        }
    }
    else if (g->vtype==VT_INDEL)
    {
        if (bcf_get_n_allele(g->v)==2)
        {
            char strand = bam_is_rev(s) ? 'R' : 'F';
            char allele = 'R';
            uint32_t qual = 30;
            uint32_t cycle = 10;
            int32_t no_mismatches = 0;

            //left align CIGAR first.

            int32_t n_cigar_op = bam_get_n_cigar_op(s);

            if (n_cigar_op)
            {
                int32_t vpos1 = g->pos1;
                int32_t cpos1 = bam_get_pos1(s);
                int32_t rpos0 = 0;

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

                            rpos0 += vpos1-cpos1;

                            std::cerr << bcf_get_allele(g->v)[0][0]  << " vs " << bam_base2char(bam_seqi(bseq, rpos0)) << "\n";

                            allele = bam_base2char(bam_seqi(bseq, rpos0)) == bcf_get_allele(g->v)[0][0] ? 'R' : 'A';
                            qual = bqual[rpos0];
                           // std::cerr << "read length " << bam_get_l_qseq(s) << " " << rpos0 << "\n";
                            cycle = strand == 'F' ? rpos0 : (bam_get_l_qseq(s) - rpos0);

                         //   break;
                        }

                        cpos1 += oplen;
                        rpos0 += oplen;
                    }
                    else if (opchr=='D' || opchr=='N')
                    {
                        cpos1 += oplen;
                    }
                    else if (opchr=='I')
                    {
                        ++no_mismatches;
                        rpos0 += oplen;
                    }
                    else if (opchr=='S')
                    {
                        rpos0 += oplen;
                    }
                }
            }

            uint8_t *md_aux;
            char* md = 0;
            ((md_aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(md_aux)));
            char* mdp = md;
            bool indel = false;
            while (*mdp)
            {
                if (isdigit(*mdp))
                {
                    indel = false;
                }
                else if (*mdp=='N')
                {
                    //ignore
                }
                else if (*mdp=='^')
                {
                    ++no_mismatches;
                    indel = true;
                }
                else //alphabet
                {
                    if (!indel) ++no_mismatches;
                }


                ++mdp;
            }

            if (false)
            {
                if (n_cigar_op)
                {
                    uint32_t *cigar = bam_get_cigar(s);
                    for (int32_t i = 0; i < n_cigar_op; ++i)
                    {
                        int32_t opchr = bam_cigar_opchr(cigar[i]);
                        int32_t oplen = bam_cigar_oplen(cigar[i]);

                        std::cerr << oplen << ((char)opchr);
                    }
                    std::cerr << "\t";
                }


                std::cerr << md << " NM=" << no_mismatches << "\n";
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
        else //multiallelic
        {

        }
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
 * Compute SNP genotype likelihoods in PHRED scale.
 */
void BCFGenotypingBufferedReader::compute_snp_pl(std::string& alleles, std::vector<uint32_t>& quals, uint32_t ploidy,  std::vector<uint32_t>& pls)
{
    if (ploidy==2)
    {
        double pRR = 1;
        double pRA = 1;
        double pAA = 1;
        double p;

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            if (alleles[i]=='R')
            {
                p = lt.pl2prob(quals[i]);
                pRR *= 1-p;
                pRA *= 0.5*(1-p+p/3);
                pAA *= p;
            }
            else
            {
                p = lt.pl2prob(quals[i])/3;
                pRR *= p;
                pRA *= 0.5*(1-p+p/3);
                pAA *= 1-p;
            }
        }

        pls[0] = -10*log10(pRR);
        pls[1] = -10*log10(pRA);
        pls[2] = -10*log10(pAA);
    }
}

/**
 * Compute Indel genotype likelihoods in PHRED scale.
 */
void BCFGenotypingBufferedReader::compute_indel_pl(std::string& alleles, std::vector<uint32_t>& quals, uint32_t ploidy,  std::vector<uint32_t>& pls)
{
    if (ploidy==2)
    {
        double pRR = 0;
        double pRA = 0;
        double pAA = 0;
        double p;
        double q;

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            if (alleles[i]=='R')
            {
                p = -((float)quals[i])/10;
                q = p - 1;
                pRR += p;
                pRA += -0.30103+lt.log10sum(p,q);
                pAA += q;
            }
            else
            {
                p = -((float)quals[i])/10;
                q = p - 1;
                pRR += q;
                pRA += -0.30103+lt.log10sum(p,q);
                pAA += p;
            }
        }

        pls[0] = -10*pRR;
        pls[1] = -10*pRA;
        pls[2] = -10*pAA;
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

        std::vector<uint32_t> pls(3);
        compute_snp_pl(g->alleles, g->quals, 2, pls);

        uint32_t min_pl = pls[0];
        uint32_t min_gt_index = 0;
        if (pls[1]<min_pl) {min_pl = pls[1]; min_gt_index = 1;}
        if (pls[2]<min_pl) {min_pl = pls[2]; min_gt_index = 2;}

        pls[0] -= min_pl;
        pls[1] -= min_pl;
        pls[2] -= min_pl;

        int32_t gts[2];
        gts[0] = min_gt_index==2 ? bcf_gt_unphased(1) : bcf_gt_unphased(0);
        gts[1] = min_gt_index==0 ? bcf_gt_unphased(0) : bcf_gt_unphased(1);

        if (g->no_nonref)
        {
            //GT
            bcf_update_genotypes(odw->hdr, v, &gts[0], 2);

            //PL
            bcf_update_format_int32(odw->hdr, v, "PL", &pls[0], 3);

            //depth
            bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);

            //ADF
            bcf_update_format_int32(odw->hdr, v, "ADF", &g->allele_depth_fwd[0], 2);

            //ADR
            bcf_update_format_int32(odw->hdr, v, "ADR", &g->allele_depth_rev[0], 2);

            //base quality
            bcf_update_format_int32(odw->hdr, v, "BQ", &g->quals[0], g->quals.size());

            //map quality
            bcf_update_format_int32(odw->hdr, v, "MQ", &g->map_quals[0], g->map_quals.size());

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
            int32_t gts[2];
            gts[0] = bcf_gt_unphased(0);
            gts[1] = bcf_gt_unphased(0);

            //GT
            bcf_update_genotypes(odw->hdr, v, &gts[0], 2);

            //bq sum
            bcf_update_format_int32(odw->hdr, v, "BQSUM", &g->base_qualities_sum, 1);

            //depth
            bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);

//            //depth
//            bcf_update_format_int32(odw->hdr, v, "DPF", &g->depth_fwd, 1);
//            bcf_update_format_int32(odw->hdr, v, "DPR", &g->depth_rev, 1);
        }

        odw->write(v);
        bcf_destroy(v);

        ++no_snps_genotyped;
    }
    else if (g->vtype==VT_INDEL)
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


        std::vector<uint32_t> pls(3);
        compute_indel_pl(g->alleles, g->quals, 2, pls);

        uint32_t min_pl = pls[0];
        uint32_t min_gt_index = 0;
        if (pls[1]<min_pl) {min_pl = pls[1]; min_gt_index = 1;}
        if (pls[2]<min_pl) {min_pl = pls[2]; min_gt_index = 2;}

        pls[0] -= min_pl;
        pls[1] -= min_pl;
        pls[2] -= min_pl;

        int32_t gts[2];
        gts[0] = min_gt_index==2 ? bcf_gt_unphased(1) : bcf_gt_unphased(0);
        gts[1] = min_gt_index==0 ? bcf_gt_unphased(0) : bcf_gt_unphased(1);


        if (g->no_nonref)
        {
            //GT
            bcf_update_genotypes(odw->hdr, v, &gts[0], 2);

            //PL
            bcf_update_format_int32(odw->hdr, v, "PL", &pls[0], 3);

            //depth
            bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);

            //ADF
            bcf_update_format_int32(odw->hdr, v, "ADF", &g->allele_depth_fwd[0], 2);

            //ADR
            bcf_update_format_int32(odw->hdr, v, "ADR", &g->allele_depth_rev[0], 2);

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
            int32_t gts[2];
            gts[0] = bcf_gt_unphased(0);
            gts[1] = bcf_gt_unphased(0);

            //GT
            bcf_update_genotypes(odw->hdr, v, &gts[0], 2);

            //qualsum
            bcf_update_format_int32(odw->hdr, v, "BQSUM", &g->base_qualities_sum, 1);

            bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);

//            //depth
//            bcf_update_format_int32(odw->hdr, v, "DPF", &g->depth_fwd, 1);
//            bcf_update_format_int32(odw->hdr, v, "DPR", &g->depth_rev, 1);


        }

        odw->write(v);
        bcf_destroy(v);

        ++no_indels_genotyped;
    }
    else if (g->vtype==VT_VNTR)
    {
        ++no_vntrs_genotyped;
    }
}