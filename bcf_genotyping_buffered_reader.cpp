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
}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 *
 * The VCF records in the buffer must never occur before
 */
void BCFGenotypingBufferedReader::process_read(bam_hdr_t *h, bam1_t *s)
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
                collect_sufficient_statistics(*i, as);

                if (beg1 <= g->beg1 && g->end1 <= end1)
                {
//                    std::cerr << "COMPLETE";
                }
                else
                {
                    //drop
//                    std::cerr << "PARTIAL";
                }
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
        g = new GenotypingRecord(odr->hdr, v, vtype);
        buffer.push_back(g);

        if (tid==g->rid)
        {
            //if (end1>=g->beg1 && pos1<=g->end1)
            if (beg1 <= g->beg1 && g->end1 <= end1)
            {
                collect_sufficient_statistics(g, as);
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
 * Collects sufficient statistics from read for variants to be genotyped.
 */
void BCFGenotypingBufferedReader::collect_sufficient_statistics(GenotypingRecord *g, AugmentedBAMRecord& as)
{
    if (g->vtype==VT_SNP)
    {
        if (bcf_get_n_allele(g->v)==2)
        {
            bam1_t *s = as.s;

            char strand = bam_is_rev(s) ? 'R' : 'F';
            int32_t allele = 0;
            uint32_t pos1 = bam_get_pos1(s);
            uint8_t* seq = bam_get_seq(s);
            uint8_t* qual = bam_get_qual(s);
            int32_t rlen = bam_get_l_qseq(s);
            uint8_t mapq = bam_get_mapq(s);
            uint32_t q = 30;
            int32_t cycle = 0;

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
                        allele = 0;
                        q = qual[rpos0];
                        cycle = rpos0<(rlen>>1) ? (rpos0+1) : -(rlen - rpos0 + 1);

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

                        allele = (aug_alt[i].at(0) == bcf_get_allele(g->v)[1][0]) ? 1 : -1;
                        q = qual[rpos0];
                        cycle = rpos0<(rlen>>1) ? (rpos0+1) : -(rlen - rpos0 + 1);

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

            if (allele==0)
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
            if (allele!=0 && q<20)
            {
                ++no_mismatches;
            }

            if (allele!=0 && no_mismatches==0)
            {
                std::cerr << "something wrong\n";
            }

            g->base_qualities_sum += q;

            ++g->depth;
            g->quals.push_back(q);
            g->map_quals.push_back(mapq);
            g->cycles.push_back(cycle);
            g->alleles.push_back(allele);
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
            if (as.beg1 <= g->beg1 && g->end1 <= as.end1)
            {
                bam1_t *s = as.s;

                char strand = bam_is_rev(s) ? 'R' : 'F';
                int32_t allele = 0;
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
                                allele = 1;
                            }
                            else
                            {
                                q = abs(len-aug_ref[i].size())*30;
                                allele = -1;
                            }

                            cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                            break;
                        }
                        else if (dlen<0 && cpos1-1==vpos1)
                        {
                            q = 30;
                            allele = -3;

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
                                allele = 1;
                            }
                            else
                            {
                                q = abs(len-aug_ref[i].size())*30;
                                allele = -1;
                            }

                            cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                            break;
                        }
                        else if (dlen>0 && cpos1-1==vpos1)
                        {
                            q = 30;
                            allele = -2;

                            cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
                            break;
                        }

                    }
                    else
                    {
                        std::cerr << "unrecognized cigar state " << opchr << "\n";
                    }
                }

                if (allele==0)
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

                if (allele!=0)
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
                g->map_quals.push_back(mapq);
                g->alleles.push_back(allele);
                g->strands.append(1, strand);
                g->no_mismatches.push_back(as.no_mismatches);
            }
            else
            {
                bam1_t *s = as.s;
                uint8_t mapq = bam_get_mapq(s);

                ++g->depth;
                g->quals.push_back(0);
                g->cycles.push_back(-1);
                g->map_quals.push_back(mapq);
                g->alleles.push_back(-1);
                g->strands.append(1, (bam_is_rev(as.s) ? 'R' : 'F'));
                g->no_mismatches.push_back(as.no_mismatches);
            }
        }
        else //multiallelic
        {
        }
    }
    else if (g->vtype==VT_VNTR)
    {
        bam1_t *s = as.s;

        char strand = bam_is_rev(s) ? 'R' : 'F';
        int32_t allele = 0;
        uint32_t pos1 = bam_get_pos1(s);
        uint8_t* seq = bam_get_seq(s);
        uint8_t* qual = bam_get_qual(s);
        uint32_t rlen = bam_get_l_qseq(s);
        uint8_t mapq = bam_get_mapq(s);

        std::vector<uint32_t>& aug_cigar = as.aug_cigar;
        std::vector<std::string>& aug_ref = as.aug_ref;
        std::vector<std::string>& aug_alt = as.aug_alt;

        //genomic bookend positions of VNTR
        int32_t vpos1 = g->beg1-1;
        int32_t vend1 = g->end1+1;

        //position with respect to read
        int32_t cpos1 = bam_get_pos1(s);
        int32_t rpos0 = 0;

        //genomic bookend positions of VNTR translated to read position
        int32_t pos0 = -1;
        int32_t end0 = -1;


        //locate repeat region on read.
        //translating genomic coordinates to read positions
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
                if (pos0==-1 && (cpos1<=vpos1 && vpos1<=(cpos1+oplen-1)))
                {
                    pos0 = rpos0 + (vpos1-cpos1+1);
                }

                if (end0==-1 && (cpos1<=vend1 && vend1<=(cpos1+oplen-1)))
                {
                    end0 = rpos0 + (vend1-cpos1+1);
                    break;
                }

                cpos1 += oplen;
                rpos0 += oplen;
            }
            else if (opchr=='X')
            {
                if (pos0==-1 && (cpos1==vpos1))
                {
                    pos0 = rpos0;
                }

                if (end0==-1 && (cpos1==vend1))
                {
                    end0 = rpos0;
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
            }
        }

        //compute repeat tract units
        float counts = 0;

//        std::cerr << "pos0,end0 = " << pos0 << "," <<end0 <<  " (" << g->motif.size() <<")\n";

        if (pos0!=-1 && end0!=-1)
        {
            counts = ((float)(end0-pos0+1))/g->motif.size();
        }

        if (counts)
        {
            //update genotyping record
            ++g->depth;
            g->counts.push_back(counts);
            g->strands.append(1, strand);
            g->no_mismatches.push_back(as.no_mismatches);
        }
    }
}

/**
 * Flush records.
 */
void BCFGenotypingBufferedReader::flush(BCFOrderedWriter* odw, bam_hdr_t *h, bam1_t *s, bool flush_all)
{
    if (flush_all)
    {
        //read all the remaining from the reference genotyping file
        bcf1_t *v = bcf_init();
        while (odr->read(v))
        {
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            GenotypingRecord* g = new GenotypingRecord(odr->hdr, v, vtype);
            buffer.push_back(g);
            v = bcf_init();
        }
        bcf_destroy(v);

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
                if (bam_get_pos1(s) > g->end1)
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

void BCFGenotypingBufferedReader::compute_snp_pl(std::vector<int32_t>& alleles, std::vector<uint32_t>& quals, uint32_t ploidy,  std::vector<uint32_t>& pls) 
{
//void BCFGenotypingBufferedReader::compute_snp_pl(std::vector<int32_t>& alleles, std::vector<uint32_t>& quals, uint32_t ploidy, uint32_t no_alleles, std::vector<uint32_t>& pls) //, float& pl_offset)
    int no_alleles = 2;
    if (ploidy==2 && no_alleles==2)
    {
//        double pRR = 1;
//        double pRA = 1;
//        double pAA = 1;
//        double p = 0;
//
//        for (uint32_t i=0; i<alleles.size(); ++i)
//        {
//            if (alleles[i]==0)
//            {
////                actual
////                p = lt.pl2prob(quals[i]);
////                pRR *= 1-p;
////                pRA *= 0.5*(1-p+p/3);
////                pAA *= p/3;
//
//                //simplified
//                p = lt.pl2prob(quals[i]);
//                pRR *= 1-p;
//                p /= 3;
//                pRA *= 0.5-p;
//                pAA *= p;
//            }
//            else if (alleles[i]==1)
//            {
////                actual
////                p = lt.pl2prob(quals[i]);
////                pRR *= p/3;
////                pRA *= 0.5*(1-p+p/3);
////                pAA *= 1-p;
//
//                //simplified
//                p = lt.pl2prob(quals[i])/3;
//                pRR *= p;
//                pRA *= 0.5-p;
//                pAA *= 1-3*p;
//            }
//            else if (alleles[i]<=-1)
//            {
////                actual
////                p = lt.pl2prob(quals[i]);
////                pRR *= p/3;
////                pRA *= 0.5*(2*p/3);
////                pAA *= p/3;
//
//                //simplified
////                p = lt.pl2prob(quals[i])/3;
////                pRR *= p;
////                pRA *= p;
////                pAA *= p;
//                //don't anything actually
//            }
//            else //deletion
//            {
//                //ignore this for the time being
//            }
//
//            double total = pRR + pRA + pAA;
//
//            pRR /= total;
//            pRA /= total;
//            pAA /= total;
//        }
//
//        pls[0] = -10*log10(pRR);
//        pls[1] = -10*log10(pRA);
//        pls[2] = -10*log10(pAA);


        double pRR = 0;
        double pRA = 0;
        double pAA = 0;
        double p = 0;

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            if (quals[i]==0) continue;
                
            if (alleles[i]==0)
            {
                p = lt.pl2prob(quals[i]);
                pRR += log10(1-p);
                p /= 3;
                pRA += log10(0.5-p);
                pAA += log10(p);


            }
            else if (alleles[i]==1)
            {
                p = lt.pl2prob(quals[i])/3;
                pRR += log10(p);
                pRA += log10(0.5-p);
                pAA += log10(1-3*p);
            }
            else
            {
                //do nothing
            }
        }
        
        pls[0] = -10*pRR;
        pls[1] = -10*pRA;
        pls[2] = -10*pAA;
  }
    else if (ploidy==2 && no_alleles>2)
    {
        int32_t no_genotypes = (no_alleles * (no_alleles+1)) >> 1;

        float pG[no_genotypes];

        for (uint32_t i = 0; i<no_genotypes; ++i)
        {
            pG[i] = 1;
        }

        float p;

        float offset = 0;

        //handle multiallelics
        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            p = lt.pl2prob(quals[i]);

            for (uint32_t g = 0; g<no_genotypes; ++g)
            {
//                if ()

                pG[i] = 1;
            }
//            pG[alleles[i]] *=  ;


        }
    }
    else //generic number of ploidy and alleles
    {
        //http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
    }
}



//void BCFGenotypingBufferedReader::compute_snp_pl(std::vector<int32_t>& alleles, std::vector<uint32_t>& quals, uint32_t ploidy,  std::vector<uint32_t>& pls) {
//    if (ploidy==2)
//    {
//        double pRR = 1;
//        double pRA = 1;
//        double pAA = 1;
//        double p = 0;
//
//        for (uint32_t i=0; i<alleles.size(); ++i)
//        {
//            if (quals[i]==0) continue;
//            
//            if (alleles[i]==0)
//            {
//                p = lt.pl2prob(quals[i]);
//                pRR *= 1-p;
//                pRA *= 0.5*(1-p+p/3);
//                pAA *= p;
//            }
//            else if (alleles[i]==1)
//            {
//                p = lt.pl2prob(quals[i])/3;
//                pRR *= p;
//                pRA *= 0.5*(1-p+p/3);
//                pAA *= 1-p;
//            }
//            else if (alleles[i]==-1)
//            {
//                p = lt.pl2prob(quals[i])/3;
//                pRR *= p;
//                pRA *= 0.5*(2*p/3);
//                pAA *= p;
//            }
//            
////            std::cerr << alleles[i] << " " << quals[i] <<  "\n" ;
////            std::cerr << "\t" << -10*log10(pRR) << " " << -10*log10(pRA) << " " << -10*log10(pAA) << "\n" ;
//        }
//
////        std::cerr << -10*log10(pRR) << " " << -10*log10(pRA) << " " << -10*log10(pAA) << "\n" ;
//        pls[0] = -10*log10(pRR);
//        pls[1] = -10*log10(pRA);
//        pls[2] = -10*log10(pAA);
//    }
//}




/**
 * Compute Indel genotype likelihoods in PHRED scale.
 */
void BCFGenotypingBufferedReader::compute_indel_pl(std::vector<int32_t>& alleles, std::vector<uint32_t>& quals, uint32_t ploidy,  std::vector<uint32_t>& pls)
{
    if (ploidy==2)
    {
        double pRR = 0;
        double pRA = 0;
        double pAA = 0;
        double p = 0;
        double q = 0;

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            if (alleles[i]==0)
            {
                p = -((float)quals[i])/10;
                q = p - 2;
                pRR += p;
                pRA += -0.30103+lt.log10sum(p,q);
                pAA += q;
            }
            else if (alleles[i]==1)
            {
                p = -((float)quals[i])/10;
                q = p - 2;
                pRR += q;
                pRA += -0.30103+lt.log10sum(p,q);
                pAA += p;
            }
            else
            {
                p = -((float)quals[i])/10;
                q = p - 2;
                p = q;
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

        //copy over values
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_snp"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_snp"));
        }
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_indel"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_indel"));
        }
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_vntr"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr"));
        }

        if (output_annotations)
        {
            char* flankseq = NULL;
            int32_t n = 0;
            if (bcf_get_info_string(odr->hdr, g->v, "FLANKSEQ", &flankseq, &n)>0)
            {
                bcf_update_info_string(odw->hdr, v, "FLANKSEQ", flankseq);
                free(flankseq);
            }
        }

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
            bcf_update_format_int32(odw->hdr, v, "AL", &g->alleles[0], g->alleles.size());

            //no of mismatches
            bcf_update_format_int32(odw->hdr, v, "NM", &g->no_mismatches[0], g->no_mismatches.size());
        }
        else
        {
            if (g->depth)
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
            }
            else
            {
                int32_t gts[2];
                gts[0] = bcf_gt_unphased(-1);
                gts[1] = bcf_gt_unphased(-1);

                //GT
                bcf_update_genotypes(odw->hdr, v, &gts[0], 2);

                //bq sum
                bcf_update_format_int32(odw->hdr, v, "BQSUM", &g->base_qualities_sum, 1);

                //depth
                bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);
            }
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

        //copy over values
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_snp"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_snp"));
        }
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_indel"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_indel"));
        }
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_vntr"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr"));
        }

        int32_t* flanks = NULL;
        int32_t n = 0;
        if (bcf_get_info_int32(odr->hdr, g->v, "FLANKS", &flanks, &n)>0)
        {
            bcf_update_info_int32(odw->hdr, v, "FLANKS", flanks, 2);
            free(flanks);
        }
        int32_t* fz_flanks = NULL;
        n = 0;
        if (bcf_get_info_int32(odr->hdr, g->v, "FZ_FLANKS", &fz_flanks, &n)>0)
        {
            bcf_update_info_int32(odw->hdr, v, "FZ_FLANKS", fz_flanks, 2);
            free(fz_flanks);
        }

        if (output_annotations)
        {
            char* motif = NULL;
            n =0;
            if (bcf_get_info_string(odr->hdr, g->v, "TR", &motif, &n)>0)
            {
                bcf_update_info_string(odw->hdr, v, "TR", motif);
                free(motif);
            }
            char* flankseq = NULL;
            n = 0;
            if (bcf_get_info_string(odr->hdr, g->v, "FLANKSEQ", &flankseq, &n)>0)
            {
                bcf_update_info_string(odw->hdr, v, "FLANKSEQ", flankseq);
                free(flankseq);
            }
            char* tr = NULL;
            n = 0;
            if (bcf_get_info_string(odr->hdr, g->v, "TR", &tr, &n)>0)
            {
                bcf_update_info_string(odw->hdr, v, "TR", tr);
                free(tr);
            }
        }

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

            //AD
            uint32_t allele_depth[2] = {g->allele_depth_fwd[0]+g->allele_depth_rev[0], g->allele_depth_fwd[1]+g->allele_depth_rev[1]};
            bcf_update_format_int32(odw->hdr, v, "AD", &allele_depth, 2);

            //ADF
            bcf_update_format_int32(odw->hdr, v, "ADF", &g->allele_depth_fwd[0], 2);

            //ADR
            bcf_update_format_int32(odw->hdr, v, "ADR", &g->allele_depth_rev[0], 2);

            //quals
            bcf_update_format_int32(odw->hdr, v, "BQ", &g->quals[0], g->quals.size());

            //map quality
            bcf_update_format_int32(odw->hdr, v, "MQ", &g->map_quals[0], g->map_quals.size());

            //cycles
            bcf_update_format_int32(odw->hdr, v, "CY", &g->cycles[0], g->cycles.size());

            //strand
            char* str = const_cast<char*>(g->strands.c_str());
            bcf_update_format_string(odw->hdr, v, "ST", const_cast<const char**>(&str), 1);

            //alleles
            bcf_update_format_int32(odw->hdr, v, "AL", &g->alleles[0], g->alleles.size());

            //no of mismatches
            bcf_update_format_int32(odw->hdr, v, "NM", &g->no_mismatches[0], g->no_mismatches.size());
        }
        else
        {
            if (g->depth)
            {
                int32_t gts[2];
                gts[0] = bcf_gt_unphased(0);
                gts[1] = bcf_gt_unphased(0);

                //GT
                bcf_update_genotypes(odw->hdr, v, &gts[0], 2);

                //qualsum
                bcf_update_format_int32(odw->hdr, v, "BQSUM", &g->base_qualities_sum, 1);

                bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);
            }
            else
            {
                int32_t gts[2];
                gts[0] = bcf_gt_unphased(-1);
                gts[1] = bcf_gt_unphased(-1);

                //GT
                bcf_update_genotypes(odw->hdr, v, &gts[0], 2);

                //bq sum
                bcf_update_format_int32(odw->hdr, v, "BQSUM", &g->base_qualities_sum, 1);

                //depth
                bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);
            }

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

        //copy over values
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_snp"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_snp"));
        }
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_indel"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_indel"));
        }
        if (bcf_has_filter(odr->hdr, g->v, const_cast<char*>("overlap_vntr"))==1)
        {
            bcf_add_filter(odw->hdr, v, bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr"));
        }

        char* motif = NULL;
        int32_t n =0;
        if (bcf_get_info_string(odr->hdr, g->v, "MOTIF", &motif, &n)>0)
        {
            bcf_update_info_string(odw->hdr, v, "MOTIF", motif);
            free(motif);
        }
        char* ru = NULL;
        n =0;
        if (bcf_get_info_string(odr->hdr, g->v, "RU", &ru, &n)>0)
        {
            bcf_update_info_string(odw->hdr, v, "RU", ru);
            free(ru);
        }
        float* fz_concordance = NULL;
        n =0;
        if (bcf_get_info_float(odr->hdr, g->v, "FZ_CONCORDANCE", &fz_concordance, &n)>0)
        {
            bcf_update_info_float(odw->hdr, v, "FZ_CONCORDANCE", fz_concordance, 1);
            free(fz_concordance);
        }
        float* fz_rl = NULL;
        n =0;
        if (bcf_get_info_float(odr->hdr, g->v, "FZ_RL", &fz_rl, &n)>0)
        {
            bcf_update_info_float(odw->hdr, v, "FZ_RL", fz_rl, 1);
            free(fz_rl);
        }
        float* fz_ll = NULL;
        n =0;
        if (bcf_get_info_float(odr->hdr, g->v, "FZ_LL", &fz_ll, &n)>0)
        {
            bcf_update_info_float(odw->hdr, v, "FZ_LL", fz_ll, 1);
            free(fz_ll);
        }

        int32_t* flanks = NULL;
        n = 0;
        if (bcf_get_info_int32(odr->hdr, g->v, "FLANKS", &flanks, &n)>0)
        {
            bcf_update_info_int32(odw->hdr, v, "FLANKS", flanks, 2);
            free(flanks);
        }
        int32_t* fz_flanks = NULL;
        n = 0;
        if (bcf_get_info_int32(odr->hdr, g->v, "FZ_FLANKS", &fz_flanks, &n)>0)
        {
            bcf_update_info_int32(odw->hdr, v, "FZ_FLANKS", fz_flanks, 2);
            free(fz_flanks);
        }
        char* flankseq = NULL;
        n = 0;


        if (output_annotations)
        {
            if (bcf_get_info_string(odr->hdr, g->v, "FLANKSEQ", &flankseq, &n)>0)
            {
                bcf_update_info_string(odw->hdr, v, "FLANKSEQ", flankseq);
                free(flankseq);
            }
        }

        if (bcf_get_info_flag(odr->hdr, g->v, "LARGE_REPEAT_REGION", NULL, 0)>0)
        {
            bcf_update_info_flag(odw->hdr, v, "LARGE_REPEAT_REGION", NULL, 0);
        }

        //CG
        std::map<float,int32_t> count_histogram;
        for (uint32_t i=0; i<g->counts.size();++i)
        {
            ++count_histogram[g->counts[i]];
        }

        float cgs[2];
        if (count_histogram.size()==1)
        {
            cgs[0] = g->counts[0];
            cgs[1] = g->counts[0];
        }
        else if (count_histogram.size()>1)
        {
            float first_ca = 0;
            int32_t first_cnt = 0;
            float second_ca = 0;
            int32_t second_cnt = 0;
            std::map<float,int32_t>::iterator i;
            for (i=count_histogram.begin(); i!=count_histogram.end();++i)
            {
                float ca = i->first;
                int32_t cnt = i->second;

                if (cnt>first_cnt)
                {
                    second_ca = first_ca;
                    second_cnt = first_cnt;
                    first_ca = ca;
                    first_cnt = cnt;
                }
                else if (cnt>second_cnt)
                {
                    second_ca = ca;
                    second_cnt = cnt;
                }
            }

            cgs[0] = first_ca;
            cgs[1] = second_ca;
        }
        else
        {
            cgs[0] = -1;
            cgs[1] = -1;
        }
        bcf_update_format_float(odw->hdr, v, "CG", &cgs[0], 2);

        //strand
        char* str = const_cast<char*>(g->strands.c_str());
        bcf_update_format_string(odw->hdr, v, "ST", const_cast<const char**>(&str), 1);

        //alleles
        bcf_update_format_float(odw->hdr, v, "CT", &g->counts[0], g->counts.size());

        //no of mismatches
        bcf_update_format_int32(odw->hdr, v, "NM", &g->no_mismatches[0], g->no_mismatches.size());

        odw->write(v);
        bcf_destroy(v);

        ++no_vntrs_genotyped;
    }
}
