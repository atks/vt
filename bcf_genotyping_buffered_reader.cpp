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
        if (bcf_get_n_allele(g->v)==2)
        {
            char strand = bam_is_rev(s) ? 'R' : 'F';
            char allele = 'R';
            uint32_t qual = 30;
            uint32_t cycle = 10;
            uint8_t *nm_aux;
            int32_t no_mismatches = 0;
            ((nm_aux=bam_aux_get(s, "NM")) &&  (no_mismatches = bam_aux2i(nm_aux)));

            if (false)
            {
                int32_t n_cigar_op = bam_get_n_cigar_op(s);
                int32_t no_mismatches1 = 0;
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
    
        //                        std::cerr << bcf_get_allele(g->v)[0][0]  << " vs " << bam_base2char(bam_seqi(bseq, rpos1)) << "\n";
    
                                char obs_allele = bam_base2char(bam_seqi(bseq, rpos1));
                                
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
                                qual = bqual[rpos1];
                               // std::cerr << "read length " << bam_get_l_qseq(s) << " " << rpos1 << "\n";
                                cycle = strand == 'F' ? rpos1 : (bam_get_l_qseq(s) - rpos1);
    
                             //   break;
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
                            ++no_mismatches1;
                            rpos1 += oplen;
                        }
                        else if (opchr=='S')
                        {
                            rpos1 += oplen;
                        }
                    }
                }
    
                uint8_t *md_aux;
                char* md = 0;
                ((md_aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(md_aux)));
                char* mdp = md;
                bool indel = false;
                uint8_t* bqual = bam_get_qual(s);
                uint32_t rpos1 = 0;
                std::string digit_string;
                uint32_t no_matches = 0;
                while (*mdp)
                {
                    if (isdigit(*mdp))
                    {
                        str2uint32(digit_string, no_matches);
                        rpos1 += no_matches;
                        indel = false;
                    }
                    else if (*mdp=='N')
                    {
                        ++rpos1;
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
                            ++rpos1;
                        }
                    }
    
                    ++mdp;
                }
            }
            else
            {
//
//            //iterate cigar
//            uint32_t n_cigar_op = bam_get_n_cigar_op(s);
//
//            char* mdp = md;
//            uint32_t cpos1 = pos1; //current 1 based genome position
//            uint32_t spos0 = 0;    //current position in read sequence
//
//            //variables for I's embedded in Matches in the MD tag
//            uint32_t md_mlen_left = 0;

//            if (n_cigar_op)
//            {
//                uint32_t *cigar = bam_get_cigar(s);
//                bool seenM = false;
//
//                if (debug>=3) pileup.print_state();
//                for (uint32_t i = 0; i < n_cigar_op; ++i)
//                {
//                    uint32_t oplen = bam_cigar_oplen(cigar[i]);
//                    char opchar = bam_cigar_opchr(cigar[i]);
//
//                    if (debug) std::cerr << "CIGAR: " << oplen << " " << opchar << "\n";
//
//                    if (opchar=='S')
//                    {
//                        if (i==n_cigar_op-1 || (i==0 && n_cigar_op>=2 && bam_cigar_opchr(cigar[1])=='M'))
//                        {
//                            //add to S evidence
//                            std::string ins = "";
//                            float mean_qual = 0;
//                            for (size_t j=0; j<oplen ; ++j)
//                            {
//                                ins += bam_base2char(bam_seqi(seq, spos0+j));
//                                mean_qual += qual[spos0+j];
//                            }
//                            mean_qual /= oplen;
//
//                            if (cpos1==pos1)
//                            {
//                                if (debug) std::cerr << "\t\t\tadding LSCLIP: " << cpos1 << "\t" << ins << " {" << mean_qual << "}\n";
//                                pileup.add_lsclip(cpos1, ins, mean_qual, strand);
//                            }
//                            else if (seenM)
//                            {
//                                if (debug) std::cerr << "\t\t\tadding RSCLIP: " << (cpos1-1) << "\t" << ins << " {" << mean_qual << "}\n";
//                                pileup.add_rsclip(cpos1-1, ins, mean_qual, strand);
//                            }
//                        }
//
//                        spos0 += oplen;
//                    }
//                    else if (opchar=='M')
//                    {
//                        uint32_t lpos1 = cpos1; // we need this because M contains matches and mismatches
//                        uint32_t sspos0 = spos0; // we need this because M contains matches and mismatches
//                        uint32_t mlen = oplen;
//                        uint32_t i = 0;
//                        seenM = true;
//
//                        if (debug) std::cerr << "\t\tmd len left : " << md_mlen_left << "\n";
//                        if (debug) std::cerr << "\t\tmlen : " << mlen << "\n";
//                        if (debug) std::cerr << "\t\tmdp : " << mdp << "\n";
//
//                        //left over MD matches to handle.
//                        if (md_mlen_left)
//                        {
//                            uint32_t ilen = md_mlen_left<=mlen ? md_mlen_left : mlen;
//                            pileup.add_ref(lpos1, sspos0, ilen, seq);
//
//                            if (debug)
//                            {
//                                uint32_t gbeg1 = lpos1;
//                                uint32_t gend1 = lpos1+ilen-1;
//
//                                std::cerr << "\t\t\tadding REF: " << gbeg1 << "-" << gend1 << ":";
//                                for (size_t i=sspos0; i<=(sspos0+ilen-1); ++i)
//                                {
//                                    std::cerr << (bam_base2char(bam_seqi(seq, i)));
//                                }
//                                std::cerr << " (" << gend1-gbeg1+1 << ") [" <<  mlen-ilen << "]\n";
//                            }
//
//                            lpos1 += ilen;
//                            sspos0 += ilen;
//
//                            if (md_mlen_left>=mlen)
//                            //yet another insertion
//                            {
//                                md_mlen_left -= ilen;
//                                cpos1 += ilen;
//                                spos0 += ilen;
//                                continue;
//                            }
//                            //a snp next
//                            else
//                            {
//                                md_mlen_left = 0;
//                                mlen -= ilen;
//                                //go to loop in the next section
//                            }
//                        }
//
//                        while (*mdp)
//                        {
//                            if (isalpha(*mdp)) //SNPs
//                            {
//                                char ref = toupper(*mdp);
//                                char alt = (bam_base2char(bam_seqi(seq, spos0+(lpos1-cpos1))));
//                                if (debug) std::cerr << "\tMD: Mismatch " << ref << "\n";
//                                if (debug) std::cerr << "\t\t\tadding SNP: " << lpos1 << ":" << ref << "/" << alt << " [" << (mlen-1)<< "]\n";
////                                if (qual[sspos0]>vf.get_snp_baseq_cutoff())
////                                {
//                                    pileup.add_snp(lpos1, ref, alt, qual[sspos0], vf.get_snp_baseq_cutoff());
////                                }
//                                ++lpos1;
//                                ++mdp;
//                                ++sspos0;
//                                --mlen;
//                            }
//                            else if (isdigit(*mdp)) //matches
//                            {
//                                char* end = 0;
//                                int32_t len = std::strtol(mdp, &end, 10);
//                                mdp = end;
//
//                                if (debug) std::cerr << "\tMD: Match " << len << "\n";
//
//                                if (len)
//                                {
//                                    uint32_t ilen = len<=mlen ? len : mlen;
//
//                                    if (debug)
//                                    {
//                                        uint32_t gbeg1 = lpos1;
//                                        uint32_t gend1 = lpos1+ilen-1;
//
//                                        //std::cerr << "\t\t\tadding REF: " << gbeg1 << "-" << gend1 << ":";
//                                        for (size_t i=sspos0; i<=(sspos0+ilen-1); ++i)
//                                        {
//                                            std::cerr << (bam_base2char(bam_seqi(seq, i)));
//                                        }
//                                        std::cerr << " (" << gend1-gbeg1+1 << ") [" <<  mlen-ilen << "]\n";
//                                    }
//
//                                    pileup.add_ref(lpos1, sspos0, ilen, seq);
//
//                                    lpos1 += ilen;
//                                    sspos0 += ilen;
//
//                                    //next up an insertion
//                                    if (len>mlen)
//                                    {
//                                        md_mlen_left = len - mlen;
//                                        break;
//                                    }
//                                    else
//                                    {
//                                        mlen -= ilen;
//                                    }
//                                }
//                            }
//                            else // deletion
//                            {
//                                break;
//                            }
//
//                            if (mlen==0)
//                            {
//                                break;
//                            }
//                        }
//
//                        //note that only insertions, matches and mismatches can only occur here.
//
//                        cpos1 += oplen;
//                        spos0 += oplen;
//                    }
//                    else if (opchar=='D')
//                    {
//                        bool is_del = false;
//
//                        if (*mdp=='0') ++mdp;
//
//                        if (*mdp!='^')
//                        {
//                            bam_print_key_values(odr->hdr, s);
//                            std::cerr << "mdp: " << mdp << "\n";
//                            std::cerr << "inconsistent MD and cigar, deletion does not occur at the right place.\n";
//                            exit(1);
//                        }
//                        else
//                        {
//                            ++mdp;
//                            std::string del = "";
//                            while (isalpha(*mdp))
//                            {
//                                del += toupper(*mdp);
//                                ++mdp;
//                            }
//
//                            if (debug) std::cerr << "\t\t\tadding DEL: " << (cpos1-1) << " " << del << "\n";
//                            pileup.add_del((cpos1-1), del);
//
//                            cpos1 += oplen;
//                        }
//                    }
//                    else if (opchar=='I')
//                    {
//                        //leading Is
//                        if (!seenM)
//                        {
////                            //add to S evidence
////                            std::string ins = "";
////                            float mean_qual = 0;
////                            for (size_t j=0; j<oplen ; ++j)
////                            {
////                                ins += bam_base2char(bam_seqi(seq, spos0+j));
////                                mean_qual += qual[spos0+j];
////                            }
////                            mean_qual /= oplen;
////
////                            if (mean_qual>sclip_mq_cutoff)
////                            {
////                                if (debug) std::cerr << "\t\t\tadding LSCLIP: " << cpos1 << "\t" << ins << " {" << mean_qual << "}\n";
////                                pileup.add_lsclip(cpos1, ins, mean_qual, strand);
////                            }
//
//                            spos0 += oplen;
//                        }
//                        //trailing Is
//                        else if (i==n_cigar_op-1 || (i+2==n_cigar_op && bam_cigar_opchr(cigar[n_cigar_op-1])=='S'))
//                        {
//                            //bam_print_key_values(odr->hdr, s);
//                            spos0 += oplen;
//                        }
//                        else
//                        {
//                            //insertions are not present in MD tags
//                            //may be handled independently of future matches
//                            std::string ins = "";
//                            for (size_t i=0; i<oplen ; ++i)
//                            {
//                                ins += bam_base2char(bam_seqi(seq, spos0+i));
//                            }
//
//                            if (debug) std::cerr << "\t\t\tadding INS: " << (cpos1-1) << " " << ins  << "\n";
//                            pileup.add_ins((cpos1-1), ins, pos1);
//
//                            spos0 += oplen;
//                        }
//                    }
//                    else
//                    {
//                        std::cerr << "never seen before state " << opchar << "\n";
//                    }
//                }
//            }
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

    //                        std::cerr << bcf_get_allele(g->v)[0][0]  << " vs " << bam_base2char(bam_seqi(bseq, rpos1)) << "\n";

                            allele = bam_base2char(bam_seqi(bseq, rpos1)) == bcf_get_allele(g->v)[0][0] ? 'R' : 'A';
                            qual = bqual[rpos1];
                           // std::cerr << "read length " << bam_get_l_qseq(s) << " " << rpos1 << "\n";
                            cycle = strand == 'F' ? rpos1 : (bam_get_l_qseq(s) - rpos1);

                         //   break;
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
                        ++no_mismatches;
                        rpos1 += oplen;
                    }
                    else if (opchr=='S')
                    {
                        rpos1 += oplen;
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
            
            //ADF
            bcf_update_format_int32(odw->hdr, v, "ADF", &g->allele_depth_fwd[0], 2);
            
            //ADR
            bcf_update_format_int32(odw->hdr, v, "ADR", &g->allele_depth_rev[0], 2);
            
            //depth
            bcf_update_format_int32(odw->hdr, v, "DP", &g->depth, 1);

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
            bcf_update_format_int32(odw->hdr, v, "BQSUM", &g->base_qualities_sum, 1);
            
            //depth
            bcf_update_format_int32(odw->hdr, v, "DPF", &g->depth_fwd, 1);
            bcf_update_format_int32(odw->hdr, v, "DPR", &g->depth_rev, 1);
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

        ++no_indels_genotyped;
    }
    else if (g->vtype==VT_VNTR)
    {
        ++no_vntrs_genotyped;
    }
}