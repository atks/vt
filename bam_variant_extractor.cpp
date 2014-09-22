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

#include "bam_variant_extractor.h"



/**
 * Constructor
 * baseq_cutoff - q value cutoff to select candidate SNPs
 */
BAMVariantExtractor::BAMVariantExtractor(uint32_t vtype,
              uint32_t evidence_allele_count_cutoff,
              double fractional_evidence_allele_count_cutoff,
              uint32_t baseq_cutoff,
              faidx_t *fai)
: vtype(vtype),
  evidence_allele_count_cutoff(evidence_allele_count_cutoff),
  fractional_evidence_allele_count_cutoff(fractional_evidence_allele_count_cutoff),
  baseq_cutoff(baseq_cutoff),
  fai(fai),
  buffer_size(800),
  X(buffer_size),
  Y(buffer_size),
  I(buffer_size),
  D(buffer_size),
  N(buffer_size,0),
  REF(buffer_size),
  ANCHOR(buffer_size),
  chrom(0),
  start(0), end(0),
  empty_buffer_space(0),
  min_empty_buffer_size(400),
  start_genome_pos0(0),
  max_used_buffer_size_threshold(buffer_size-min_empty_buffer_size),
  max_indel_length(50),
  debug(false)
{
    alleles = {0,0,0};
    read_seq = {0,0,0};
    qual = {0,0,0};
    cigar = {0,0,0};
};

/**
 * Transfer read into a buffer for processing later
 */
void BAMVariantExtractor::process_read(bam_hdr_t *h, bam1_t *s)
{
    //extract relevant information from sam record
    const char* chrom = bam_get_chrom(h, s);
    int32_t pos0 = bam_get_pos0(s);
    bam_get_seq_string(s, &read_seq);
    bam_get_cigar_expanded_string(s, &cigar);
    bam_get_qual_string(s, &qual);
    int32_t ref_len;
    uint32_t q;

    //reset buffer if necessary
    if (this->chrom)
    {
        //change in chromosome
        if (strcmp(this->chrom, chrom))
        {
            extract_candidate_variants(this->chrom, UINT_MAX);
            free(this->chrom);
            this->chrom = strdup(chrom);
        }
        else
        {
            extract_candidate_variants(this->chrom, pos0);
        }
    }
    else //first read
    {
        this->chrom = strdup(chrom);
    }

    //basically equivalent to emptying the buffer
    extract_candidate_variants(chrom, pos0);
    char* genome_seq = faidx_fetch_seq(fai, chrom, pos0-1, pos0+cigar.l, &ref_len);
    uint32_t genome_seq_pos0 = 1;
    uint32_t read_seq_pos0 = 0;
    uint32_t cur_pos0 = get_cur_pos0(pos0); //current buffer index
    bool last_position_had_snp = false;
    bool mnp_allele_construction_in_progress = false;
    bool ins_init = true;
    bool del_init = true;
    uint32_t last_snp_pos = 0;
    uint32_t mnp_init_pos = 0;
    char mnp_init_base = 'N';
    uint32_t ins_init_pos = 0;
    uint32_t del_init_pos = 0;
    if (0)
    {
        std::cerr << "===============\n";
        std::cerr << "ADD READ\n";
        std::cerr << "pos1                : " << (pos0+1) << "\n";
        std::cerr << "start_genome_pos0   : " << start_genome_pos0 << "\n";
        std::cerr << "end_genome_pos0     : " << start_genome_pos0 + diff(end,start)  << "\n";
        std::cerr << "start               : " << start << "\n";
        std::cerr << "end                 : " << end << "\n";
        std::cerr << "cur_pos0            : " << cur_pos0 << "\n";
        std::cerr << "chrom               : " << chrom << "\n";
        std::cerr << "genome sequence     : " << genome_seq << "\n";
        std::cerr << "read sequence       : " << read_seq.s << "\n";
        std::cerr << "qual                : " << qual.s << "\n";
        std::cerr << "cigar               : " << cigar.s << "\n";
        std::cerr << "cigar length        : " << cigar.l << "\n";
    }

    if (is_empty())
        start_genome_pos0 = pos0;

    //cycle through the cigar string, this cigar string has essentially been expanded
    for (uint32_t cigar_pos0=0; cigar_pos0<cigar.l; ++cigar_pos0)
    {
        char a = cigar.s[cigar_pos0];

        if (a=='M')
        {
            REF[cur_pos0] = genome_seq[genome_seq_pos0];
            q = qual.s[read_seq_pos0]-33;

            if (genome_seq[genome_seq_pos0]!=read_seq.s[read_seq_pos0] && q>=baseq_cutoff)
            {
                X[cur_pos0].push_back(read_seq.s[read_seq_pos0]);

                //initialize mnp
                if (last_position_had_snp && !mnp_allele_construction_in_progress)
                {
                    mnp_allele_construction_in_progress  = true;
                    mnp_init_pos = last_snp_pos;
                    last_position_had_snp = false;

                    Y[mnp_init_pos].push_back("");
                    Y[mnp_init_pos].back().append(1, mnp_init_base);
                }

                if (mnp_allele_construction_in_progress)
                {
                   Y[mnp_init_pos].back().append(1, read_seq.s[read_seq_pos0]);
                }

                last_position_had_snp = true;
                last_snp_pos = cur_pos0;
                mnp_init_base = read_seq.s[read_seq_pos0];
            }
            else
            {
                last_position_had_snp = false;
                mnp_allele_construction_in_progress = false;
            }

            ins_init = true;
            del_init = true;

            ++N[cur_pos0];
            add(cur_pos0);
            ++genome_seq_pos0;
            ++read_seq_pos0;
        }
        else if (a=='I')
        {
            if (ins_init)
            {
                ins_init_pos = cur_pos0;
                I[ins_init_pos].push_back("");
                I[ins_init_pos].back().append(1, (read_seq_pos0!=0?read_seq.s[read_seq_pos0-1]:genome_seq[genome_seq_pos0-1]));
                ANCHOR[ins_init_pos] = genome_seq[genome_seq_pos0-1];
            }

            last_position_had_snp =false;
            mnp_allele_construction_in_progress = false;
            I[ins_init_pos].back().append(1, read_seq.s[read_seq_pos0]);
            ins_init = false;
            del_init = true;

            //helps maintain count as I's can be 3' hanging
            if (cigar_pos0==cigar.l-1 || cigar.s[cigar_pos0+1]=='S')
            {
                ++N[ins_init_pos];
            }
            ++read_seq_pos0;
        }
        else if (a=='D')
        {
            REF[cur_pos0] = genome_seq[genome_seq_pos0];
            if (del_init)
            {
                del_init_pos = cur_pos0;
                D[del_init_pos].push_back("");
                D[del_init_pos].back().append(1, (read_seq_pos0!=0?read_seq.s[read_seq_pos0-1]:genome_seq[genome_seq_pos0-1]));
                ++N[del_init_pos];
                ANCHOR[del_init_pos] = genome_seq[genome_seq_pos0-1];
            }
            D[del_init_pos].back().append(1, genome_seq[genome_seq_pos0]);

            last_position_had_snp =false;
            mnp_allele_construction_in_progress = false;
            del_init = false;
            ins_init = true;
            add(cur_pos0);
            ++genome_seq_pos0;
        }
        else //S, H and others
        {
            ++read_seq_pos0;
        }
    }

    if (ref_len>0) free(genome_seq);

    if (0)
    {
        std::cout << "final cur_pos0        : " << cur_pos0 << "\n";
        std::cout << "final start_genome_pos: " << start_genome_pos0 << "\n";
        std::cout << "final endGenomePos  : " << start_genome_pos0 + diff(end,start)  << "\n";
        std::cout << "final start         : " << start << "\n";
        std::cout << "final end           : " << end << "\n";
        std::cout << "====================\n";
        printBuffer();
    }
};

/**
 * Processes buffer to pick up variants
 */
void BAMVariantExtractor::extract_candidate_variants()
{
    extract_candidate_variants(chrom, 0, true);
};

/**
 * Processes buffer to pick up variants
 * Empty buffer to recover space.
 * @chrom  - remove variants on chrom
 * @pos1   - remove variants up to pos1
 * @flush  - remove all variants
 */
void BAMVariantExtractor::extract_candidate_variants(const char* chrom, uint32_t pos1, bool flush=false)
{
    //variable to tell when to stop flushing
    uint32_t stop = 0;
    if (flush)
    {
        stop = end;
    }
    else if (is_empty())
    {
        return;
    }
    //extract when separated
    else if (start_genome_pos0+(diff(end,start))<pos1)
    {
        if (debug)
        {
            std::cout << "***********************\n";
            std::cout << "flush buffer segregated\n";

        }
        //flush buffer completely
        stop = end;
    }
    //extract when overlapping
    else
    {
        if (pos1-start_genome_pos0>max_used_buffer_size_threshold)
        {
            stop = add(start, pos1-start_genome_pos0);

            if (debug)
            {
                std::cout << "************************\n";
                std::cout << "flush buffer overlapping\n";
            }

        }
        else
        {
            return;
        }
    }

//  if (debug)
//  {
//      std::cout << "pos1               : " << pos1 << "\n";
//      std::cout << "stop              : " << stop << "\n";
//      std::cout << "usedBufferSize    : " << diff(end,start) << "\n";
//      std::cout << "buffer_size        : " << buffer_size << "\n";
//      std::cout << "min_empty_buffer_size: " << min_empty_buffer_size << "\n";
//      std::cout << "start_genome_pos0    : " << start_genome_pos0 << "\n";
//      std::cout << "endGenomePos      : " << start_genome_pos0 + diff(end,start)  << "\n";
//      std::cout << "start             : " << start << "\n";
//      std::cout << "end               : " << end << "\n";
//
//      std::cout << pos1 << "," << diff(end,start) << "," <<  buffer_size-min_empty_buffer_size << "\n";
//  }

    std::map<char, int32_t> snp_alts;
    std::map<std::string, int32_t> mnp_alts;
    std::map<std::string, int32_t> indel_alts;
    char anchor, ref;
    int32_t ref_len;
    //print out candidate variants
    while (start!=stop)
    {
        //assayed position
        if (N[start]>=1)
        {
            if (vtype&VT_INDEL)
            {
                //handling insertions
                indel_alts.clear();
                anchor = ANCHOR[start];

                if (I[start].size()!=0)
                {
                    for (uint32_t i=0; i<I[start].size(); ++i)
                    {
                        if (indel_alts.find(I[start][i])==indel_alts.end())
                        {
                            indel_alts[I[start][i]] = 1;
                        }
                        else
                        {
                            ++indel_alts[I[start][i]];
                        }
                    }

                    for (std::map<std::string, int32_t>::iterator i =indel_alts.begin(); i!=indel_alts.end(); ++i)
                    {
                        //make sure that we do not output alleles with N bases.
                        if (i->second>= evidence_allele_count_cutoff &&
                            ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff &&
                            anchor!='N' && (i->first).find_first_of('N')==std::string::npos)
                        {
                            alleles.l = 0;
                            kputc(anchor, &alleles);
                            kputc(',', &alleles);
                            kputs(i->first.c_str(), &alleles);
//                                bcf_update_alleles_str(odw->hdr, v, alleles.s);
//                                bcf_update_format_int32(odw->hdr, v, "E", &i->second, 1);
//                                bcf_update_format_int32(odw->hdr, v, "N", &N[start], 1);
//                                odw->write(v);
                        }
                    }
                }

                //handling deletions
                indel_alts.clear();
                if (D[start].size()!=0)
                {
                    for (uint32_t i=0; i<D[start].size(); ++i)
                    {
                        if (indel_alts.find(D[start][i])==indel_alts.end())
                        {
                            indel_alts[D[start][i]] = 1;
                        }
                        else
                        {
                            ++indel_alts[D[start][i]];
                        }
                    }

                    for (std::map<std::string, int32_t>::iterator i = indel_alts.begin(); i!= indel_alts.end(); ++i)
                    {
                        //make sure that we do not output alleles with N bases.
                        if (i->second>= evidence_allele_count_cutoff &&
                            ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff &&
                            anchor!='N' && (i->first).find_first_of('N')==std::string::npos)
                        {
                            alleles.l = 0;
                            const char* deletedAllele = i->first.c_str();
                            char replacement_anchor = deletedAllele[0];
                            kputc(anchor, &alleles);
                            ++deletedAllele;
                            kputs(deletedAllele, &alleles);
                            kputc(',', &alleles);
                            kputc(replacement_anchor, &alleles);
                        }
                    }
                }
            }

            if (vtype&VT_SNP)
            {
                //handling SNPs
                snp_alts.clear();

                ref = REF[start];

                if (X[start].size()!=0)
                {
                    for (uint32_t i=0; i<X[start].size(); ++i)
                    {
                        if (snp_alts.find(X[start][i])==snp_alts.end())
                        {
                            snp_alts[X[start][i]] = 1;
                        }
                        else
                        {
                            ++snp_alts[X[start][i]];
                        }
                    }

                    for (std::map<char, int32_t>::iterator i =snp_alts.begin(); i!=snp_alts.end(); ++i)
                    {
                        //make sure that we do not output alleles with N bases.
                        if (i->second>= evidence_allele_count_cutoff &&
                            ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff &&
                            ref!='N' && (i->first)!='N')
                        {
                            alleles.l = 0;
                            kputc(ref, &alleles);
                            kputc(',', &alleles);
                            kputc(i->first, &alleles);
                        }
                    }
                }
            }
        }

        Y[start].clear();
        X[start].clear();
        I[start].clear();
        D[start].clear();
        N[start] = 0;

        add(start);
        ++start_genome_pos0;
    }

//          //clean up final position too
//          if (is_empty())
//          {
//              X[start].clear();
//              I[start].clear();
//              D[start].clear();
//              N[start] = 0;
//          }

//  if (debug)
//  {
//      std::cout << "final start : " << start << "\n";
//      std::cout << "final end   : " << end << "\n";
//      std::cout << "*************\n";
//      printBuffer();
//  }
};
