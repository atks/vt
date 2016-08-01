/* The MIT License

   Copyright (c) 2016 Adrian Tan <atks@umich.edu>

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

#include "read_filter.h"

/**
 * Constructor.
 */
ReadFilter::ReadFilter(uint32_t read_mapq_cutoff, uint16_t read_exclude_flag, bool ignore_overlapping_read)
{
    this->read_mapq_cutoff = read_mapq_cutoff;
    this->read_exclude_flag = read_exclude_flag;
    this->ignore_overlapping_read = ignore_overlapping_read;    
    
    //for tracking overlapping reads
    reads = kh_init(rdict);

    chrom = "";
    tid = -1;
    
    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_reads = 0;
    no_passed_reads = 0;
    no_overlapping_reads = 0;
    no_passed_reads = 0;
    no_exclude_flag_reads = 0;
    no_low_mapq_reads = 0;
    no_unaligned_cigars = 0;
    no_malformed_del_cigars = 0;
    no_malformed_ins_cigars = 0;
    no_salvageable_ins_cigars = 0;
}

/**
 * Filter reads.
 *
 * Returns true if read is failed.
 */
bool ReadFilter::filter_read(bam_hdr_t* h, bam1_t *s)
{
    if (ignore_overlapping_read)
    {
        //check to see that hash should be cleared when encountering new contig.
        //some bams may not be properly formed and contain orphaned reads that
        //is retained in the hash
        if (bam_get_tid(s)!=tid)
        {
            clear_reads();
            tid = bam_get_tid(s);
        }
        
        //this read is part of a mate pair on the same contig
        if (bam_get_mpos1(s) && (bam_get_tid(s)==bam_get_mtid(s)))
        {
            khiter_t k;
            int32_t ret;
    
            //first mate
            if (bam_get_mpos1(s)>bam_get_pos1(s))
            {
                //overlapping
                if (bam_get_mpos1(s)<=(bam_get_pos1(s) + bam_get_l_qseq(s) - 1))
                {
                    //add read that has overlapping
                    //duplicate the record and perform the stitching later
                    char* qname = strdup(bam_get_qname(s));
                    k = kh_put(rdict, reads, qname, &ret);
                    if (!ret)
                    {
                        //already present
                        free(qname);
                    }
                    kh_val(reads, k) = {bam_get_pos1(s), bam_get_pos1(s)+bam_get_l_qseq(s)-1};
                }
            }
            else
            {
                //check overlap
                if((k = kh_get(rdict, reads, bam_get_qname(s)))!=kh_end(reads))
                {
                    if (kh_exist(reads, k))
                    {
                        free((char*)kh_key(reads, k));
                        kh_del(rdict, reads, k);
                        ++no_overlapping_reads;
                    }
                    //set this on to remove overlapping reads.
                    return false;
                }
            }
        }
    }

    if(bam_get_flag(s) & read_exclude_flag)
    {
        //1. unmapped
        //2. secondary alignment
        //3. not passing QC
        //4. PCR or optical duplicate
        ++no_exclude_flag_reads;
        return false;
    }

    if (bam_get_mapq(s) < read_mapq_cutoff)
    {
        //filter short aligments and those with too many indels (?)
        ++no_low_mapq_reads;
        return false;
    }

    //*****************************************************************
    //should we have an assertion on the correctness of the bam record?
    //Is, Ds not sandwiched in M
    //leading and trailing Is - convert to S
    //no Ms!!!!!
    //*****************************************************************
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(s);
        bool seenM = false;
        int32_t last_opchr = '^';

        for (int32_t i = 0; i < n_cigar_op; ++i)
        {
            int32_t opchr = bam_cigar_opchr(cigar[i]);
            int32_t oplen = bam_cigar_oplen(cigar[i]);
            if (opchr=='S')
            {
                if (i!=0 && i!=n_cigar_op-1)
                {
                    std::cerr << "S issue\n";
                    bam_print_key_values(h, s);
                    //++malformed_cigar;
                }
            }
            else if (opchr=='M')
            {
              seenM = true;
            }
            else if (opchr=='D')
            {
                if (last_opchr!='M' || (i<=n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M'))
                {
                    std::cerr << "D issue\n";
                    ++no_malformed_del_cigars;
                    bam_print_key_values(h, s);
                }
            }
            else if (opchr=='I')
            {
                if (last_opchr!='M' || (i<n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M'))
                {
                    if (last_opchr!='M')
                    {
                        if (last_opchr!='^' && last_opchr!='S')
                        {
                            std::cerr << "leading I issue\n";
                            bam_print_key_values(h, s);
                            ++no_malformed_ins_cigars;
                        }
                        else
                        {
                            ++no_salvageable_ins_cigars;
                        }
                    }
                    else if (i==n_cigar_op-1)
                    {
                        ++no_salvageable_ins_cigars;
                    }
                    else if (i==n_cigar_op-2 && (bam_cigar_opchr(cigar[i+1])=='S'))
                    {
                        ++no_salvageable_ins_cigars;
                    }
                    else
                    {
                        std::cerr << "trailing I issue\n";
                        bam_print_key_values(h, s);
                        ++no_malformed_ins_cigars;
                    }
                }
            }

            last_opchr = opchr;
        }

        if (!seenM)
        {
            std::cerr << "NO! M issue\n";
            bam_print_key_values(h, s);
            ++no_unaligned_cigars;
        }
    }

    return true;
}

/**
 * Clear reads from hash.
 */
void ReadFilter::clear_reads()
{
    for (khiter_t k = kh_begin(reads); k != kh_end(reads); ++k)
    {
        if (kh_exist(reads, k))
        {
            free((char*)kh_key(reads, k));
            kh_del(rdict, reads, k);
        }
    }
}

/**
 * Print BAM for debugging purposes.
 */
void ReadFilter::bam_print_key_values(bam_hdr_t *h, bam1_t *s)
{
    const char* chrom = bam_get_chrom(h, s);
    uint32_t pos1 = bam_get_pos1(s);
    kstring_t seq = {0,0,0};
    bam_get_seq_string(s, &seq);
    uint32_t len = bam_get_l_qseq(s);
    kstring_t qual = {0,0,0};
    bam_get_qual_string(s, &qual);
    kstring_t cigar_string = {0,0,0};
    bam_get_cigar_string(s, &cigar_string);
    kstring_t cigar_expanded_string = {0,0,0};
    bam_get_cigar_expanded_string(s, &cigar_expanded_string);
    uint16_t flag = bam_get_flag(s);
    uint32_t mapq = bam_get_mapq(s);

    uint8_t *aux;
    char* md = NULL;
    (aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(aux));

    std::cerr << "##################" << "\n";
    std::cerr << "chrom:pos: " << chrom << ":" << pos1 << "\n";
    std::cerr << "read     : " << seq.s << "\n";
    std::cerr << "qual     : " << qual.s << "\n";
    std::cerr << "cigar_str: " << cigar_string.s << "\n";
    std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
    std::cerr << "len      : " << len << "\n";
    std::cerr << "mapq     : " << mapq << "\n";
    std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
    std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";
    std::cerr << "md       : " << (aux?md:"") << "\n";
    std::cerr << "##################" << "\n";

    if (seq.m) free(seq.s);
    if (qual.m) free(qual.s);
    if (cigar_string.m) free(cigar_string.s);
    if (cigar_expanded_string.m) free(cigar_expanded_string.s);
}