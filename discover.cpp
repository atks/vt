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

#include "discover.h"

namespace
{

typedef struct
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

#define SNP 1
#define MNP 2
#define INDEL 4

/**
 * Class for mining candidate variants.
 *
 * Processes an align read and records the variants observed against the reference
 */
class VariantHunter
{
    public:
    /**
     * Constructor
     * baseq_cutoff - q value cutoff to select candidate SNPs
     */
    VariantHunter(uint32_t vtype,
                  uint32_t evidence_allele_count_cutoff,
                  double fractional_evidence_allele_count_cutoff,
                  uint32_t baseq_cutoff,
                  faidx_t *fai,
                  BCFOrderedWriter *odw)
    : vtype(vtype),
      evidence_allele_count_cutoff(evidence_allele_count_cutoff),
      fractional_evidence_allele_count_cutoff(fractional_evidence_allele_count_cutoff),
      baseq_cutoff(baseq_cutoff),
      fai(fai),
      odw(odw),
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
    void process_read(bam_hdr_t *h, bam1_t *s)
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
    void extract_candidate_variants()
    {
        extract_candidate_variants(chrom, 0, true);
    };

    private:

    uint32_t buffer_size;
    std::vector<std::vector<char> > X; // contains read bases that differ from the genome
    std::vector<std::vector<std::string> > Y; // contains multiple consecutive read bases that differ from the genome
    std::vector<std::vector<std::string> > I; //contains inserted bases
    std::vector<std::vector<std::string> > D; //contains reference bases that are deleted
    std::vector<int32_t> N; // number of evidences observed here - combination of X, I and D
    std::vector<char> REF;
    std::vector<char> ANCHOR;
    std::vector<std::string> ALT;
    char* chrom;

    //key control variables for circular buffer
    uint32_t start, end;
    uint32_t empty_buffer_space;
    uint32_t min_empty_buffer_size;
    uint32_t start_genome_pos0;
    uint32_t max_used_buffer_size_threshold;
    uint32_t max_indel_length;
    uint32_t baseq_cutoff;
    uint32_t evidence_allele_count_cutoff;
    double fractional_evidence_allele_count_cutoff;
    faidx_t *fai;
    uint32_t vtype;
    kstring_t s;
    kstring_t alleles;
    kstring_t read_seq;
    kstring_t qual;
    kstring_t cigar;

    bcf1_t *v;
    BCFOrderedWriter *odw;
    bool debug;

    /**
     * Processes buffer to pick up variants
     * Empty buffer to recover space.
     * @chrom  - remove variants on chrom
     * @pos1   - remove variants up to pos1
     * @flush  - remove all variants
     */
    void extract_candidate_variants(const char* chrom, uint32_t pos1, bool flush=false)
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
                if (vtype&INDEL)
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
                                v = odw->get_bcf1_from_pool();
                                bcf_set_chrom(odw->hdr, v, chrom);
                                bcf_set_pos1(v, start_genome_pos0);
                                alleles.l = 0;
                                kputc(anchor, &alleles);
                                kputc(',', &alleles);
                                kputs(i->first.c_str(), &alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles.s);
                                bcf_update_format_int32(odw->hdr, v, "E", &i->second, 1);
                                bcf_update_format_int32(odw->hdr, v, "N", &N[start], 1);
                                odw->write(v);
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
                                v = odw->get_bcf1_from_pool();
                                bcf_set_chrom(odw->hdr, v, chrom);
                                bcf_set_pos1(v, start_genome_pos0);

                                alleles.l = 0;
                                const char* deletedAllele = i->first.c_str();
                                char replacement_anchor = deletedAllele[0];
                                kputc(anchor, &alleles);
                                ++deletedAllele;
                                kputs(deletedAllele, &alleles);
                                kputc(',', &alleles);
                                kputc(replacement_anchor, &alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles.s);
                                bcf_update_format_int32(odw->hdr, v, "E", &i->second, 1);
                                bcf_update_format_int32(odw->hdr, v, "N", &N[start], 1);
                                odw->write(v);
                            }
                        }
                    }
                }

                if (vtype&SNP)
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
                                v = odw->get_bcf1_from_pool();
                                bcf_set_chrom(odw->hdr, v, chrom);
                                bcf_set_pos1(v, start_genome_pos0+1);
                                alleles.l = 0;
                                kputc(ref, &alleles);
                                kputc(',', &alleles);
                                kputc(i->first, &alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles.s);
                                bcf_update_format_int32(odw->hdr, v, "E", &i->second, 1);
                                bcf_update_format_int32(odw->hdr, v, "N", &N[start], 1);
                                odw->write(v);
                            }
                        }
                    }
                }

                if (vtype&MNP)
                {
                    //handling MNPs
                    mnp_alts.clear();
                    if (Y[start].size()!=0)
                    {
                        for (uint32_t i=0; i<Y[start].size(); ++i)
                        {
                            if (mnp_alts.find(Y[start][i])==mnp_alts.end())
                            {
                                mnp_alts[Y[start][i]] = 1;
                            }
                            else
                            {
                                ++mnp_alts[Y[start][i]];
                            }
                        }

                        for (std::map<std::string, int32_t>::iterator i =mnp_alts.begin(); i!=mnp_alts.end(); ++i)
                        {
                            char* seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), start_genome_pos0, start_genome_pos0+i->first.size()-1, &ref_len);

                            //make sure that we do not output alleles with N bases.
                            if (i->second>= evidence_allele_count_cutoff &&
                                ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff &&
                                !strchr(seq, 'N') && (i->first).find_first_of('N')==std::string::npos)
                            {
                                v = odw->get_bcf1_from_pool();
                                bcf_set_chrom(odw->hdr, v, chrom);
                                bcf_set_pos1(v, start_genome_pos0+1);
                                alleles.l = 0;
                                kputc(anchor, &alleles);
                                kputc(',', &alleles);
                                kputs(i->first.c_str(), &alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles.s);
                                bcf_update_format_int32(odw->hdr, v, "E", &i->second, 1);
                                bcf_update_format_int32(odw->hdr, v, "N", &N[start], 1);
                                odw->write(v);
                            }

                            if (ref_len>0) free(seq);
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

    /**
     * Checks if buffer is empty
     */
    bool is_empty()
    {
        return start==end;
    };

    /**
     *Increments buffer index i by 1.
     */
    void add(uint32_t& i)
    {
        if (i>=buffer_size)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
            exit(1);
        }

        uint32_t temp = (i+1)%buffer_size;
        i = end==i ? (end=temp) : temp;
    };

    /**
     * Increments buffer index i by j.
     */
    uint32_t add(uint32_t i, uint32_t j)
    {
        if (i>=buffer_size)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
            exit(1);
        }

        return (i+j)%buffer_size;
    };

    /**
     * Decrements buffer index i by j.
     */
    uint32_t minus(uint32_t& i, uint32_t j)
    {
        if (i>=buffer_size)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
            exit(1);
        }

        return (i>=j ? i-j : buffer_size-(j-i));
    };

    /**
     * Decrements buffer index i by 1.
     */
    void minus(uint32_t& i)
    {
        if (i>=buffer_size)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
            exit(1);
        }

        i = (i>=1 ? i-1 : buffer_size-1);
    };

    /**
     * Returns the difference between 2 buffer positions
     */
    uint32_t diff(uint32_t i, uint32_t j)
    {
        return (i>=j ? i-j : buffer_size-(j-i));
    };

    /**
     * Gets the position in the buffer that corresponds to
     * the genome position indicated by pos.
     */
    uint32_t get_cur_pos0(uint32_t genome_pos0)
    {
        //when buffer is empty
        if (is_empty())
        {
            start_genome_pos0 = genome_pos0;
            return start;
        }
        else
        {
            if (genome_pos0-start_genome_pos0>buffer_size)
            {
                std::cerr << "overflow buffer\n" ;
                //should allow for unbuffering here

            }
            return (start + (genome_pos0-start_genome_pos0))%buffer_size;
        }
    };

    /**
     * Print buffer contents for debugging purpose
     */
    void printBuffer()
    {
        std::cout << "PRINT BUFFER" << "\n";
        std::cout << "usedBufferSize: " << diff(end,start) << "\n";
        uint32_t cur_pos0 = start;
        uint32_t genome_pos = start_genome_pos0;

        while (cur_pos0!=end)
        {
            std::cout << genome_pos << "\t" << cur_pos0 << "\t" << REF[cur_pos0] << "\t";

            for (uint32_t j=0; j<I[cur_pos0].size(); ++j)
            {
                std::cout << I[cur_pos0][j] << ",";
            }

            for (uint32_t j=0; j<D[cur_pos0].size(); ++j)
            {
                std::cout << D[cur_pos0][j] << ",";
            }

            std::cout << "\t" <<  N[cur_pos0] << "\n";

            add(cur_pos0);
            ++genome_pos;
        }
    };
};

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::vector<GenomeInterval> intervals;
    std::string output_vcf_file;
    std::string input_bam_file;
    std::string ref_fasta_file;
    std::string sample_id;
    uint32_t mapq_cutoff;
    uint32_t baseq_cutoff;
    //takes on snps, mnps, indels
    std::string variant_type;
    uint32_t evidence_allele_count_cutoff;
    double fractional_evidence_allele_count_cutoff;

    uint16_t exclude_flag;

    ///////
    //i/o//
    ///////
    BAMOrderedReader *odr;
    bam1_t *s;

    BCFOrderedWriter *odw;
    bcf1_t *v;

    /////////
    //stats//
    /////////
    uint32_t no_reads;
    uint32_t no_overlapping_reads;
    uint32_t no_passed_reads;
    uint32_t no_exclude_flag_reads;
    uint32_t no_low_mapq_reads;

    /////////
    //tools//
    /////////
    VariantHunter *variantHunter;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Discovers variants from reads in a BAM file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_sample_id("s", "s", "sample ID", true, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_mapq_cutoff("m", "m", "MAPQ cutoff for alignments [20]", false, 20, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_baseq_cutoff("q", "q", "base quality cutoff for bases [13]", false, 13, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_evidence_allele_count_cutoff("e", "e", "evidence count cutoff for candidate allele [2]", false, 2, "int", cmd);
            TCLAP::ValueArg<double> arg_fractional_evidence_allele_count_cutoff("f", "f", "fractional evidence cutoff for candidate allele [0.1]", false, 0.1, "float", cmd);
            TCLAP::ValueArg<std::string> arg_variant_type("v", "v", "variant types [snps,mnps,indels]", false, "snps,mnps,indels", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_bam_file("b", "b", "input BAM file", true, "", "string", cmd);

            cmd.parse(argc, argv);

            input_bam_file = arg_input_bam_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            sample_id = arg_sample_id.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
            mapq_cutoff = arg_mapq_cutoff.getValue();
            baseq_cutoff = arg_baseq_cutoff.getValue();
            variant_type = arg_variant_type.getValue();
            evidence_allele_count_cutoff = arg_evidence_allele_count_cutoff.getValue();
            fractional_evidence_allele_count_cutoff = arg_fractional_evidence_allele_count_cutoff.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////
        exclude_flag = 0x0704;

        odr = new BAMOrderedReader(input_bam_file, intervals);
        s = bam_init1();

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bam_hdr_transfer_contigs_to_bcf_hdr(odr->hdr, odw->hdr);
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=E,Number=1,Type=Integer,Description=\"Number of reads containing evidence of the alternate allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Total number of reads at a candidate locus with reads that contain evidence of the alternate allele\">");
        bcf_hdr_add_sample(odw->hdr, sample_id.c_str());
        bcf_hdr_add_sample(odw->hdr, NULL);
        //bcf_hdr_set_n_sample(odw->hdr, 1);
                                
        v = NULL;

        std::vector<std::string> variant_types;
        split(variant_types, ",", variant_type);
        uint32_t vtype = 0;
        for (uint32_t i = 0; i<variant_types.size(); ++i)
        {
            if (variant_types[i] == "snps")
            {
                vtype |= SNP;
            }
            else if (variant_types[i] == "mnps")
            {
                vtype |= MNP;
            }
            else if (variant_types[i] == "indels")
            {
                vtype |= INDEL;
            }
            else if (variant_types[i] == "all")
            {
                vtype = SNP|MNP|INDEL;
            }
        }

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_reads = 0;
        no_overlapping_reads = 0;
        no_passed_reads = 0;
        no_exclude_flag_reads = 0;
        no_low_mapq_reads = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        faidx_t *fai = fai_load(ref_fasta_file.c_str());
        if (fai==NULL) 
        {
            fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
            exit(1);
        }
        variantHunter = new VariantHunter(vtype,
                                    evidence_allele_count_cutoff,
                                    fractional_evidence_allele_count_cutoff,
                                    baseq_cutoff,
                                    fai,
                                    odw);
    }

    void discover()
    {
        odw->write_hdr();

        //for tracking overlapping reads
        khash_t(rdict) *reads = kh_init(rdict);
        khiter_t k;
        int32_t ret;

        while (odr->read(s))
        {   
            ++no_reads;

            //this read is the first of the pair
            if (bam_get_mpos1(s) && (bam_get_tid(s)==bam_get_mtid(s)))
            {
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
                    //todo: perform stitching in future
                    if((k = kh_get(rdict, reads, bam_get_qname(s)))!=kh_end(reads))
                    {
                        if (kh_exist(reads, k))
                        {
                            free((char*)kh_key(reads, k));
                            kh_del(rdict, reads, k);
                            ++no_overlapping_reads;
                        }
                        //continue;
                    }
                }
            }

            if(bam_get_flag(s) & exclude_flag)
            {
                //1. unmapped
                //2. secondary alignment
                //3. not passing QC
                //4. PCR or optical duplicate
                ++no_exclude_flag_reads;
                continue;
            }

            if (bam_get_mapq(s) < mapq_cutoff)
            {
                //filter short aligments and those with too many indels (?)
                ++no_low_mapq_reads;
                continue;
            }

            if (0)
            {
               bam_print(s);
            }

            variantHunter->process_read(odr->hdr, s);

//          if (no_reads%100000==0) std::cerr << no_reads << "\n";

            ++no_passed_reads;
        }

        odw->close();

    };

    void bam_print(bam1_t *s)
    {
        const char* chrom = bam_get_chrom(odr->hdr, s);
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

        std::cerr << "##################" << "\n";
        std::cerr << "read no  : " << no_reads << "\n";
        std::cerr << "chrom-pos: " << chrom << "-" << pos1 << "\n";
        std::cerr << "read     : " << seq.s << "\n";
        std::cerr << "qual     : " << qual.s << "\n";
        std::cerr << "cigar_str: " << cigar_string.s << "\n";
        std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
        std::cerr << "len      : " << len << "\n";
        std::cerr << "mapq     : " << mapq << "\n";
        std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
        std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";

        if (seq.m) free(seq.s);
        if (qual.m) free(qual.s);
        if (cigar_string.m) free(cigar_string.s);
        if (cigar_expanded_string.m) free(cigar_expanded_string.s);
    }

    void print_options()
    {
        std::clog << "discover v" << version << "\n\n";

        std::clog << "options: [b] input BAM File               " << input_bam_file << "\n";
        std::clog << "         [o] output VCF File              " << output_vcf_file << "\n";
        std::clog << "         [s] sample ID                    " << sample_id << "\n";
        std::clog << "         [r] reference FASTA File         " << ref_fasta_file << "\n";
        std::clog << "         [m] MAPQ cutoff                  " << mapq_cutoff << "\n";
        std::clog << "         [q] base quality cutoff          " << baseq_cutoff << "\n";
        std::clog << "         [v] variant type(s)              " << variant_type << "\n";
        std::clog << "         [e] evidence cutoff              " << evidence_allele_count_cutoff << "\n";
        std::clog << "         [f] fractional evidence cutoff   " << fractional_evidence_allele_count_cutoff<< "\n";
        print_int_op("         [i] intervals                    ", intervals);
        std::clog << "\n";

    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. reads              : " << no_reads << "\n";
        std::clog << "       no. overlapping reads  : " << no_overlapping_reads << "\n";
        std::clog << "       no. low mapq reads     : " << no_low_mapq_reads << "\n";
        std::clog << "       no. passed reads       : " << no_passed_reads << "\n";
        std::clog << "       no. exclude flag reads : " << no_exclude_flag_reads << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void discover(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.discover();
    igor.print_stats();
};
