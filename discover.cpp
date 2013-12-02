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

typedef struct 
{ 
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

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
      bufferSize(800),
      X(bufferSize),
      Y(bufferSize),
      I(bufferSize), 
      D(bufferSize),
      N(bufferSize,0),
      REF(bufferSize),
      start(0), end(0),
      empty_buffer_space(0),
      min_empty_buffer_size(400),
      start_genome_pos1(0),
      max_used_buffer_size_threshold(bufferSize-min_empty_buffer_size),
      max_indel_length(50),
      debug(false)
    {   
        alleles = {0,0,0};
        read_seq_string = {0,0,0};
        qual_string = {0,0,0};
        cigar_string = {0,0,0};
        v = bcf_init();
    };
    
    /**
     * Transfer read into a buffer for processing later
     */
    void process_read(bam_hdr_t *h, bam1_t *s)
    {   
        //extract relevant information from sam record
        const char* chrom = bam_get_chrom(h, s);
        uint32_t pos1 = bam_get_pos1(s);
        bam_get_seq_string(s, &read_seq_string);
        const char* read_seq = read_seq_string.s;
        bam_get_cigar_string(s, &cigar_string);
        const char* cigar = cigar_string.s;
        bam_get_qual_string(s, &qual_string);
        const char* qual = qual_string.s;
        char a;
        int32_t ref_len;
        
        uint32_t q;
                                
        //basically equivalent to emptying the buffer
        extractCandidateVariants(chrom, pos1);
        
        const char* genome_seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), pos1-1, pos1+cigar_string.l-1, &ref_len);
                    
        //myRefSeq->getString(genome_seqString, chromNo, (genomeIndex_t) pos-1, cigarString.size()+2);
        //const char* genome_seq = genome_seqString.c_str();
        
        uint32_t genome_seq_pos1 = 1;
        uint32_t read_seq_pos = 0;
        uint32_t cur_pos = getcur_pos(pos1); //current buffer index
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
            std::cerr << "start_genome_pos    : " << start_genome_pos1 << "\n";
            std::cerr << "endGenomePos        : " << start_genome_pos1 + diff(end,start)  << "\n"; 
            std::cerr << "start               : " << start << "\n";
            std::cerr << "end                 : " << end << "\n";
            std::cerr << "cur_pos             : " << cur_pos << "\n";
            std::cerr << "chrom               : " << chrom << "\n";
            std::cerr << "genome sequence     : " << genome_seq << "\n";
        }
        
    //          if(pos<239411 && pos+strlen(read_seq)>239411)
    //          {
    //              std::cerr << "=====================\n";
    //              std::cerr << "CIGAR " <<  cigar << "\n";
    //              std::cerr << "POS " <<  pos << "\n";    
    //              std::cerr << "SEQ " <<  read_seq << "\n";    
    //          
    //                  
    //          }
    
    //    if (0 && pos>239400 && pos<239421)
    //    {
    //        std::cerr << "=====================\n";
    //        std::cerr << "POS   " <<  pos << "\n";
    //        std::cerr << "CIGAR " <<  cigar << "\n";
    //        std::cerr << "CIGAR STRING " <<  record.getCigar() << "\n";
    //        std::cerr << "SEQ   " <<  read_seq << "\n";    
    //    }
        
        if (isEmpty())
            start_genome_pos1 = pos1;
        
        //cycle through the cigar string, this cigar string has essentially been expanded
        for (uint32_t cigar_pos=0; cigar_pos<cigar_string.l; ++cigar_pos)
        {
            a = cigar[cigar_pos];
                            
            if (a=='M')
            {
                REF[cur_pos] = genome_seq[genome_seq_pos1];
                //get quality
                q = qual[read_seq_pos]-33;
                
    //          if ((bq[read_seq_pos]-64)>(qual[read_seq_pos]-33))
    //          {   
    //              std::cerr << "===========\n"; 
    //              std::cerr << "baq " << q << "\n";
    //              std::cerr << "q "<< qual[read_seq_pos]-33 << "\n";
    //              std::cerr << "bq "<< bq[read_seq_pos]-64 << "\n";
    //          }
                    
                if (genome_seq[genome_seq_pos1]!=read_seq[read_seq_pos] && q>=baseq_cutoff)
                {
                    X[cur_pos].push_back(read_seq[read_seq_pos]);
                    
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
                       Y[mnp_init_pos].back().append(1, read_seq[read_seq_pos]);
                    }
                    
                    last_position_had_snp = true;
                    last_snp_pos = cur_pos;
                    mnp_init_base = read_seq[read_seq_pos];
                }
                else
                {
                    last_position_had_snp = false;
                    mnp_allele_construction_in_progress = false;
                }
                
                ins_init = true;
                del_init = true;
    
                ++N[cur_pos];                    
                add(cur_pos);
                ++genome_seq_pos1;
                ++read_seq_pos;
            }
            else if (a=='I')
            {
                if (ins_init)
                {   
                    ins_init_pos = minus(cur_pos,1);
                    I[ins_init_pos].push_back("");
                    //place anchor
                    I[ins_init_pos].back().append(1, (read_seq_pos!=0?read_seq[read_seq_pos-1]:genome_seq[genome_seq_pos1-1]));
                                    
                    //if (read_seq_pos==0)
                    //    ++N[cur_pos];
                }
                
                last_position_had_snp =false;
                mnp_allele_construction_in_progress = false;
                I[ins_init_pos].back().append(1, read_seq[read_seq_pos]);
                ins_init = false;
                del_init = true;
                ++read_seq_pos;                   
            }
            else if (a=='D')
            {               
                REF[cur_pos] = genome_seq[genome_seq_pos1];  
                if (del_init)
                {
                    D[cur_pos].push_back("");
                    D[cur_pos].back().append(1, (read_seq_pos!=0?read_seq[read_seq_pos-1]:genome_seq[genome_seq_pos1-1]));
                    ++N[cur_pos];
                    del_init_pos = cur_pos;
                    
                    //if (read_seq_pos==0)
                    //    ++N[cur_pos];
                }
                D[del_init_pos].back().append(1, genome_seq[genome_seq_pos1]);
                
                last_position_had_snp =false;
                mnp_allele_construction_in_progress = false;
                del_init = false;
                ins_init = true;
                add(cur_pos);
                ++genome_seq_pos1;
            }
            else //S, H and others
            {
                ++read_seq_pos;
            }
        }
        
    //  if (debug && (pos>239300 && pos<239500))
    //  {
    //      std::cout << "final cur_pos        : " << cur_pos << "\n";
    //      std::cout << "final start_genome_pos: " << start_genome_pos1 << "\n";
    //      std::cout << "final endGenomePos  : " << start_genome_pos1 + diff(end,start)  << "\n"; 
    //      std::cout << "final start         : " << start << "\n";
    //      std::cout << "final end           : " << end << "\n";
    //      std::cout << "====================\n";  
    //      printBuffer();
    //  }
                
    };

    /**
     * Processes buffer to pick up variants
     */ 
    void extractCandidateVariants()
    {
        extractCandidateVariants(chrom, 0, true);
    };

    /**
     * Reset buffer for new chromosome
     */
    void reset(const char* chrom)
    {
        extractCandidateVariants(chrom, UINT_MAX);
        if (!this->chrom) 
        {
          //  free(this->chrom);
        }
        this->chrom = strdup(chrom);        
    };
        
    private:

    uint32_t bufferSize;
    std::vector<std::vector<char> > X; // contains read bases that differ from the genome
    std::vector<std::vector<std::string> > Y; // contains multiple consecutive read bases that differ from the genome
    std::vector<std::vector<std::string> > I; //contains inserted bases
    std::vector<std::vector<std::string> > D; //contains reference bases that are deleted
    std::vector<int32_t> N; // number of evidences observed here - combination of X, I and D
    std::vector<char> REF;
    std::vector<std::string> ALT;
    const char* chrom;
    uint32_t chromNo;
    uint32_t start, end;
    uint32_t empty_buffer_space;
    uint32_t min_empty_buffer_size;
    uint32_t start_genome_pos1;
    uint32_t max_used_buffer_size_threshold;
    uint32_t max_indel_length;
    uint32_t baseq_cutoff;
    uint32_t evidence_allele_count_cutoff;
    double fractional_evidence_allele_count_cutoff;
    faidx_t *fai;
    uint32_t vtype;
    kstring_t s;
    kstring_t alleles;
    kstring_t read_seq_string;
    kstring_t qual_string;
    kstring_t cigar_string;
    
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
    void extractCandidateVariants(const char* chrom, uint32_t pos1, bool flush=false)
    {
        //variable to tell when to stop flushing
        uint32_t stop = 0;
        if (flush)
        {
            stop = end;
        }
        else if (isEmpty())
        {
            return;             
        }
        //extract when separated
        else if (start_genome_pos1+(diff(end,start))<pos1)
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
            if (pos1-start_genome_pos1>max_used_buffer_size_threshold)
            {
                stop = add(start, pos1-start_genome_pos1);
    
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
            
    //              stop = add(start, pos1)
    //                              
    //              if (diff(end,start)>max_used_buffer_size_threshold)
    //              {
    //                  stop = add(start, (pos1>max_indel_length+start_genome_pos1 ? (pos1-max_indel_length)-start_genome_pos1 : 0));
    //              }
        }
                    
    //  if (pos1==160141 && debug)
    //  {
    //      std::cout << "pos1               : " << pos1 << "\n";
    //      std::cout << "stop              : " << stop << "\n";
    //      std::cout << "usedBufferSize    : " << diff(end,start) << "\n";
    //      std::cout << "bufferSize        : " << bufferSize << "\n";
    //      std::cout << "min_empty_buffer_size: " << min_empty_buffer_size << "\n";
    //      std::cout << "start_genome_pos1    : " << start_genome_pos1 << "\n";
    //      std::cout << "endGenomePos      : " << start_genome_pos1 + diff(end,start)  << "\n";
    //      std::cout << "start             : " << start << "\n";
    //      std::cout << "end               : " << end << "\n";
    //          
    //      std::cout << pos1 << "," << diff(end,start) << "," <<  bufferSize-min_empty_buffer_size << "\n";
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
                if (vtype&4)
                {
                    //handling insertions
                    indel_alts.clear();                 
                    char* seq = faidx_fetch_seq(fai, chrom, start_genome_pos1-1, start_genome_pos1-1, &ref_len);
                    anchor = seq[0];
                    if (ref_len<1) free(seq);
                    
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
                                bcf_set_pos1(v, start_genome_pos1);
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
                        
                        char* seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), start_genome_pos1-1, start_genome_pos1-1, &ref_len);
                        anchor = seq[0];
                        if (ref_len<1) free(seq);
                        //myRefSeq->getString(seq, chromNo, (genomeIndex_t) (start_genome_pos1-1), 1);
                        //anchor = seq.at(0);
                        
                        for (std::map<std::string, int32_t>::iterator i = indel_alts.begin(); i!= indel_alts.end(); ++i)
                        {   
                            //make sure that we do not output alleles with N bases.
                            if (i->second>= evidence_allele_count_cutoff && 
                                ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff && 
                                anchor!='N' && (i->first).find_first_of('N')==std::string::npos)
                            {
                                v = odw->get_bcf1_from_pool();
                                bcf_set_chrom(odw->hdr, v, chrom);
                                bcf_set_pos1(v, start_genome_pos1);
                                alleles.l = 0;
                                kputc(anchor, &alleles);
                                kputc(',', &alleles);
                                const char* deletedAllele = i->first.c_str();
                                ++deletedAllele;
                                kputs(deletedAllele, &alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles.s);                                
                                bcf_update_format_int32(odw->hdr, v, "E", &i->second, 1);
                                bcf_update_format_int32(odw->hdr, v, "N", &N[start], 1);
                                odw->write(v);
                            }
                        }                   
                    }
                }
                
                if (vtype&1)
                {
                    //handling SNPs
                    snp_alts.clear();
                    
                    //myRefSeq->getString(seq, chromNo, (genomeIndex_t) (start_genome_pos1), 1);
                    //ref = seq.at(0);
                    char* seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), start_genome_pos1-1, start_genome_pos1-1, &ref_len);
                    ref = seq[0];
                    if (ref_len<1) free(seq);
                        
                    
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
                                bcf_set_pos1(v, start_genome_pos1);
                                alleles.l = 0;
                                kputc(anchor, &alleles);
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
                
                if (vtype&2)
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
                            char* seq = faidx_fetch_seq(fai, const_cast<char*>(chrom), start_genome_pos1-1, start_genome_pos1+i->first.size()-1, &ref_len);
//                          myRefSeq->getString(seq, chromNo, (genomeIndex_t) (start_genome_pos1), i->first.size());
                            
                            //make sure that we do not output alleles with N bases.
                            if (i->second>= evidence_allele_count_cutoff && 
                                ((double)i->second/(double) N[start]) >= fractional_evidence_allele_count_cutoff && 
                                strchr(seq, 'N') && (i->first).find_first_of('N')==std::string::npos)
                            {
                                                                bcf_set_chrom(odw->hdr, v, chrom);
                                bcf_set_pos1(v, start_genome_pos1);
                                alleles.l = 0;
                                kputc(anchor, &alleles);
                                kputc(',', &alleles);
                                kputs(i->first.c_str(), &alleles);
                                bcf_update_alleles_str(odw->hdr, v, alleles.s);
                                bcf_update_format_int32(odw->hdr, v, "E", &i->second, 1);
                                bcf_update_format_int32(odw->hdr, v, "N", &N[start], 1);
                                odw->write(v);
                            }
                            
                            if (ref_len<1) free(seq);
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
            ++start_genome_pos1;   
        }
        
    //          //clean up final position too
    //          if (isEmpty())
    //          {
    //              X[start].clear();
    //              I[start].clear();
    //              D[start].clear();
    //              N[start] = 0;
    //          }
        
    //  if (pos==160141  && debug)
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
    bool isEmpty()
    {
        return start==end;          
    };
    
    /**
     *Increments buffer index i by 1.
     */
    void add(uint32_t& i)
    {
        if (i>=bufferSize)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
            exit(1);
        }
        
        uint32_t temp = (i+1)%bufferSize;
        i = end==i ? (end=temp) : temp; 
    };
    
    /**
     * Increments buffer index i by j.
     */
    uint32_t add(uint32_t i, uint32_t j)
    {
        if (i>=bufferSize)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
            exit(1);
        }
        
        return (i+j)%bufferSize;
    };

    /**
     * Decrements buffer index i by j.
     */
    uint32_t minus(uint32_t& i, uint32_t j)
    {
        if (i>=bufferSize)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
            exit(1);
        }
        
        return (i>=j ? i-j : bufferSize-(j-i));
    };
    
    /**
     * Decrements buffer index i by 1.
     */
    void minus(uint32_t& i)
    {
        if (i>=bufferSize)
        {
            std::cerr << "Unaccepted buffer index: " << i << " ("  << bufferSize << ")\n";
            exit(1);
        }
        
        i = (i>=1 ? i-1 : bufferSize-1);
    };
    
    /**
     * Returns the difference between 2 buffer positions
     */
    uint32_t diff(uint32_t i, uint32_t j)
    {
        return (i>=j ? i-j : bufferSize-(j-i));
    };
    
    /**
     * Gets the position in the buffer that corresponds to 
     * the genome position indicated by pos.
     */
    uint32_t getcur_pos(uint32_t genomePos)
    {
        //when buffer is empty
        if (isEmpty())
        {
            start_genome_pos1 = genomePos;
            return start;
        }
        else
        {
            if (genomePos-start_genome_pos1>bufferSize)
            {
                std::cerr << "overflow buffer\n" ;
                //should allow for unbuffering here
                
            }
            return (start + (genomePos-start_genome_pos1))%bufferSize;
        }
    };
    
    /**
     * Print buffer contents for debugging purpose
     */
    void printBuffer()
    {
        std::cout << "PRINT BUFFER" << "\n";
        std::cout << "usedBufferSize: " << diff(end,start) << "\n";
        uint32_t cur_pos = start;
        uint32_t genome_pos = start_genome_pos1;
            
        while (cur_pos!=end)
        {
            std::cout << genome_pos << "\t" << cur_pos << "\t" << REF[cur_pos] << "\t";
                
            for (uint32_t j=0; j<I[cur_pos].size(); ++j)
            {
                std::cout << I[cur_pos][j] << ","; 
            }
            
            for (uint32_t j=0; j<D[cur_pos].size(); ++j)
            {
                std::cout << D[cur_pos][j] << ","; 
            }
            
            std::cout << "\t" <<  N[cur_pos] << "\n";
            
            add(cur_pos);
            ++genome_pos;
        }
    };
};

namespace
{

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
    
    uint16_t excludeFlag;   
    
    ///////
    //i/o//
    ///////
    BAMOrderedReader *odr;
    BCFOrderedWriter *odw;
    bam1_t *s;
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
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_sample_id("s", "ss", "sample ID", true, "", "str", cmd);
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
        excludeFlag = 0x0704;                   

        odr = new BAMOrderedReader(input_bam_file, intervals);
        s = bam_init1();
   
        odw = new BCFOrderedWriter(output_vcf_file, 1);
        //move contigs over from BAM
        bam_hdr_transfer_contigs_to_bcf_hdr(odr->hdr, odw->hdr);
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=E,Number=1,Type=Integer,Description=\"Number of reads containing evidence of the alternate allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Total number of reads at a candidate locus with reads that contain evidence of the alternate allele\">");
        bcf_hdr_add_sample(odw->hdr, sample_id.c_str());
        v = bcf_init();
        
        std::vector<std::string> variant_types;
        split(variant_types, ",", variant_type);
        uint32_t vtype = 0;
        for (uint32_t i = 0; i<variant_types.size(); ++i)
        {
            if (variant_types[i] == "snps")
            {
                vtype |= 1;
            }
            else if (variant_types[i] == "mnps")
            {
                vtype |= 2;
            }
            else if (variant_types[i] == "indels")
            {
                vtype |= 4;
            }
            else if (variant_types[i] == "all")
            {
                vtype = 7;
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
        
        variantHunter = new VariantHunter(vtype,
                                    evidence_allele_count_cutoff, 
                                    fractional_evidence_allele_count_cutoff, 
                                    baseq_cutoff, 
                                    fai, 
                                    odw);
    }

    /**
     * expands cigar string
     */
    void expand_cigar(kstring_t *cigar, kstring_t *cigarString)
    {
        cigar->l = 0;
        int32_t lastIndex = cigarString->l;
        int32_t i = 0;
        kstring_t token = {0,0,0};
        
        if (lastIndex<0)
        {
            return;
        }
        char c;
        bool seenM = false;
        
        while (i<=lastIndex)
        {
            c = cigarString->s[i];
            
            //captures the numeric count
            if (c<'A')
            {
                kputc(c, &token);
            }
    
            if (c>'A' || i==lastIndex) 
            {
                //it is possible for I's to be observed before the first M's in the cigar string
                //in this case, we treat them as 'S'
                if (!seenM)
                {
                    if (c=='I')
                    {
                        c = 'S';
                    }
                    else if (c=='M')
                    {
                        seenM = true;
                    }    
                }
                
                int32_t count = atoi(token.s);
                for (uint32_t j=0; j<count; ++j)
                    kputc(c, cigar);
                token.l = 0;;
            }
            
            ++i;
        }
        
        if (token.m) free(token.s);
            
    };
    
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
            if (bam_get_mpos1(s))
            {
                //overlapping
                if ((bam_get_tid(s)==bam_get_mtid(s)) && bam_get_mpos1(s)<=(bam_get_pos1(s) + bam_get_l_qseq(s) - 1))
                {
                    //add read that has overlapping 
                    //duplicate the record and perform the stitching later
                    k = kh_put(rdict, reads, bam_get_qname(s), &ret);
                    kh_val(reads, k) = {0,0};
                    kh_val(reads, k).start1 = bam_get_pos1(s);
                    kh_val(reads, k).end1 = bam_get_pos1(s) + bam_get_l_qseq(s) - 1;
                }
            }
            else 
            {
                //check overlap
                //todo: perform stitching in future
                if((k = kh_get(rdict, reads, bam_get_qname(s)))!=kh_end(reads))
                {
                    int32_t start1 = bam_get_pos1(s);
                    int32_t end1 = bam_get_pos1(s)+bam_get_l_qseq(s);
                    ++no_overlapping_reads;
                    //std::cerr << "deleting overlapping read\n";
                    //std::cerr << "[" << kh_val(reads, k).start1 << "," << kh_val(reads, k).end1 << "] [" <<  start1 << "," << end1 << "]\n";
                    kh_del(rdict, reads, k); 
                    
                    //skip for the time being
                    continue;
                } 
            }
            
            uint16_t flag = bam_get_flag(s);
            if(flag & excludeFlag)
            {
                ++no_exclude_flag_reads;
                //1. unmapped
                //2. secondary alignment
                //3. not passing QC
                //4. PCR or optical duplicate
                continue;
            }
            
            uint32_t mapq = bam_get_mapq(s);
            if (mapq<mapq_cutoff)
            {
                ++no_low_mapq_reads;
                //filter short aligments and those with too many indels (?)
                continue;
            }

            if (0)
            {
                const char* chrom = bam_get_chrom(odr->hdr, s);
                uint32_t pos1 = bam_get_pos1(s);
                kstring_t seq = {0,0,0};
                bam_get_seq_string(s, &seq);
                uint32_t len = bam_get_l_qseq(s);
                kstring_t qual = {0,0,0};
                bam_get_qual_string(s, &qual);
                kstring_t cigar = {0,0,0};
                bam_get_cigar_string(s, &cigar);
                kstring_t cigar_string = {0,0,0};
                expand_cigar(&cigar_string, &cigar);
                mapq = bam_get_mapq(s);
                
                //if (strchr(cigar.s, 'I') || strchr(cigar.s, 'D'))
                //if (strchr(seq.s, 'N') && strchr(cigar.s, 'D'))
                {
                    std::cerr << "##################" << "\n";
                    std::cerr << "read no  : " << no_reads << "\n";        
                    std::cerr << "chrom-pos: " << chrom << "-" << pos1 << "\n";
                    std::cerr << "read     : " << seq.s << "\n";
                    std::cerr << "qual     : " << qual.s << "\n";
                    std::cerr << "cigar_str: " << cigar_string.s << "\n";
                    std::cerr << "cigar    : " << cigar.s << "\n";
                    std::cerr << "len      : " << len << "\n";
                    std::cerr << "mapq     : " << mapq << "\n";
                    std::cerr << "rnext    : " << bam_get_mtid(s) << "\n";
                }
//              std::cout << "len  :" << cigar.size() << "\n";
//              std::cout << "cigar:" << cigar << "\n";
//              std::cout << "ref  :" << genome_seq << "\n";
            
                if (seq.m) free(seq.s);
                if (qual.m) free(qual.s);
                if (cigar.m) free(cigar.s);
                if (cigar_string.m) free(cigar_string.s);
            }
            
//            if (!strcmp(chrom, current_chrom))
//            {
//                current_chrom = chrom;
//                variantHunter->reset(current_chrom);
//                kh_clear(rdict, reads);
//            }
            
            variantHunter->process_read(odr->hdr, s);
            
            if (no_reads%100000==0) std::cerr << no_reads << "\n";
           
            ++no_passed_reads;
        }
    };

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
        std::clog << "         [i] intervals                    ";
        for (uint32_t i=0; i<std::min((uint32_t)intervals.size(),(uint32_t)5); ++i)
        {
            if (i) std::clog << ", ";
            std::clog << intervals[i].to_string();
        }
        if (intervals.size()>=5)
        {
            std::clog << "  and " << (intervals.size()-5) <<  " other intervals\n";
        }    
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
