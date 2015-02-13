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

#include "read_alignment.h"

/**
 * Constructor
 */
ReadAlignment::ReadAlignment(){};

/**
 * Sets the read object with a new read alignment.
 */
void ReadAlignment::set(bam1_t *s, uint8_t* seq, uint8_t* qual, int32_t l_qseq, uint32_t* cigar, char* md)
{
    this->seq = seq;
    this->qual = qual;
    this->l_qseq = l_qseq;
    this->cigar = cigar;
    this->md = md;

    kstring_t str = {0,0,0};

    std::cerr << "*************" << "\n";
    int32_t pos1 = bam_get_pos1(s);
    std::cerr << "\tpos1     : " << pos1 << "\n";
    bam_get_seq_string(s, &str);
    std::cerr << "\tread     : " << str.s << "\n";
    bam_get_qual_string(s, &str);
    std::cerr << "\tqual     : " << str.s << "\n";
    bam_get_cigar_string(s, &str);
    std::cerr << "\tcigar_str: " << str.s << "\n";
    bam_get_cigar_expanded_string(s, &str);
    std::cerr << "\tcigar    : " << str.s << "\n";
    std::cerr << "\tlen      : " << l_qseq << "\n";
    uint8_t *aux;
    char* md1;
    (aux=bam_aux_get(s, "MD")) &&  (md1 = bam_aux2Z(aux));
    std::cerr << "\tmd       : " << md1 << "\n";

    if (str.m) free(str.s);


    //iterate cigar
    int32_t n_cigar_op = bam_get_n_cigar_op(s);

    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(s);
        for (int32_t i = 0; i < n_cigar_op; ++i)
        {
            std::cerr << bam_cigar_oplen(cigar[i]) << " " << bam_cigar_opchr(cigar[i]) << "\n";
        }
    }

    char* mdp = md;
    bool is_del = false;
    while (*mdp)
    {
        //M
        if (isdigit(*mdp))
        {
            char* end = 0;
            int32_t len = std::strtol(mdp, &end, 10);
            mdp = end;
            is_del = false;
            std::cerr << "Match " << len << "\n";
        }
        //reference base
        else if (isalpha(*mdp))
        {
            if (is_del)
            {
                std::cerr << "Deletion " << *mdp << "\n";
                ++mdp;
            }
            else
            {
                std::cerr << "Mismatch " << *mdp << "\n";
                ++mdp;
            }
        }
        //deletion
        else if (*mdp=='^')
        {
            ++mdp;
            is_del = true;
        }
    }

    std::cerr << "************ " << "\n";

    bool process_cigar = true;
    int32_t c_len = 0;
    int32_t md_len = 0;
    std::string del = "";
    mdp = md;
    int32_t cpos1 = pos1;
    int32_t spos0 = 0;
    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(s);
        for (int32_t i = 0; i < n_cigar_op; ++i)
        {
            int32_t oplen = bam_cigar_oplen(cigar[i]);
            char opchar = bam_cigar_opchr(cigar[i]);

            std::cerr << oplen << " " << opchar << "\n";

            if (opchar=='S')
            {
                //add to S evidence
                //do nothing
                
                std::string ins = "";
                float mean_qual = 0;   
                for (size_t i=0; i<oplen ; ++i)
                {
                    ins += bam_base2char(bam_seqi(seq, spos0+i));
                    mean_qual += qual[spos0+i];
                }
                
                mean_qual /= oplen;
                
                if (cpos1==pos1)
                {
                    std::cerr << "\t\t\tadding sclip left: " << cpos1 << "\t" << ins << " (" << mean_qual << ")\n";
                }
                else
                {
                    std::cerr << "\t\t\tadding sclip right: " << cpos1 << "\t" << ins << " (" << mean_qual << ")\n";
                }
                
                spos0 += oplen;
            }
            else if (opchar=='M')
            {
                int32_t lpos1 = cpos1;
                while (*mdp)
                {
                    if (isdigit(*mdp))
                    {
                        char* end = 0;
                        int32_t len = std::strtol(mdp, &end, 10);
                        mdp = end;
    
                        std::cerr << "\tMatch " << len << "\n";
                        
                        if (len)
                        {
                            std::cerr << "\t\t\tadding matches: " << lpos1 << "-" << (lpos1+len-1) << "\n";
                            lpos1 += len;
                        }
                    }
                    else if (isalpha(*mdp))
                    {
                        std::cerr << "\tMismatch " << *mdp << "\n";
                        std::cerr << "\t\t\tadding mismatch: " << lpos1 << "\t" << *mdp << "\n";
                    
                    
                        ++lpos1;
                        ++mdp;
                    }
                    else
                    {
                        break;
                    }
                }
                
                //note that only Insertions , matches and mismatches can only occur here.
                
                cpos1 += oplen;
                spos0 += oplen;
            }
            else if (opchar=='D')
            {
                bool is_del = false;

                if (*mdp!='^')
                {
                    std::cerr << "inconsistent MD and cigar\n";
                    exit(1);
                }
                else
                {
                    ++mdp;
                    std::string del = "";
                    while (isalpha(*mdp))
                    {
                        del += *mdp;
                        ++mdp;
                    }

                    std::cerr << "\t\t\tadding deletion: " << cpos1 << " " << del << "\n";
                        
                    cpos1 += oplen;
                    spos0 += oplen;
                }
            }
            else if (opchar=='I')
            {
                //may be handled independently of future matches
                std::cerr << "insertion " << opchar << "\n";

                std::string ins = "";
                for (size_t i=0; i<oplen ; ++i)
                {
                    ins += bam_base2char(bam_seqi(seq, spos0+i));
                }

                std::cerr << "\t\t\tadding insertion: " << cpos1 << " " << ins  << "\n";
                
                spos0 += oplen;    
                //extract sequence from read



                //go back an subtract one evidence from the array for indels.
            }
            else
            {
                std::cerr << "never seen before state " << opchar << "\n";
            }

        }
    }

    //return variants to populate the buffer
    //
    //start-end to fill data
    //actual variants
    //
}

/**
 * Gets the next feature.
 */
void ReadAlignment::get_next_feature(feature_t& f)
{
    //if end, show end

    //get next cigar op cigar



    //get next md op




    //update cigar op
    //update md op
};

void ReadAlignment::next_op(feature_t& f)
{
    //M
    if (isdigit(*md))
    {
        f.type = M;

        const char* start = md;
        char *end = 0;
        f.len = std::strtol(md, &end, 10);
    }
    //reference base
    else if (isalpha(*md))
    {
        f.type = X;
    }
    //deletion
    else if (*md=='^')
    {
        f.type = D;
    }
};