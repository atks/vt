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

#include "augmented_bam_record.h"

/**
 * Constructor.
 */
AugmentedBAMRecord::AugmentedBAMRecord()
{
    clear();
}

/**
 * Constructor.
 */
AugmentedBAMRecord::AugmentedBAMRecord(bam1_t* s)
{
    clear();
    initialize(s);
}

/**
 * Initialize.
 */
void AugmentedBAMRecord::initialize(bam1_t* s)
{
    this->s = s;

    uint32_t *cigar = bam_get_cigar(s);
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    char opchr;
    int32_t oplen;

    aug_cigar.clear();

    //get MD tag
    uint8_t *md_aux;
    char* md = 0;
    ((md_aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(md_aux)));
    
    //this points to the part of md that is yet to be processed.
    char* mdp = md;

    uint32_t cpos1 = pos1; //current 1 based genome position
    uint32_t spos0 = 0;    //current position in read sequence

    //variables for I's embedded in Matches in the MD tag
    uint32_t md_mlen_left = 0;
    bool seenM = false;

    for (int32_t i = 0; i < n_cigar_op; ++i)
    {
        opchr = bam_cigar_opchr(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);

        //variables for I's embedded in Matches in the MD tag
        uint32_t md_mlen_left = 0;

        if (opchr=='S')
        {
            spos0 += oplen;
            aug_cigar.push_back(cigar[i]);
            aug_seq.push_back(NULL);
        }
        else if (opchr=='M')
        {
            uint32_t lpos1 = cpos1; // we need this because M contains matches and mismatches
            uint32_t sspos0 = spos0; // we need this because M contains matches and mismatches
            uint32_t mlen = oplen;
            uint32_t i = 0;
            seenM = true;

            //left over MD matches to handle.
            if (md_mlen_left)
            {
                uint32_t ilen = md_mlen_left<=mlen ? md_mlen_left : mlen;

                lpos1 += ilen;
                sspos0 += ilen;

                //yet another insertion
                if (md_mlen_left>=mlen)
                {
                    md_mlen_left -= ilen;
                    cpos1 += ilen;
                    spos0 += ilen;
                    continue;
                }
                //a snp next
                else
                {
                    md_mlen_left = 0;
                    mlen -= ilen;
                    //go to loop in the next section
                }
            }

            while (*mdp)
            {
                if (isalpha(*mdp)) //SNPs
                {
                    char ref = toupper(*mdp);
                    char alt = (bam_base2char(bam_seqi(seq, spos0+(lpos1-cpos1))));

                    ++lpos1;
                    ++mdp;
                    ++sspos0;
                    --mlen;
                }
                else if (isdigit(*mdp)) //matches
                {
                    char* end = 0;
                    int32_t len = std::strtol(mdp, &end, 10);
                    mdp = end;

                    if (len)
                    {
                        uint32_t ilen = len<=mlen ? len : mlen;

                        lpos1 += ilen;
                        sspos0 += ilen;

                        //next up an insertion
                        if (len>mlen)
                        {
                            md_mlen_left = len - mlen;
                            break;
                        }
                        else
                        {
                            mlen -= ilen;
                        }
                    }
                }
                else // deletion
                {
                    break;
                }

                if (mlen==0)
                {
                    break;
                }

                //note that only insertions, matches and mismatches can only occur here.

                cpos1 += oplen;
                spos0 += oplen;
            }
        }
        else if (opchr=='D')
        {
            bool is_del = false;

            if (*mdp=='0') ++mdp;

            if (*mdp!='^')
            {
                std::cerr << "mdp: " << mdp << "\n";
                std::cerr << "inconsistent MD and cigar, deletion does not occur at the right place.\n";
                exit(1);
            }

            ++mdp;
            std::string del = "";
            while (isalpha(*mdp))
            {
                del += toupper(*mdp);
                ++mdp;
            }

            cpos1 += oplen;
        }
        else if (opchr=='I')
        {
            //leading Is
            if (!seenM)
            {
                spos0 += oplen;
            }
            //trailing Is
            else if (i==n_cigar_op-1 || (i+2==n_cigar_op && bam_cigar_opchr(cigar[n_cigar_op-1])=='S'))
            {
                spos0 += oplen;
            }
            else
            {
                //insertions are not present in MD tags
                //may be handled independently of future matches
                std::string ins = "";
                for (size_t i=0; i<oplen ; ++i)
                {
                    ins += bam_base2char(bam_seqi(seq, spos0+i));
                }

                spos0 += oplen;
            }
        }
        else
        {
            std::cerr << "never seen before state " << opchr << "\n";
            exit(1);
        }

       // aug_cigar
        std::cerr << oplen << ((char)opchr);
    }
}

/**
 * left_align augmented cigar.
 */
bool AugmentedBAMRecord::left_align()
{
    return true;
}

/**
 * Clear.
 */
void AugmentedBAMRecord::clear()
{
    s = NULL;
    seq = NULL;
    cigar = NULL;
    md = NULL;
    pos1 = 0;
    aug_cigar.clear();
    aug_seq.clear();
}

/**
 * Prints alignment of record.
 */
void AugmentedBAMRecord::print()
{
    int32_t oplen;
    char opchr;
    
    std::string ref;
    std::string align;
    std::string seq;
    int32_t spos0 = 0;     
    
    for (uint32_t i=0; i<aug_cigar.size(); ++i)
    {
        oplen = bam_cigar_oplen(aug_cigar[i]);
        opchr = bam_cigar_opchr(aug_cigar[i]);
        
        if (opchr=='S')
        {
            ref.append('-', oplen);
            
            for (uint32_t j=0; j<oplen; ++j)
            {
                seq.append(bam_base2char(bam_seqi(this->seq, spos0+j)), 1);
            }
            
            align.append('S', oplen);
        }    
        else if (opchr=='=')
        {
            for (uint32_t j=0; j<oplen; ++j)
            {
                ref.append(bam_base2char(bam_seqi(this->seq, spos0+j)), 1);
                seq.append(bam_base2char(bam_seqi(this->seq, spos0+j)), 1);
            }
        
            align.append('=', oplen);
            
            spos0 += oplen;
        }
        else if (opchr=='X')
        {
            for (uint32_t j=0; j<oplen; ++j)
            {
                ref.append(aug_seq[j], 1);
                seq.append(bam_base2char(bam_seqi(this->seq, spos0+j)), 1);
            }
        
            align.append('X', oplen);
            
            spos0 += oplen;
        } 
        else if (opchr=='I')
        {
            ref.append('-', oplen);
            
            for (uint32_t j=0; j<oplen; ++j)
            {
                seq.append(bam_base2char(bam_seqi(this->seq, spos0+j)), 1);
            }
            
            align.append('I', oplen);
            
            spos0 += oplen;
        }
        else if (opchr=='D')
        {
            for (uint32_t j=0; j<oplen; ++j)
            {
                ref.append(aug_seq[j], oplen);
            }
            
            seq.append('-', oplen);
            
            align.append('D', oplen);
        }
        else
        {
            std::cerr << "unrecognized cigar state " << opchr << "\n";
//            exit(1);
        }
        
    }
    
    std::cerr << ref << "\n";
    std::cerr << align << "\n";
    std::cerr << seq << "\n";        
}