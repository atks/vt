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

    uint8_t* seq = bam_get_seq(s);

    //CIGAR related variables
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    uint32_t *cigar = bam_get_cigar(s);
    char opchr;
    int32_t oplen;

    //get MD tag
    uint8_t *md_aux;
    char* md = 0;
    ((md_aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(md_aux)));
    char* mdp = md; //pointer to md


    //variables for keep track of CIGAR and MD tag
    uint32_t cpos1 = bam_get_pos1(s); //current 1-based genome position
    uint32_t spos0 = 0;                //current 0-based position in read sequence and qual field
    uint32_t md_mlen_left = 0;         //unprocessed MD matches that embeds insertions
    bool seenM = false;                //to check read is aligned
    bool expectedI = true;             //true when the next state expected is I, for error checking

    for (int32_t i = 0; i < n_cigar_op; ++i)
    {
        opchr = bam_cigar_opchr(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);

        if (expectedI && opchr!='I')
        {
            fprintf(stderr, "[%s:%d %s] Inconsistent CIGAR and MD.\n", __FILE__, __LINE__, __FUNCTION__);
            //print warning
        }
        
        if (opchr=='S')
        {
            spos0 += oplen;
            aug_cigar.push_back(cigar[i]);
            aug_ref.push_back("");
            aug_alt.push_back("");
        }
        else if (opchr=='M')
        {
            seenM = true;

            uint32_t lpos1 = cpos1; // we need this because M contains matches and mismatches
            uint32_t sspos0 = spos0; // we need this because M contains matches and mismatches
            
            uint32_t mlen = oplen; //match length
    
            //left over MD matches to handle.
            //this occurs when I is embedded.
            if (md_mlen_left)
            {
                //cigar I embedded in MD matches
                //CIGAR: 6M4I6M
                //MD   : 12
                //to process the next insertion or SNP
                if (md_mlen_left>mlen)
                {
                    md_mlen_left -= mlen;
                    cpos1 += mlen;
                    spos0 += mlen;
                    
                    //no need to process MD tag since we are in the midst of perfect matches
                    //we skip to process the next insertion cigar operation
                    expectedI = true;
                    continue;
                }
                //this is very important
                //CIGAR: 12M
                //MD   : 6A5
                //need to process MD tag
                //to process the next insertion or SNP
                else
                {
                    md_mlen_left = 0;
                    mlen -= md_mlen_left;
                    //go to loop in the next section
                }
            }

            //might have multiple mismatches
            while (*mdp)
            {
                if (isalpha(*mdp)) //mismatches
                {
                    aug_cigar.push_back(bam_cigar_gen(1, BAM_CDIFF));
                    aug_ref.push_back(std::string(1, toupper(*mdp)));
                    aug_alt.push_back(std::string(1, bam_base2char(bam_seqi(seq, spos0))));

                    ++mdp;
                    ++spos0;
                    --mlen;
                }
                else if (isdigit(*mdp)) //matches
                {
                    char* end = 0;
                    int32_t len = std::strtol(mdp, &end, 10);
                    mdp = end;

                    if (len)
                    {
                        //another I
                        if (len>mlen)
                        {
                            md_mlen_left = len - mlen;
                            expectedI = true;
                            mlen = 0;
                            
                            aug_cigar.push_back(bam_cigar_gen(mlen, BAM_CEQUAL));
                            aug_ref.push_back("");
                            aug_alt.push_back("");
                            
                            break;
                        }
                        //another mismatch
                        else
                        {
                            mlen -= len;
                            
                            aug_cigar.push_back(bam_cigar_gen(len, BAM_CEQUAL));
                            aug_ref.push_back("");
                            aug_alt.push_back("");
                            
                            //continue processing MD tag
                        }
                    }
                }
                else // deletion, to be handled in the cigar's D operation
                {
                    expectedI = false;
                    break;
                }
            }
        }
        else if (opchr=='D')
        {
            bool is_del = false;

            if (*mdp=='0') ++mdp;

            if (*mdp!='^')
            {
                fprintf(stderr, "[%s:%d %s] Inconsistent CIGAR and MD where deletion is expected\n", __FILE__, __LINE__, __FUNCTION__);
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

            aug_cigar.push_back(cigar[i]);
            aug_ref.push_back(del);
            aug_alt.push_back("");
        }
        else if (opchr=='I')
        {
            //leading Is
            if (!seenM)
            {
                //convert to Ss
                aug_cigar.push_back(cigar[i]);
                
                
                aug_cigar.push_back(bam_cigar_gen(oplen, BAM_CSOFT_CLIP));
                aug_ref.push_back("");
                aug_alt.push_back("");
                spos0 += oplen;
            }
            //trailing Is
            else if (i==n_cigar_op-1 || (i+2==n_cigar_op && bam_cigar_opchr(cigar[n_cigar_op-1])=='S'))
            {
                aug_cigar.push_back(bam_cigar_gen(oplen, BAM_CSOFT_CLIP));
                aug_ref.push_back("");
                aug_alt.push_back("");
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

                aug_cigar.push_back(bam_cigar_gen(oplen, BAM_CINS));
                aug_ref.push_back("");
                aug_alt.push_back(ins);
            
                spos0 += oplen;
            }
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cigar state not handled %c\n", __FILE__, __LINE__, __FUNCTION__, opchr);
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
    aug_cigar.clear();
    aug_ref.clear();
    aug_alt.clear();
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
                seq.append(bam_base2char(bam_seqi(bam_get_seq(this->s), spos0+j)), 1);
            }
            
            align.append('S', oplen);
        }    
        else if (opchr=='=')
        {
            for (uint32_t j=0; j<oplen; ++j)
            {
                ref.append(bam_base2char(bam_seqi(bam_get_seq(this->s), spos0+j)), 1);
                seq.append(bam_base2char(bam_seqi(bam_get_seq(this->s), spos0+j)), 1);
            }
            align.append('=', oplen);
            
            spos0 += oplen;
        }
        else if (opchr=='X')
        {
            //assume oplen is always 1.
            
            ref.append(aug_ref[i]);
            seq.append(aug_alt[i]);
            align.append('X', 1);
            
            spos0 += 1;
        } 
        else if (opchr=='I')
        {
            ref.append('-', oplen);
            seq.append(aug_alt[i]);
            align.append('I', oplen);
            
            spos0 += oplen;
        }
        else if (opchr=='D')
        {
            ref.append(aug_ref[i]);
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