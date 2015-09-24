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
AugmentedBAMRecord::AugmentedBAMRecord(bam_hdr_t* h, bam1_t* s)
{
    clear();
    initialize(h, s);
}

/**
 * Initialize.
 */
void AugmentedBAMRecord::initialize(bam_hdr_t* h, bam1_t* s)
{
    clear();

    this->h = h;
    this->s = s;
    beg1  = bam_get_pos1(s);
    end1  = beg1;

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
    uint32_t rpos0 = 0;                //current 0-based position in read sequence and qual field
    uint32_t md_mlen_left = 0;         //unprocessed MD matches that embeds insertions
    bool seenM = false;                //to check read is aligned
    bool expectedI = false;             //true when the next state expected is I, for error checking

    for (int32_t i = 0; i < n_cigar_op; ++i)
    {
        opchr = bam_cigar_opchr(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);

        if (expectedI && opchr!='I')
        {
            print();
            fprintf(stderr, "[%s:%d %s] Inconsistent CIGAR and MD.\n", __FILE__, __LINE__, __FUNCTION__);
        }

        if (opchr=='S')
        {
            aug_cigar.push_back(cigar[i]);
            aug_ref.push_back("");
            aug_alt.push_back("");

            rpos0 += oplen;
        }
        else if (opchr=='M')
        {
            seenM = true;
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
                    aug_cigar.push_back(bam_cigar_gen(mlen, BAM_CEQUAL));
                    aug_ref.push_back("");
                    aug_alt.push_back("");

                    md_mlen_left -= mlen;
                    rpos0 += mlen;

                    //no need to process MD tag since we are in the midst of perfect matches
                    //we skip to process the next insertion cigar operation
                    expectedI = true;
                    mlen = 0;
                    continue;
                }
                //this is very important
                //CIGAR: 12M
                //MD   : 6A5
                //need to process MD tag
                //to process the next insertion or SNP
                else
                {
                    aug_cigar.push_back(bam_cigar_gen(md_mlen_left, BAM_CEQUAL));
                    aug_ref.push_back("");
                    aug_alt.push_back("");

                    mlen -= md_mlen_left;
                    rpos0 += md_mlen_left;

                    md_mlen_left = 0;
                    //go to loop in the next section
                }
            }

            //might have multiple mismatches
            while (*mdp)
            {
                if (isalpha(*mdp)) //mismatches
                {
                    aug_cigar.push_back(bam_cigar_gen(1, BAM_CDIFF));

//                    std::cerr << "ADDING " << ((char)toupper(*mdp)) << " " << bam_base2char(bam_seqi(seq, rpos0)) << "\n";

                    aug_ref.push_back(std::string(1, toupper(*mdp)));
                    aug_alt.push_back(std::string(1, bam_base2char(bam_seqi(seq, rpos0))));

                    ++mdp;
                    ++rpos0;
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
                            aug_cigar.push_back(bam_cigar_gen(mlen, BAM_CEQUAL));
                            aug_ref.push_back("");
                            aug_alt.push_back("");

                            md_mlen_left = len - mlen;
                            rpos0 += mlen;
                            mlen = 0;

                            expectedI = true;
                            break;
                        }
                        //another mismatch
                        else
                        {
                            aug_cigar.push_back(bam_cigar_gen(len, BAM_CEQUAL));
                            aug_ref.push_back("");
                            aug_alt.push_back("");

                            mlen -= len;
                            rpos0 += len;

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
            expectedI = false;

            //leading Is
            if (!seenM)
            {
                //convert to Ss
                aug_cigar.push_back(bam_cigar_gen(oplen, BAM_CSOFT_CLIP));
                aug_ref.push_back("");
                aug_alt.push_back("");

                rpos0 += oplen;
            }
            //trailing Is
            else if (i==n_cigar_op-1 || (i+2==n_cigar_op && bam_cigar_opchr(cigar[n_cigar_op-1])=='S'))
            {
                aug_cigar.push_back(bam_cigar_gen(oplen, BAM_CSOFT_CLIP));
                aug_ref.push_back("");
                aug_alt.push_back("");

                rpos0 += oplen;
            }
            else
            {
                //insertions are not present in MD tags
                //may be handled independently of future matches
                std::string ins = "";
                for (size_t i=0; i<oplen ; ++i)
                {
                    ins += bam_base2char(bam_seqi(seq, rpos0+i));
                }

                aug_cigar.push_back(bam_cigar_gen(oplen, BAM_CINS));
                aug_ref.push_back("");
                aug_alt.push_back(ins);

                rpos0 += oplen;
            }
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cigar state not handled %c\n", __FILE__, __LINE__, __FUNCTION__, opchr);
            exit(1);
        }
    }

    //count number of mismatches
    //compute end1 position of alignment
    uint8_t* qual = bam_get_qual(s);
    rpos0 = 0;
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
            rpos0 += oplen;
            end1 +=oplen;
        }
        else if (opchr=='X')
        {
            if (qual[rpos0]>=20)
            {
                ++no_mismatches;
            }
            ++rpos0;
            ++end1;
        }
        else if (opchr=='I')
        {
            rpos0 += oplen;
            ++no_mismatches;
        }
        else if (opchr=='D')
        {
            ++no_mismatches;
            end1 +=oplen;
        }
        else
        {
            std::cerr << "unrecognized cigar state " << opchr << "\n";
        }
    }
    
    --end1;

     print();
}

/**
 * Left align indels in augmented cigar.
 */
bool AugmentedBAMRecord::left_align()
{
    //chec if any indel is not left aligned!


    return true;
}

/**
 * Right align indels in augmented cigar.
 */
bool AugmentedBAMRecord::right_align()
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
    no_mismatches = 0;
}

/**
 * Prints alignment of record.
 */
void AugmentedBAMRecord::print()
{
    return;
    bool has_indels = false;

    if (true)
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

        if (strchr(md, '^'))
        {
            has_indels = true;
        }

        for (uint32_t i=0; i<aug_cigar.size(); ++i)
        {
            int32_t oplen = bam_cigar_oplen(aug_cigar[i]);
            char opchr = bam_cigar_opchr(aug_cigar[i]);

            if (opchr=='I' || opchr=='D')
            {
                has_indels = true;
                break;
            }

        }

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
        std::cerr << "nm       : " << no_mismatches << "\n";

        int32_t nm = 0;
        (aux=bam_aux_get(s, "NM")) &&  (nm = bam_aux2i(aux));
        std::cerr << "bam nm   : " << nm << "\n";

        std::cerr << "##################" << "\n";

        if (seq.m) free(seq.s);
        if (qual.m) free(qual.s);
        if (cigar_string.m) free(cigar_string.s);
        if (cigar_expanded_string.m) free(cigar_expanded_string.s);


        if (strchr(md, '^'))
        {
            has_indels = true;
        }

        for (uint32_t i=0; i<aug_cigar.size(); ++i)
        {
            int32_t oplen = bam_cigar_oplen(aug_cigar[i]);
            char opchr = bam_cigar_opchr(aug_cigar[i]);

            if (opchr=='I' || opchr=='D')
            {
                has_indels = true;
                break;
            }
        }
    }

    std::cerr << "##################" << "\n";
    std::cerr << "Augmented CIGAR" << "\n";
    std::cerr << "##################" << "\n";

    int32_t oplen;
    char opchr;

    std::cerr << "AUG_CIGAR : ";
    for (uint32_t i=0; i<aug_cigar.size(); ++i)
    {
        oplen = bam_cigar_oplen(aug_cigar[i]);
        opchr = bam_cigar_opchr(aug_cigar[i]);

        std::cerr << oplen << opchr;

    }
    std::cerr << "\n";

    std::string ref;
    std::string align;
    std::string seq;
    std::string quals;
    int32_t spos0 = 0;

    uint8_t* qual = bam_get_qual(s);
    for (uint32_t i=0; i<aug_cigar.size(); ++i)
    {
        oplen = bam_cigar_oplen(aug_cigar[i]);
        opchr = bam_cigar_opchr(aug_cigar[i]);

        if (opchr=='S')
        {
            ref.append(oplen, '-');

            for (uint32_t j=0; j<oplen; ++j)
            {
                seq.append(1, bam_base2char(bam_seqi(bam_get_seq(this->s), spos0+j)));
                quals.append(1, qual[spos0+j]+33);
            }

            align.append(oplen, 'S');
        }
        else if (opchr=='=')
        {
            for (uint32_t j=0; j<oplen; ++j)
            {
                ref.append(1, bam_base2char(bam_seqi(bam_get_seq(this->s), spos0+j)));
                seq.append(1, bam_base2char(bam_seqi(bam_get_seq(this->s), spos0+j)));
                quals.append(1, qual[spos0+j]+33);
            }
            align.append(oplen, '=');

            spos0 += oplen;
        }
        else if (opchr=='X')
        {
            //assume oplen is always 1.

            ref.append(aug_ref[i]);
            seq.append(aug_alt[i]);
            align.append(1, 'X');
            quals.append(1, qual[spos0]+33);

            spos0 += 1;
        }
        else if (opchr=='I')
        {
            ref.append(oplen, '-');
            seq.append(aug_alt[i]);
            align.append( oplen, 'I');

            for (uint32_t j=0; j<oplen; ++j)
            {
                quals.append(1, qual[spos0+j]+33);
            }

            spos0 += oplen;
        }
        else if (opchr=='D')
        {
            ref.append(aug_ref[i]);
            seq.append(oplen, '-');
            align.append(oplen, 'D');
        }
        else
        {
            std::cerr << "unrecognized cigar state " << opchr << "\n";
        }
    }

    std::cerr << "REF   : "<< ref << "\n";
    std::cerr << "ALIGN : "<< align << "\n";
    std::cerr << "READ  : "<< seq << "\n";
    std::cerr << "QUAL  : "<< quals << "\n";
}