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

#include "pileup.h"

/**
 * Clears the soft clipped information.
 */
void SoftClipInfo::clear()
{
    no = 0;
    mean_quals.clear();
    strands.clear();
}

/**
 * Constructor.
 */
PileupPosition::PileupPosition()
{
    clear();
};

/**
 * Clears pileup position.
 */
void PileupPosition::clear()
{
    R = 'X';
    X[1] = 0;
    X[2] = 0;
    X[4] = 0;
    X[8] = 0;
    X[15] = 0;
    D.clear();
    I.clear();
    J.clear();
    K.clear();
    N = 0;
    F = 0;
    E = 0;
    REF_Q.clear();
    ALT_Q.clear();
};

/**
 * Returns true if pileup record is cleared.
 */
bool PileupPosition::is_cleared()
{
    return (R=='X' &&
    X[1] == 0 &&
    X[2] == 0 &&
    X[4] == 0 &&
    X[8] == 0 &&
    X[15] == 0 &&
    D.empty() &&
    I.empty() &&
    J.empty() &&
    K.empty() &&
    N == 0 &&
    F == 0 &&
    E == 0 &&
    REF_Q.empty()&&
    ALT_Q.empty());
};

/**
 * Prints pileup position.
 */
void PileupPosition::print()
{
    std::cerr << R << ":" << N << "+" << E << "\n";

    if (X[1]+X[2]+X[4]+X[8]+X[15])
    {
        std::cerr << "\tSNP: ";
        if (X[1])
        {
            std::cerr << X[1] << "A";
        }
        if (X[2])
        {
            std::cerr << X[2] << "C";
        }
        if (X[4])
        {
            std::cerr << X[4] << "G";
        }
        if (X[8])
        {
            std::cerr << X[8] << "T";
        }
        if (X[15])
        {
            std::cerr << X[15] << "N";
        }
        std::cerr << "\n";
    }

    if (D.size()!=0)
    {
        std::cerr << "\tDEL: ";
        for (std::map<std::string, uint32_t>::iterator i = D.begin(); i!=D.end(); ++i)
        {
            std::cerr << i->first << " (" << i->second << ")";
        }
        std::cerr << "\n";
    }
    
    if (I.size()!=0)
    {
        std::cerr << "\tINS: ";
        for (std::map<std::string, uint32_t>::iterator i = I.begin(); i!=I.end(); ++i)
        {
            std::cerr << i->first << " (" << i->second << ")";
        }

        std::cerr << "\n";
    }
    
    if (J.size()!=0)
    {
        std::cerr << "\tRSCLIP: ";
        for (std::map<std::string, SoftClipInfo>::iterator i = J.begin(); i!=J.end(); ++i)
        {
            std::cerr << i->first << " (" << i->second.no << ")";
        }

        std::cerr << "\n";
    }
    
    if (K.size()!=0)
    {
        std::cerr << "\tLSCLIP: ";
        for (std::map<std::string, SoftClipInfo>::iterator i = K.begin(); i!=K.end(); ++i)
        {
            std::cerr << i->first << " (" << i->second.no << ")";
        }

        std::cerr << "\n";
    }
}

/**
 * Prints pileup position.
 */
void PileupPosition::print(uint32_t gpos1)
{
    std::cerr << gpos1 << ":";
    print();
}

/**
 * Constructor.
 *
 * @k - size of pileup is 2^k
 */
Pileup::Pileup(uint32_t k, uint32_t window_size)
{
    //Buffer size is a power of 2^k.
    buffer_size = 1 << k;
    //this provides a cheaper way to do modulo operations for a circular array.
    buffer_size_mask = (0xFFFFFFFF >> (32-k));
    this->window_size = window_size;

    P.resize(buffer_size);

    tid = -1;
    beg0 = end0 = 0;
    gbeg1 = 0;

    debug = 0;
};

/**
 * Overloads subscript operator for accessing pileup positions.
 */
PileupPosition& Pileup::operator[] (const int32_t i)
{
    return P[i];
}

/**
 * Returns the maximum size of the pileup.
 */
uint32_t Pileup::max_size()
{
    return buffer_size - 1;
}

/**
 * Returns the size of the pileup.
 */
uint32_t Pileup::size()
{
    return (end0>=beg0 ? end0-beg0 : buffer_size-(beg0-end0));
}

/**
 * Checks if buffer is empty.
 */
bool Pileup::is_empty()
{
    return beg0==end0;
};

/**
 * Set reference fasta file.
 */
void Pileup::set_reference(std::string& ref_fasta_file)
{
    if (ref_fasta_file!="")
    {
        fai = fai_load(ref_fasta_file.c_str());
        if (fai==NULL)
        {
            fprintf(stderr, "[%s:%d %s] cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
            exit(1);
        }
    }
};

/**
 * Set debug.
 */
void Pileup::set_debug(int32_t debug)
{
    this->debug = debug;
};

/**
 * Sets tid.
 */
void Pileup::set_tid(uint32_t tid)
{
    this->tid = tid;
}

/**
 * Gets tid.
 */
uint32_t Pileup::get_tid()
{
    return this->tid;
}

/**
 * Sets chrom.
 */
void Pileup::set_chrom(std::string& chrom)
{
    this->chrom = chrom;
}

/**
 * Gets chrom.
 */
std::string Pileup::get_chrom()
{
    return chrom;
}

/**
 * Gets window_size.
 */
uint32_t Pileup::get_window_size()
{
    return window_size;
}

/**
 * Converts gpos1 to index in P.
 * If P is empty, initialize first position as gpos1.
 */
uint32_t Pileup::g2i(uint32_t gpos1)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        return beg0;
    }
    else
    {
        if (gpos1<gbeg1)
        {
            fprintf(stderr, "[%s:%d %s] buffer out of extent: gpos1 %d < gbeg1 %d\n", __FILE__, __LINE__, __FUNCTION__, gpos1, gbeg1);
            abort();
        }
        return (beg0 + (gpos1-gbeg1)) & buffer_size_mask;
    }
}

/**
 * Increments i by 1 circularly.
 */
uint32_t Pileup::inc(uint32_t i)
{
    return (i+1) & buffer_size_mask;
};

/**
 * Sets gbeg1.
 */
void Pileup::set_gbeg1(uint32_t gbeg1)
{
    this->gbeg1 = gbeg1;
}

/**
 * Gets gbeg1.
 */
uint32_t Pileup::get_gbeg1()
{
    return gbeg1;
}

/**
 * Gets gend1.
 */
uint32_t Pileup::get_gend1()
{
    if (is_empty())
    {
        return 0;
    }
    else
    {
        return gbeg1 + diff(end0, beg0) - 1;
    }
}

/**
 * Sets beg0.
 */
void Pileup::set_beg0(uint32_t beg0)
{
    this->beg0 = beg0;
}

/**
 * Sets end0.
 */
void Pileup::set_end0(uint32_t end0)
{
    this->end0 = end0;
}

/**
 * Gets the index of the first element.
 */
uint32_t Pileup::begin()
{
    return beg0;
}

/**
 * Gets the index of the last element.
 */
uint32_t Pileup::end()
{
    return end0;
}

/**
 * Returns the difference between 2 buffer positions
 */
uint32_t Pileup::diff(uint32_t i, uint32_t j)
{
    return (i>=j ? i-j : buffer_size-(j-i));
};

/**
 * Increments beg0 by 1.
 */
void Pileup::inc_beg0()
{
    beg0 = (beg0+1) & buffer_size_mask;
};

/**
 * Increments end0 by 1.
 */
void Pileup::inc_end0()
{
    end0 = (end0+1) & buffer_size_mask;
};

/**
 * Increments index i by j cyclically.
 */
uint32_t Pileup::inc(uint32_t i, uint32_t j)
{
    return (i+j) & buffer_size_mask;
};
/**
 * Updates the last aligned base in a read.
 *
 * @gpos1 - starting genome position
 */
void Pileup::update_read_end(uint32_t gpos1)
{
    if (debug) std::cerr << "update_read_end_base: gpos1" << gpos1 << "\n";
    uint32_t i = g2i(gpos1);
    if (P[i].N==0)
    {
        std::cerr << "update_read_end_base: gpos1" << gpos1 << "\n";
        std::cerr << "OOPSs...\n";
    }
    else
    {
        --P[i].N;
        ++P[i].E;
    }
}

/**
 * Inserts a stretch of aligned bases identified by M in the cigar string.
 */
void Pileup::add_M(uint32_t mgpos1, uint32_t spos0, uint32_t len, uint8_t* seq, uint8_t* qual, uint32_t snp_baseq_cutoff)
{
    add_3prime_padding(mgpos1);

    uint32_t gend1 = get_gend1();
    uint32_t mgend1 = mgpos1 + len - 1;
    uint32_t i = g2i(mgpos1);
    uint32_t j = 0;

//    std::cerr << "addM: \n";
//    std::cerr << " mgpos1: " << mgpos1 << "\n";
//    std::cerr << "    len: " << len << "\n";
//    std::cerr << "  gend1: " << gend1 << "\n";
//    std::cerr << " mgend1: " << mgend1 << "\n";
//    std::cerr << "      i: " << i << "\n";

    //existing reference
    for (uint32_t k=mgpos1; k<=std::min(mgend1,gend1); ++k)
    {
        ++P[i].N;

        uint8_t q = qual[spos0+j];
        if (q>snp_baseq_cutoff)
        {
            char alt = (bam_base2char(bam_seqi(seq, spos0+j)));

            if (alt!=P[i].R)
            {
                ++P[i].X[base2index(alt)];
                P[i].ALT_Q.push_back(q);
            }
            else
            {
                P[i].REF_Q.push_back(q);
            }
        }
        else
        {
            ++P[i].F;
        }

        i = inc(i);
        ++j;
    }

    //nonexisting reference
    if (mgend1>gend1)
    {
        uint32_t alen = len-j;
        char* ref = get_sequence(chrom, mgpos1+j, alen);
        uint32_t l = 0;

        for (uint32_t k=mgpos1+j; k<=mgend1; ++k)
        {
            P[i].R = ref[l];
            ++P[i].N;

            uint8_t q = qual[spos0+j];
            if (q>snp_baseq_cutoff)
            {
                char alt = (bam_base2char(bam_seqi(seq, spos0+j)));

                if (alt!=P[i].R)
                {
                    ++P[i].X[base2index(alt)];
                    P[i].ALT_Q.push_back(q);
                }
                else
                {
                    P[i].REF_Q.push_back(q);
                }
            }
            else
            {
                ++P[i].F;
            }

            inc_end0();
            i = inc(i);
            ++j;
            ++l;
        }

        free(ref);
    }
}

/**
 * Updates a stretch of deleted bases identified by D in the cigar string.
 */
void Pileup::add_D(uint32_t gpos1, uint32_t len)
{
    //there should never be a need to perform 3' padding for deletions
    //add_3prime_padding(gpos1);

    //check if the base exists.
    if (gpos1>get_gend1())
    {
        //change to a ignore, return?
        std::cerr << "anchor base not present for deletion update " << chrom << ":" << gpos1 << " (gpos1) >" << get_gend1() << " (gend1), gbeg1 = " << get_gbeg1() << "\n";
        abort();
    }

    char* ref = get_sequence(chrom, gpos1+1, len);
    std::string del(ref);
    free(ref);

    uint32_t i = g2i(gpos1);
    if (!is_normalized(P[i].R, del))
    {
        uint32_t a_gpos1 = gpos1;
        std::string a_ref(1, P[i].R);
        a_ref.append(del);
        std::string a_alt(1, P[i].R);
        normalize(chrom, a_gpos1, a_ref, a_alt);

        if (a_gpos1<gbeg1)
        {
            fprintf(stderr, "[%s:%d %s] deletion left aligned to beyond the bounds of the pileup: %s:%d<%d\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), gpos1, gbeg1);
            add_5prime_padding(a_gpos1);
        }

        if (debug>=2)  std::cerr << "\t\t\tdeletion left aligned : " << chrom << ":" << gpos1 << ":" << P[i].R << del << "/" << P[i].R  << " => " << chrom << ":" << a_gpos1 << ":" << a_ref << "/" << a_alt << "\n";
        uint32_t j = g2i(a_gpos1);
        ++P[j].D[a_ref.substr(1)];
    }
    else
    {
        ++P[i].D[del];
    }

    i = g2i(gpos1+1);
    for (uint32_t j = 0; j<del.size(); ++j)
    {
        P[i].R = del[j];
        if (i==end0) inc_end0();
        i = inc(i);
    }
}

/**
 * Updates an occurence of an insertion.
 */
void Pileup::add_I(uint32_t gpos1, std::string& ins, uint32_t rpos1)
{
    //there should never be a need to perform 3' padding for insertions
    //add_3prime_padding(gpos1);

    uint32_t i = g2i(gpos1);

    if (!is_normalized(P[i].R, ins))
    {
        uint32_t a_gpos1 = gpos1;
        std::string a_ref(1, P[i].R);
        std::string a_alt(1, P[i].R);
        a_alt.append(ins);
        normalize(chrom, a_gpos1, a_ref, a_alt);

        if (a_gpos1<gbeg1)
        {
            fprintf(stderr, "[%s:%d %s] insertion left aligned to beyond the bounds of the pileup: %s:%d<%d\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), a_gpos1, gbeg1);
            add_5prime_padding(a_gpos1);
        }

        if (debug>=2)  std::cerr << "\t\t\tinsertion left aligned : " << chrom << ":" << gpos1 << ":" << P[i].R << "/" << P[i].R << ins << " => " << chrom << ":" << a_gpos1 << ":" << a_ref << "/" << a_alt << "\n";
        uint32_t j = g2i(a_gpos1);
        ++P[j].I[a_alt.substr(1)];

        //for insertions shifted beyond the edge of a read alignment
        if (a_gpos1 < rpos1)
        {
            ++P[j].N;
        }
    }
    else
    {
        ++P[i].I[ins];
    }
}

/**
 * Updates the occurence of a left soft clip.
 */
void Pileup::add_LSC(uint32_t gpos1, std::string& alt, float mean_qual, char strand)
{
    add_3prime_padding(gpos1);

    uint32_t i = g2i(gpos1);
    SoftClipInfo& info = P[i].J[alt];
    ++info.no;
    info.mean_quals.push_back(mean_qual);
    info.strands.push_back(strand);
    if (i==end0)
    {
        P[i].R = get_base(chrom, gpos1);
        inc_end0();
    }
}

/**
 * Updates the occurence of a right soft clip.
 */
void Pileup::add_RSC(uint32_t gpos1, std::string& alt, float mean_qual, char strand)
{
    //there should never be a need to perform 3' padding for right soft clips
    //add_3prime_padding(gpos1);

    uint32_t i = g2i(gpos1);
    SoftClipInfo& info = P[i].K[alt];
    ++info.no;
    info.mean_quals.push_back(mean_qual);
    info.strands.push_back(strand);
    if (i==end0) inc_end0();
}

/**
 * Inserts a stretch of reference padding bases at the 5' prime end of the buffer from gpos1 if gpos1 is behind the start of the pileup.
 *
 * @gpos1 - 1 based genome position
 */
void Pileup::add_5prime_padding(uint32_t gpos1)
{
    if (is_empty())
    {
        return;
    }

    if (gpos1<gbeg1)
    {
        if (max_size()-size()<gbeg1-gpos1)
        {
//            std::cerr << "max size: " << max_size() << "\n";
//            std::cerr << "    size: " <<size() << "\n";
//            std::cerr << "   gbeg1: " <<gbeg1 << "\n";
//            std::cerr << "   gpos1: " <<gpos1 << "\n";
            fprintf(stderr, "[%s:%d %s] buffer overflow. requires %d but maximum size is %d\n", __FILE__, __LINE__, __FUNCTION__, size()+gbeg1-gpos1, max_size());
            abort();
        }

        uint32_t o_gbeg1 = gbeg1;

        beg0 = diff(beg0, gbeg1-gpos1);
        gbeg1 = gpos1;

        if (debug) std::cerr << "add_5prime_padding: (" << o_gbeg1 << "," << get_gend1() << ") extended to (" << gbeg1 << "," << get_gend1() << "\n";
  }
}

/**
 * Inserts a stretch of reference padding bases at the 3' prime end of the buffer to gpos1 if gpos1 is ahead of end of pileup.
 *
 * @gpos1 - 1 based genome position
 */
void Pileup::add_3prime_padding(uint32_t gpos1)
{
    if (is_empty())
    {
        return;
    }

    uint32_t cgend1 = get_gend1();
    if (gpos1>cgend1+1)
    {
        //std::cerr << "adding " << gpos1 << "-" << cgend1 << " " << (gpos1-cgend1-1) << "\n";
        end0 = inc(end0, gpos1-cgend1-1);
        if (debug) std::cerr << "add_3prime_padding: gpos1=" << cgend1 << ", gend1=" << get_gend1() << "\n";
    }
}

/**
 * Inserts a stretch of reference bases from read.
 *
 * requirement: stretch of bases on read is detected to be the same as the reference
 *              using CIGAR and MD tag.
 *              pileup is empty
 *              pileup is filled and partially overlaps with this stretch
 *
 * @gpos1 - starting genome position
 * @spos0 - starting position in the sequence
 * @len   - length of sequence to be inserted
 * @seq   - bam encoded sequence
 * @end   - if the last base is considered the end of a read alignment
 */
void Pileup::add_ref(uint32_t gpos1, uint32_t spos0, uint32_t len, uint8_t* seq)
{
    add_3prime_padding(gpos1);

    if (debug) std::cerr << "add_ref: gpos1=" << gpos1 << ", spos0=" << spos0 << ", len=" << len << "\n";

    uint32_t i = g2i(gpos1);

    for (uint32_t j = 0; j<len; ++j)
    {
        P[i].R = (bam_base2char(bam_seqi(seq, spos0+j)));
        ++P[i].N;
        if (i==end0) inc_end0();
        i = inc(i);
    }
}

/**
 * Updates an occurence of a SNP.
 */
void Pileup::add_snp(uint32_t gpos1, char ref, char alt, uint8_t qual, uint32_t baseq_cutoff)
{
    add_3prime_padding(gpos1);

    uint32_t i = g2i(gpos1);
    P[i].R = ref;

    ++P[i].N;

    if (qual>baseq_cutoff)
    {
        if (alt!=P[i].R)
        {
            ++P[i].X[base2index(alt)];
            P[i].ALT_Q.push_back(qual);
        }
        else
        {
            P[i].REF_Q.push_back(qual);
        }
    }
    else
    {
        ++P[i].F;
    }

    if (i==end0) inc_end0();
}

/**
 * Updates an occurence of a deletion.
 */
void Pileup::add_del(uint32_t gpos1, std::string& del)
{
    //there should never be a need to perform 3' padding for deletions
    //add_3prime_padding(gpos1);

    //std::cerr << "adding del " << gpos1 << " " << g2i(gpos1) << "\n";
    uint32_t i = g2i(gpos1);

    if (!is_normalized(P[i].R, del))
    {
        uint32_t a_gpos1 = gpos1;
        std::string a_ref(1, P[i].R);
        a_ref.append(del);
        std::string a_alt(1, P[i].R);
        normalize(chrom, a_gpos1, a_ref, a_alt);

        if (a_gpos1<gbeg1)
        {
            fprintf(stderr, "[%s:%d %s] deletion left aligned to beyond the bounds of the pileup: %s:%d<%d\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), a_gpos1, gbeg1);
            add_5prime_padding(a_gpos1);
        }

        if (debug>=2)  std::cerr << "\t\t\tdeletion left aligned : " << chrom << ":" << gpos1 << ":" << P[i].R << del << "/" << P[i].R << " => " << chrom << ":" << a_gpos1 << ":" << a_ref << "/" << a_alt << "\n";
        uint32_t j = g2i(a_gpos1);
        ++P[j].D[a_ref.substr(1)];
    }
    else
    {
        ++P[i].D[del];
    }

    //fill up for reference too.
    i = g2i(gpos1+1);
    for (uint32_t j = 0; j<del.size(); ++j)
    {
        //std::cerr << "adding to ref under del " << alt[j] << "\n";
        P[i].R = del[j];
        if (i==end0) inc_end0();
        i = inc(i);
    }
}

/**
 * Updates an occurence of an insertion.
 */
void Pileup::add_ins(uint32_t gpos1, std::string& ins, uint32_t rpos1)
{
    //there should never be a need to perform 3' padding for insertions
    //add_3prime_padding(gpos1);

    uint32_t i = g2i(gpos1);

    if (!is_normalized(P[i].R, ins))
    {
        uint32_t a_gpos1 = gpos1;
        std::string a_ref(1, P[i].R);
        std::string a_alt(1, P[i].R);
        a_alt.append(ins);
        normalize(chrom, a_gpos1, a_ref, a_alt);

        if (a_gpos1<gbeg1)
        {
            fprintf(stderr, "[%s:%d %s] insertion left aligned to beyond the bounds of the pileup: %s:%d<%d\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), a_gpos1, gbeg1);
            add_5prime_padding(a_gpos1);
        }

        if (debug>=2) std::cerr << "\t\t\tinsertion left aligned : " << chrom << ":" << gpos1 << ":" << P[i].R << "/" << P[i].R << ins << " => " << chrom << ":" << a_gpos1 << ":" << a_ref << "/" << a_alt << "\n";
        uint32_t j = g2i(a_gpos1);
        ++P[j].I[a_alt.substr(1)];

        //for insertions shifted beyond the edge of a read alignment
        if (a_gpos1 < rpos1)
        {
            ++P[j].N;
        }
    }
    else
    {
        ++P[i].I[ins];
    }
}

/**
 * Updates the occurence of a left soft clip.
 */
void Pileup::add_lsclip(uint32_t gpos1, std::string& alt, float mean_qual, char strand)
{
    add_3prime_padding(gpos1);

    uint32_t i = g2i(gpos1);
    SoftClipInfo& info = P[i].J[alt];
    ++info.no;
    info.mean_quals.push_back(mean_qual);
    info.strands.push_back(strand);
    if (i==end0) inc_end0();
}

/**
 * Updates the occurence of a right soft clip.
 */
void Pileup::add_rsclip(uint32_t gpos1, std::string& alt, float mean_qual, char strand)
{
    //add_3prime_padding(gpos1);

    uint32_t i = g2i(gpos1);
    SoftClipInfo& info = P[i].K[alt];
    ++info.no;
    info.mean_quals.push_back(mean_qual);
    info.strands.push_back(strand);
}

/**
 * Checks if an indel is normalized.
 */
bool Pileup::is_normalized(char ref, std::string& indel)
{
    return ref != indel.at(indel.size()-1);
}

/**
 * Normalize a biallelic variant.
 *
 * If N exists in either of the alleles, the normalization does not proceed.
 */
void Pileup::normalize(std::string& chrom, uint32_t& pos1, std::string& ref, std::string& alt)
{
    if (ref.find_first_of("Nn", 0)!=std::string::npos || alt.find_first_of("Nn", 0)!=std::string::npos) return;

    while (ref.at(ref.size()-1)==alt.at(alt.size()-1))
    {
        --pos1;
        char base = get_base(chrom, pos1);
        ref.insert(0, 1, base);
        alt.insert(0, 1, base);
        ref.erase(ref.size()-1, 1);
        alt.erase(alt.size()-1, 1);
    }
};

/**
 * Converts base to bam encoding.
 */
uint32_t Pileup::base2index(char base)
{
    switch(base)
    {
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 4;
        case 'T':
            return 8;
        case 'N':
            return 15;
        default:
            return 15;
    }
}

/**
 * Get a base.
 */
char Pileup::get_base(std::string& chrom, uint32_t& pos1)
{
    int ref_len = 0;
    char *refseq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos1-1, pos1-1, &ref_len);
    if (!refseq)
    {
        fprintf(stderr, "[%s:%d %s] failure to extrac base from fasta file: %s:%d: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos1-1);
        exit(1);
    }
    char base = refseq[0];
    free(refseq);

    return base;
};

/**
 * Get a sequence.  User have to free the char* returned.
 */
char* Pileup::get_sequence(std::string& chrom, uint32_t pos1, uint32_t len)
{
    int ref_len = 0;
    char* seq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos1-1, pos1+len-2, &ref_len);
    if (!seq || ref_len!=len)
    {
        fprintf(stderr, "[%s:%d %s] failure to extract sequence from fasta file: %s:%d: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos1-1);
        exit(1);
    }

    return seq;
};

/**
 * Print pileup state.
 */
void Pileup::print_state()
{
    std::cerr << "******************" << "\n";
    std::cerr << "gindex   : " << gbeg1 << "-" << get_gend1() << "\n";
    std::cerr << "index   : " << beg0 << "-" << end0 << " (" << size() << ")\n";
    std::cerr << "******************" << "\n";
    uint32_t k = 0;
    for (uint32_t i=beg0; i!=end0; i=inc(i))
    {
        P[i].print(gbeg1+k);
        ++k;
    }
    std::cerr << "******************" << "\n";
}

/**
 * Print pileup state extent.
 */
void Pileup::print_state_extent()
{
    std::cerr << "******************" << "\n";
    std::cerr << "gindex   : " << gbeg1 << "-" << get_gend1() << "\n";
    std::cerr << "index   : " << beg0 << "-" << end0 << " (" << size() << ")\n";
    std::cerr << "******************" << "\n";
}