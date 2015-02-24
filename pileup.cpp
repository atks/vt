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
    E = 0;
};

/**
 * Prints pileup position.
 */
void PileupPosition::print()
{
    std::cerr << R << ":";
    std::cerr << "ACGTN " << X[1] << "|" << X[2] << "|" << X[4] << "|" << X[8] << "|"<< X[15] << "\n";
    std::cerr << "N :" << N << ":";
    std::cerr << "E :" << E << "\n:";

    if (D.size()!=0)
    {
        for (std::map<std::string, uint32_t>::iterator i = D.begin(); i!=D.end(); ++i)
        {
            std::cerr << i->first << " (" << i->second << ")";
        }

        std::cerr << "\n";
    }
    if (I.size()!=0)
    {
        std::cerr << "INS: ";
        for (std::map<std::string, uint32_t>::iterator i = I.begin(); i!=I.end(); ++i)
        {
            std::cerr << i->first << " (" << i->second << ")";
        }

        std::cerr << "\n";
    }
    if (J.size()!=0)
    {
        std::cerr << "R: ";
        for (std::map<std::string, SoftClipInfo>::iterator i = J.begin(); i!=J.end(); ++i)
        {
            std::cerr << i->first << " (" << i->second.no << ")";
        }

        std::cerr << "\n";
    }
    if (K.size()!=0)
    {
        std::cerr << "L: ";
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
    std::cerr << gpos1 << ":" << R << ":" << N << "+" << E << "\n";

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
 * Constructor.
 *
 * @k - size of pileup is 2^k
 */
Pileup::Pileup(uint32_t k)
{
    //Buffer size is a power of 2^k.
    buffer_size = 1 << k;
    //this provides a cheaper way to do modulo operations for a circular array.
    buffer_size_mask = (0xFFFFFFFF >> (32-k));

    //std::cerr << buffer_size << ":" << buffer_size_mask  << "\n";
    P.resize(buffer_size);

    tid = -1;
    beg0 = end0 = 0;
    gbeg1 = 0;

    debug = false;
//    std::cerr << gbeg1 << "\n";
//    std::cerr << beg0 << "-" << end0 << "\n";
};

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
 * Sets tid.
 */
void Pileup::set_tid(uint32_t tid)
{
    this->tid = tid;
}

/**
 * Sets chrom.
 */
void Pileup::set_chrom(const char* chrom)
{
    this->chrom.assign(chrom);
}

/**
 * Gets chrom.
 */
std::string Pileup::get_chrom()
{
    return chrom;
}

/**
 * Converts gpos1 to index in P.
 * should this be responisble for updating values????!?!?!?!?!?!??!?!?!
 */
uint32_t Pileup::g2i(uint32_t gpos1)
{
    //not empty
    if (beg0!=end0)
    {
        //assumes connectedness!
        return (beg0 + (gpos1-gbeg1)) & buffer_size_mask;
    }
    else
    {
        //disconnected
        gbeg1 = gpos1;
        end0 = inc(beg0);
        return beg0;
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
    if (end0==beg0)
    {
        return 0;
    }    
    else
    {    
        return gbeg1 + diff(end0, beg0) -1;
    }
}

/**
 * Sets gbeg1.
 */
void Pileup::set_gbeg1(uint32_t gbeg1)
{
    this->gbeg1 = gbeg1;
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
 * Updates the last aligned base in a read.
 *
 * @gpos1 - starting genome position
 */
void Pileup::update_read_end(uint32_t gpos1)
{
    if (debug) std::cerr << "update_read_end_base: gpos1" << gpos1 << "\n";
    uint32_t i = g2i(gpos1);
    --P[i].N;
    ++P[i].E;
}

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
 * Updates an occurence of a SNP.
 */
void Pileup::add_snp(uint32_t gpos1, char ref, char alt)
{
    uint32_t i = g2i(gpos1);
    //std::cerr << "i in adding SNP " << i << " " << end0 << "\n";
    P[i].R = ref;
    ++P[i].X[base2index(alt)];
    ++P[i].N;
    if (i==end0) inc_end0(); 
}

/**
 * Updates an occurence of a deletion.
 */
void Pileup::add_del(uint32_t gpos1, std::string& del)
{
    //std::cerr << "adding del " << gpos1 << " " << g2i(gpos1) << "\n";
    
    uint32_t i = g2i(gpos1);
    
    if (!is_normalized(P[i].R, del))
    {
        char ref = P[i].R;
        std::string d = del;
        std::cerr << "not normalized " << chrom << ":" << gpos1 << ":" << ref  << d << ":" <<  d << "\n";
        normalize(chrom, gpos1, ref, d);
        std::cerr << "normalized " << chrom << ":" << gpos1 << ":" << ref  << d << ":" <<  d << "\n";
    }
    
    ++P[i].D[del];
    
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
void Pileup::add_ins(uint32_t gpos1, std::string& ins)
{
    uint32_t i = g2i(gpos1);
 
    if (!is_normalized(P[i].R, ins))
    {
        char ref = P[i].R;
        std::string d = ins;
        std::cerr << "not normalized " << chrom << ":" << gpos1 << ":" << ref  << d << ":" <<  d << "\n";
        normalize(chrom, gpos1, ref, d);
        std::cerr << "normalized " << chrom << ":" << gpos1 << ":" << ref  << d << ":" <<  d << "\n";
    }
    
    ++P[i].I[ins];
}

/**
 * Inserts a reference base at pos0 into the buffer.
 */
void Pileup::add_lsclip(uint32_t gpos1, std::string& alt, float mean_qual, char strand)
{
    uint32_t i = g2i(gpos1);
    SoftClipInfo& info = P[i].J[alt];
    ++info.no;
    info.mean_quals.push_back(mean_qual);
    info.strands.push_back(strand);    
    if (i==end0) inc_end0();
}

/**
 * Inserts a reference base at pos0 into the buffer.
 */
void Pileup::add_rsclip(uint32_t gpos1, std::string& alt, float mean_qual, char strand)
{
    uint32_t i = g2i(gpos1);
    SoftClipInfo& info = P[i].K[alt];
    ++info.no;
    info.mean_quals.push_back(mean_qual);
    info.strands.push_back(strand);  
}

/**
 * Overloads subscript operator for accessing pileup positions.
 */
PileupPosition& Pileup::operator[] (const int32_t i)
{
    return P[i];
}

/**
 * Returns the difference between 2 buffer positions
 */
inline uint32_t Pileup::diff(uint32_t i, uint32_t j)
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
   

/**
 * Checks if an indel is normalized.
 */
bool Pileup::is_normalized(char ref, std::string& indel)
{
    return ref == indel.at(indel.size()-1);
}

/**
 * Normalize a biallelic variant.
 */
void Pileup::normalize(std::string& chrom, uint32_t& pos1, char& ref, std::string& indel)
{
    while (indel.at(indel.size()-1)==ref)
    {
        --pos1;
        int ref_len = 0;
        char *refseq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos1-1, pos1-1, &ref_len);
        if (!refseq)
        {
            fprintf(stderr, "[%s:%d %s] failure to extrac base from fasta file: %s:%d: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos1-1);
            exit(1);
        }
        char base = refseq[0];
        free(refseq);
    
        ref = base;
        indel.insert(0, 1, base);
        indel.erase(indel.size()-1, 1);
    }
};