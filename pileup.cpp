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
 * Clears pileup position.
 */
void PileupPosition::clear()
{
    R = 'N';
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
    std::cerr << "ref: " << R << "\n";
    std::cerr << "A: " << X[1] << "\n";
    std::cerr << "C: " << X[2] << "\n";
    std::cerr << "G: " << X[4] << "\n";
    std::cerr << "T: " << X[8] << "\n";
    std::cerr << "N: " << X[15] << "\n";
    for (std::map<std::string, uint32_t>::iterator i = D.begin(); i!=D.end(); ++i)
    {
        std::cerr << i->first << ": " << i->second << "\n";
    }
    for (std::map<std::string, uint32_t>::iterator i = I.begin(); i!=D.end(); ++i)
    {
        std::cerr << i->first << ": " << i->second << "\n";
    }
    for (std::map<std::string, uint32_t>::iterator i = J.begin(); i!=D.end(); ++i)
    {
        std::cerr << i->first << ": " << i->second << "\n";
    }
    for (std::map<std::string, uint32_t>::iterator i = K.begin(); i!=D.end(); ++i)
    {
        std::cerr << i->first << ": " << i->second << "\n";
    }
    std::cerr << "#evidences: " << N << "\n";
    std::cerr << "#tail evidences: " << E << "\n";
}

/**
 * Constructor.
 */
Pileup::Pileup(uint32_t k)
{
    //Buffer size must be a power of 2^k.
    buffer_size = 1 << k;
    buffer_size_mask = 0xFFFF >> (32-k);

    P.resize(buffer_size);

    beg0 = 0;
    end0 = 0;

    gbeg1 = 0;
    gend1 = 0;
};

/**
 * Converts gpos1 to index in P.
 */
size_t Pileup::g2i(uint32_t gpos1)
{
    return add(beg0, gpos1-gbeg1); 
}

/**
 * Checks if the position is present.
 */
bool Pileup::position_is_present(int32_t tid, uint32_t gpos1)
{
    return (tid==this->tid) && (gpos1>=this->gbeg1) && (gpos1<=this->gend1);
}

/**
 * Inserts a stretch of reference bases.
 */
void Pileup::add_ref(uint32_t gpos1, uint32_t spos0, uint32_t len, uint8_t* seq, bool end)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        gend1 = gpos1;
    }
    
    for (uint32_t i=0; i<len; ++i)
    {
        if (gpos1+i>gend1)
        {
           P[g2i(gpos1+i)].R = (bam_base2char(bam_seqi(seq, spos0+len))); 
           inc();        
        }
        
        ++P[g2i(gpos1+i)].N;
    }        

}

/**
 * Updates an occurence of a SNP.
 */
void Pileup::add_snp(uint32_t gpos1, char ref, char alt, bool end)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        gend1 = gpos1;
    }
    
    P[g2i(gpos1)].R = ref; 
    ++P[g2i(gpos1)].X[alt];                 
}

/**
 * Updates an occurence of a deletion.
 */
void Pileup::add_del(uint32_t gpos1, std::string& alt)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        gend1 = gpos1;
    }
    
    ++P[g2i(gpos1)].D[alt];  
}

/**
 * Updates an occurence of an insertion.
 */
void Pileup::add_ins(uint32_t gpos1, std::string& alt)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        gend1 = gpos1;
    }
    
    ++P[g2i(gpos1)].I[alt];  
}

/**
 * Inserts a reference base at pos0 into the buffer.
 */
void Pileup::add_lsclip(uint32_t gpos1, std::string& alt)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        gend1 = gpos1;
    }
    
    ++P[g2i(gpos1)].J[alt];  
}

/**
 * Inserts a reference base at pos0 into the buffer.
 */
void Pileup::add_rsclip(uint32_t gpos1, std::string& alt)
{
    if (is_empty())
    {
        gbeg1 = gpos1;
        gend1 = gpos1;
    }
    
    ++P[g2i(gpos1)].K[alt];  
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
inline size_t Pileup::diff(size_t i, size_t j)
{
    return (i>=j ? i-j : buffer_size-(j-i));
};

/**
 * Checks if buffer is empty.
 */
inline bool Pileup::is_empty()
{
    return beg0==end0;
};

/**
 *Increments buffer index i by 1.
 */
void Pileup::inc()
{
    ++gend1;
    end0 = (end0+1) & buffer_size_mask;
};

/**
 * Increments buffer index i by j.
 */
inline uint32_t Pileup::add(uint32_t i, uint32_t j)
{
    return (i+j) & buffer_size_mask;
};

/**
 * Decrements buffer index i by j.
 */
size_t Pileup::minus(size_t& i, size_t j)
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
void Pileup::minus(size_t& i)
{
    if (i>=buffer_size)
    {
        std::cerr << "Unaccepted buffer index: " << i << " ("  << buffer_size << ")\n";
        exit(1);
    }

    i = (i>=1 ? i-1 : buffer_size-1);
};

/**
 * Returns the size of the pileup.
 */
inline uint32_t Pileup::size()
{
    return (end0>=beg0 ? end0-beg0 : buffer_size-(beg0-end0));
}

/**
 * Gets the position in the buffer that corresponds to
 * the genome position indicated by pos.
 */
size_t Pileup::get_cur_pos0(size_t gpos1)
{
    //when buffer is empty
    if (is_empty())
    {
        //start_genome_pos0 = genome_pos0;
        return beg0;
    }
    else
    {
        if (gpos1-gbeg1-1>buffer_size)
        {
            std::cerr << "overflow buffer\n" ;
            //should allow for unbuffering here
        }

        return (beg0 + (gpos1-gbeg1-1))%buffer_size;
    }
};

/**
 * Print buffer contents for debugging purpose
 */
void Pileup::printBuffer()
{
    std::cout << "PRINT BUFFER" << "\n";
};

/**
 * Checks if a variant is normalized.
 */
bool Pileup::is_biallelic_normalized(std::string& ref, std::string& alt)
{
    bool last_base_same = ref.at(ref.size()-1) == alt.at(alt.size()-1);
    bool exists_len_one_allele = ref.size()==1 || alt.size()==1;
    bool first_base_same = ref.at(0) == alt.at(0);

    if (last_base_same || (!exists_len_one_allele && first_base_same))
    {
        return false;
    }
    else
    {
        return true;
    }
}

/**
 * Normalize a biallelic variant.
 */
void Pileup::normalize_biallelic(size_t pos0, std::string& ref, std::string& alt)
{
    if (!is_biallelic_normalized(ref, alt))
    {
        bool to_right_trim = true;
        bool to_left_extend = false;

        while (to_right_trim || to_left_extend)
        {
            //checks if right trimmable or left extendable
            if (!ref.empty() && !alt.empty())
            {
                if (ref.at(ref.size()-1) == alt.at(alt.size()-1))
                {
                    to_right_trim = true;
                    to_left_extend = false;
                }

                if (pos0==0 && (ref.size()==1||alt.size()==1))
                {
                    to_right_trim = false;
                    to_left_extend = false;
                }
            }
            else
            {
                to_right_trim = false;
                to_left_extend = true;
            }

            if (to_right_trim)
            {
                ref.erase(ref.size()-1);
                alt.erase(alt.size()-1);
            }

            if (to_left_extend)
            {
                --pos0;
                int ref_len = 0;
                char *refseq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos0, pos0, &ref_len);
                if (!refseq)
                {
                    fprintf(stderr, "[%s:%d %s] failure to extrac base from fasta file: %s:%zu: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos0);
                    exit(1);
                }
                char base = refseq[0];
                free(refseq);

                ref.insert(0, 1, base);
                alt.insert(0, 1, base);
            }
        }

        bool to_left_trim =  true;

        while (to_left_trim)
        {
            //checks if left trimmable.
            if (ref.size()==1 || ref.at(0)!=alt.at(0))
            {
                to_left_trim = false;
            }

            if (to_left_trim)
            {
                ref.erase(0, 1);
                alt.erase(0, 1);
                ++pos0;
            }
        }
    }
};