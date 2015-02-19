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

#ifndef PILEUP_H
#define PILEUP_H

#include <map>
#include <vector>
#include "utils.h"
#include "hts_utils.h"
#include "variant_manip.h"
#include "htslib/faidx.h"

class PileupPosition
{
    public:
    //reference base
    char R;
    //alternative bases
    uint32_t X[16];
    //for deletions and insertions with anchor R
    std::map<std::string, uint32_t> D;
    std::map<std::string, uint32_t> I;
    std::map<std::string, uint32_t> J;
    std::map<std::string, uint32_t> K;
    //occurence of all observations for internal bases
    uint32_t N;
    //occurence of all observations for bases on the ends of reads
    uint32_t E;

    //to count evidences
    //SNPs - X[A] / (N+E)
    //INDELs - #indel / N
    //SCLIPS -


    /**
     * Constructor.
     * Buffer size must be a power of 2.
     */
    PileupPosition(){};

    /**
     * Clears pileup position.
     */
    void clear();

    /**
     * Prints pileup position.
     */
    void print();
};

/**
 * A pileup for variant discovery.
 */
class Pileup
{
    private:
    uint32_t buffer_size;
    uint32_t buffer_size_mask;
    std::vector<PileupPosition> P;

    std::string chrom;
    int32_t tid;
    uint32_t beg0, end0; // index of P for the start and end position in 0 base coordinates

    //genome position at the beginning of the pileup
    uint32_t gbeg1, gend1;

    bool debug;

    faidx_t *fai;

    public:

    /**
     * Constructor.
     */
    Pileup(uint32_t k=10);

    /**
     * Overloads subscript operator for accessing pileup positions.
     */
    PileupPosition& operator[] (const int32_t i);

    /**
     * Check if flushable.
     */
    bool flushable(int32_t tid, uint32_t gpos1);

    /**
     * Gets gbeg1.
     */
    uint32_t get_gbeg1();

    /**
     * Gets gend1.
     */
    uint32_t get_gend1();
    
    /**
     * Sets gbeg1.
     */
    void set_gbeg1(uint32_t gbeg1);

    /**
     * Sets gend1.
     */
    void set_gend1(uint32_t gend1);

    /**
     * Gets the index of the first element.
     */
    uint32_t begin();

    /**
     * Gets the index of the last element.
     */
    uint32_t end();

    /**
     * Sets beg0.
     */
    void set_beg0(uint32_t gbeg0);

    /**
     * Sets end0.
     */
    void set_end0(uint32_t end0);
    
    /**
     * Inserts a stretch of reference bases.
     */
    void add_ref(uint32_t gpos1, uint32_t spos0, uint32_t len, uint8_t* seq, bool end);

    /**
     * Updates an occurence of a SNP.
     */
    void add_snp(uint32_t gpos1, char ref, char alt, bool end);

    /**
     * Updates an occurence of a deletion.
     */
    void add_del(uint32_t gpos1, std::string& alt);

    /**
     * Updates an occurence of an insertion.
     */
    void add_ins(uint32_t gpos1, std::string& alt);

    /**
     * Inserts a reference base at pos0 into the buffer.
     */
    void add_lsclip(uint32_t gpos1, std::string& alt);

    /**
     * Inserts a reference base at pos0 into the buffer.
     */
    void add_rsclip(uint32_t gpos1, std::string& alt);

    /**
     * Returns the size of the pileup.
     */
    uint32_t size();

    /**
     * Checks if buffer is empty
     */
    bool is_empty();

    /**
     * Increments i by 1 circularly.
     */
    uint32_t inc(uint32_t i);

    /**
     * Increments beg0 by 1.
     */
    void inc_beg0();

    /**
     * Increments end0 by 1.
     */
    void inc_end0();

    /**
     * Increments index i by j cyclically.
     */
    uint32_t inc(uint32_t i, uint32_t j);

    /**
     * Converts gpos1 to index in P.
     */
    uint32_t g2i(uint32_t gpos1);

    /**
     * Print pileup state.
     */
    void print_state();
        
    private:

    /**
     * Returns the difference between 2 buffer positions
     */
    size_t diff(size_t i, size_t j);

    /**
     * Print buffer contents for debugging purpose
     */
    void printBuffer();
        
    /**
     * Checks if a variant is normalized.
     */
    bool is_biallelic_normalized(std::string& ref, std::string& alt);

    /**
     * Normalize a biallelic variant.
     */
    void normalize_biallelic(size_t pos0, std::string& ref, std::string& alt);

};

#endif