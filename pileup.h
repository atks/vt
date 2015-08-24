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

/**
 * Contains sufficient statistic for a position in the pileup.
 */
class SoftClipInfo
{
    public:
    uint32_t no;
    std::vector<float> mean_quals;
    std::vector<char> strands;

    /**
     * Constructor.
     */
    SoftClipInfo() { clear();};

    /**
     * Clears the soft clipped information.
     */
    void clear();
};

/**
 * Contains sufficient statistic for a position in the pileup.
 */
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
    std::map<std::string, SoftClipInfo> J;
    std::map<std::string, SoftClipInfo> K;
    //occurence of all observations for internal bases
    uint32_t N;
    //number of bases that fail quality cutoff
    uint32_t F;
    //occurence of all observations for bases on the ends of reads
    //this is important as indels have to be contained within a read
    //so if this is not distinguished, the relative proportion of
    //the reads containing the indel out of reads that could contain
    //that indel is will be underestimated.
    uint32_t E;
    //base qualities for reference allele
    std::vector<uint32_t> REF_Q;
    //base qualities for alternative allele
    std::vector<uint32_t> ALT_Q;

    //to count evidences
    //SNPs - X[A] / (N+E)
    //INDELs - #indel / N
    //SCLIPS -

    /**
     * Constructor.
     */
    PileupPosition();

    /**
     * Clears pileup position.
     */
    void clear();

    /**
     * Returns true if pileup record is cleared.
     */
    bool is_cleared();

    /**
     * Prints pileup position.
     */
    void print();

    /**
     * Prints pileup position.
     */
    void print(uint32_t gpos1);
};

/**
 * A pileup for variant discovery.
 * This only contains the pileup structure with
 * methods to access and modify the pileup.
 *
 * This pileup only works for a single chromosome
 * the user has to update the pileup and clear it
 * if the chromosome is changed.
 *
 * Transferring the pileup information to a VCF
 * file is delegated to a class that uses this
 * structure.
 */
class Pileup
{
    public:
    uint32_t buffer_size;
    uint32_t buffer_size_mask;
    uint32_t window_size;
    std::vector<PileupPosition> P;

    int32_t tid;
    uint32_t beg0, end0; // index of P for the start and end position in 0 base coordinates


    //  beg0 points to the beginning of the pileup
    //  end0 points to the position after the end of the pileup
    //
    //      end0
    //       |
    //       V
    //  +--+--+--+--+--+--+--+--+--+--+--+--+
    //  |  |  |  |  |  |  |  |  |  |  |  |  |
    //  +--+--+--+--+--+--+--+--+--+--+--+--+
    //       ^
    //       |
    //      beg0
    //
    //  when beg0==end0, pileup is empty
    //  gbeg1 - is the coordinate of the genome position that beg0 represents.
    //
    //  invariance:  pileup is always continuous
    //
    //  how do we check if
    //

    //for use if we need to update the reference.
    std::string chrom;

    //where should we check if a read is from another chromosome?

    //genome position at the beginning of the pileup
    uint32_t gbeg1;

    int32_t debug;

    faidx_t *fai;

    public:

    /**
     * Constructor.
     *
     * @k - size of pileup is 2^k
     */
    Pileup(uint32_t k=10, uint32_t window_size=256);

    /**
     * Overloads subscript operator for accessing pileup positions.
     */
    PileupPosition& operator[] (const int32_t i);

    /**
     * Returns the maximum size of the pileup.
     */
    uint32_t max_size();

    /**
     * Returns the size of the pileup.
     */
    uint32_t size();

    /**
     * Checks if buffer is empty
     */
    bool is_empty();

    /**
     * Set reference fasta file.
     */
    void set_reference(std::string& ref_fasta_file);

    /**
     * Set debug.
     */
    void set_debug(int32_t debug);

    /**
     * Sets tid.
     */
    void set_tid(uint32_t tid);

    /**
     * Gets tid.
     */
    uint32_t get_tid();

    /**
     * Sets chrom.
     */
    void set_chrom(std::string& chrom);

    /**
     * Gets chrom.
     */
    std::string get_chrom();

    /**
     * Gets window_size.
     */
    uint32_t get_window_size();

    /**
     * Sets gbeg1.
     */
    void set_gbeg1(uint32_t gbeg1);

    /**
     * Gets gbeg1.
     */
    uint32_t get_gbeg1();

    /**
     * Gets gend1.
     */
    uint32_t get_gend1();

    /**
     * Sets beg0.
     */
    void set_beg0(uint32_t beg0);

    /**
     * Sets end0.
     */
    void set_end0(uint32_t end0);

    /**
     * Gets the index of the first element.
     */
    uint32_t begin();

    /**
     * Gets the index of the last element.
     */
    uint32_t end();

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
     * Returns the difference between 2 buffer positions
     */
    uint32_t diff(uint32_t i, uint32_t j);

    /**
     * Updates a stretch of aligned bases identified by M in the cigar string.
     */
    void add_M(uint32_t gpos1, uint32_t spos0, uint32_t len, uint8_t* seq, uint8_t* qual, uint32_t snp_baseq_cutoff);

    /**
     * Updates a stretch of deleted bases identified by D in the cigar string.
     */
    void add_D(uint32_t gpos1, uint32_t len);

    /**
     * Updates an occurence of an insertion.
     */
    void add_I(uint32_t gpos1, std::string& ins, uint32_t rpos1);

    /**
     * Updates the occurence of a left soft clip.
     */
    void add_LSC(uint32_t gpos1, std::string& alt, float mean_qual, char strand);

    /**
     * Updates the occurence of a right soft clip.
     */
    void add_RSC(uint32_t gpos1, std::string& alt, float mean_qual, char strand);

    /**
     * Inserts a stretch of reference padding bases at the 5' prime end of the buffer from gpos1 if gpos1 is behind the start of the pileup.
     *
     * @gpos1 - 1 based genome position
     */
    void add_5prime_padding(uint32_t gpos1);

    /**
     * Inserts a stretch of reference padding bases at the 3' prime end of the buffer to gpos1 if gpos1 is ahead of end of pileup.
     *
     * @gpos1 - 1 based genome position
     */
    void add_3prime_padding(uint32_t gpos1);

    /**
     * Updates a stretch of reference bases.
     */
    void add_ref(uint32_t gpos1, uint32_t spos0, uint32_t len, uint8_t* seq);

    /**
     * Updates an occurence of a SNP.
     */
    void add_snp(uint32_t gpos1, char ref, char alt, uint8_t qual, uint32_t baseq_cutoff);

    /**
     * Updates an occurence of a deletion.
     */
    void add_del(uint32_t gpos1, std::string& del);

    /**
     * Updates an occurence of an insertion.
     */
    void add_ins(uint32_t gpos1, std::string& ins, uint32_t rpos1);

    /**
     * Updates the occurence of a left soft clip.
     */
    void add_lsclip(uint32_t gpos1, std::string& alt, float mean_qual, char strand);

    /**
     * Updates the occurence of a right soft clip.
     */
    void add_rsclip(uint32_t gpos1, std::string& alt, float mean_qual, char strand);

    /**
     * Updates the last aligned base in a read.
     *
     * @gpos1 - starting genome position
     */
    void update_read_end(uint32_t gpos1);

    /**
     * Checks if a variant is normalized.
     */
    bool is_normalized(char ref, std::string& indel);

    /**
     * Normalize a biallelic variant.
     *
     * If N exists in either of the alleles, the normalization does not proceed.
     */
    void normalize(std::string& chrom, uint32_t& pos1, std::string& ref, std::string& alt);

    /**
     * Converts base to bam encoding.
     */
    uint32_t base2index(char base);

    /**
     * Get a base.
     */
    char get_base(std::string& chrom, uint32_t& pos1);

    /**
     * Get a sequence.  User have to free the char* returned.
     */
    char* get_sequence(std::string& chrom, uint32_t pos1, uint32_t len);

    /**
     * Print pileup state.
     */
    void print_state();

    /**
     * Print pileup state extent.
     */
    void print_state_extent();
};

#endif