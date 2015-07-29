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

#ifndef MOTIF_TREE_H
#define MOTIF_TREE_H

#include <cstdint>
#include <cstring>
#include <vector>
#include <iostream>
#include <queue>
#include <map>
#include "motif_map.h"

/**
 * Macros for accessing node.
 */

#define A 0
#define C 1
#define G 2
#define T 3

/**
 * struct for storing sequence content.
 */
typedef struct
{
    uint32_t base[4]; //counts of bases
    uint32_t n; //total number of bases
} scontent;

/**
 * struct for encoding a node.
 */
typedef struct
{
    uint32_t index; //index
    uint32_t cindex; //index of the cannonical form
    uint32_t count;  //count of occurence of this motif
    uint32_t len;    //length of the motif
    uint32_t seq;    //sequence encoded in 2 bits per base, this allows for a maximum length of 16 bases

} node;

#define index2base(i) ("ACGT"[(i)])

/**
 * Candidate Motif.
 */
class CandidateMotif
{
    public:
    std::string motif;
    float score;
    uint32_t len;
    float fit;

    CandidateMotif(std::string motif, float score, uint32_t len, float fit)
    {
        this->motif = motif;
        this->score = score;
        this->len = len;
        this->fit = fit;
    }
};

/**
 * Comparator for Candidate Motif for use in priority_queue.
 */
class CompareCandidateMotif
{
    public:
    bool operator()(CandidateMotif& a, CandidateMotif& b)
    {
        if (a.score!=b.score)
        {
            return a.score < b.score;
        }
        else
        {
            return a.len > b.len;
        }
    }
};

/**
 * Motif Suffix Tree for selecting candidate motifs.
 *
 * Each node is represented by int64_t
 *
 * the first 20 bits indexes to the cannonical form
 * the next 12 bits gives the count
 * the last 32 bits provides the sequence represented by 2 bits per base
 */
class MotifTree
{
    public:
    node* tree;
    MotifMap *mm;
    uint32_t max_len;
    std::map<uint32_t, uint32_t> cm;
    std::vector<uint32_t> lc;   // for counting the number of motifs of length x.
    std::priority_queue<CandidateMotif, std::vector<CandidateMotif>, CompareCandidateMotif > pcm;
    uint32_t cmax_len; //candidate maximum length
    bool debug;

    /**
     * Constructor.
     */
    MotifTree(uint32_t max_len, bool debug=false);

    /**
     * Destructor.
     */
    ~MotifTree();

    /**
     * Clear the suffix tree.
     */
    void clear();

    /**
     * Gets subsequence from a C string.
     */
    uint32_t get_sub_seq(char* seq, uint32_t len, uint32_t spos0, uint32_t& s);

    /**
     * Inserts prefix s into tree.
     */
    void insert_prefix(uint32_t s, uint32_t len);

    /**
     * Consolidate motif counts.
     */
    void consolidate_motif_counts();

    /**
     * Consolidate motif counts.
     */
    void consolidate_motif_counts(node* n);

    /**
     * Construct suffix tree based on sequence up to max_motif_len.
     */
    void set_sequence(char* seq, uint32_t len);

    /**
     * Shifts a string.
     */
    std::string shift_str(std::string& seq, uint32_t i);

    /**
     * Checks if two copies of a motif exists in a seq.
     */
    bool exist_two_copies(std::string& seq, std::string& motif);

    /**
     * Compute fit of expected content of nucleotides.
     */
    float compute_fit(uint32_t index, scontent* sc);

    /**
     * Detects candidate motifs from seq of length len.
     */
    void detect_candidate_motifs(std::string& seq);

    /**
     * Detects candidate motifs from seq of length len.
     */
    void detect_candidate_motifs(char* seq, uint32_t len);

    /**
     * Gets index of first child.
     */
    uint32_t get_first_child(uint32_t seq, uint32_t len);


    private:

    /**
     * Converts base to index.
     */
    int32_t base2index(char base);

    /**
     * Print node.
     */
    void print_node(node* n);

    /**
     * Print tree.
     */
    void print_tree();

    /**
     * Print tree.
     */
    void print_tree(node* n);

};

#undef A
#undef C
#undef G
#undef T
#undef N

#endif