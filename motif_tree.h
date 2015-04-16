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
#include <list>
#include <iostream>
#include <queue>
#include "candidate_motif.h"
#include "motif_map.h"

/**
 * Macros for accessing node.
 */

#define A 0
#define C 1
#define G 2
#define T 3

/**
 * struct for encoding a node.
 */
typedef struct 
{
    uint32_t cindex; //index of the cannonical form
    uint32_t count;  //count of occurence of this motif
    uint32_t len;    //length of the motif
    uint32_t seq;    //sequence encoded in 2 bits per base, this allows for a maximum length of 16 bases
    
} node;

#define index2base(i) ("ACGT"[(i)])

/**
 * Motif Suffix Tree for selecting candidate motifs.
 *
 * Each node is represented by int64_t
 *
 * the first 20 bits indexes to the cannonical form
 * the next 12 bits gives the count
 * the last 32 bits provides the sequence represented by 2 bits per base
 * the 
 *  
 *
 */
class MotifTree
{
    public:
    node* tree;
    MotifMap *mm;

    /**
     * Constructor.
     */
    MotifTree(uint32_t max_len);

    /**
     * Destructor.
     */
    ~MotifTree();

    /**
     * Clear the suffix tree.
     */
    void clear();

    /**
     * Construct suffix tree based on sequence.
     */
    void set_sequence(char* sequence);

    /**
     * Construct suffix tree based on sequence up to max_motif_len.
     */
    void set_sequence(char* sequence, int32_t max_motif_len);

    /**
     * Gets candidate motifs up to max_motif_len.
     */
    void get_candidate_motifs(std::vector<CandidateMotif>& candidate_motifs);
    
    /**
     * Gets index of child.
     */
    uint32_t get_first_child(uint32_t index);
    

    private:

    /**
     * Adds a suffix of sequence from start to end.
     */
    void add_suffix(char* sequence, int32_t start, int32_t end);

    /**
     * Converts base to index.
     */
    int32_t base2index(char base);

    /**
     * Print node.
     */
    void print_node(node* n);

};

#undef A
#undef C
#undef G
#undef T
#undef N

#endif