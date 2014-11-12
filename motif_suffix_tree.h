/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#ifndef MOTIF_SUFFIX_TREE_H
#define MOTIF_SUFFIX_TREE_H

#include <cstdint>
#include <cstring>
#include <vector>
#include <list>
#include <iostream>
#include <queue>
#include "candidate_motif.h"

#define A 0
#define C 1
#define G 2
#define T 3

/**
 * A node for the Motif Suffix Tree.
 */
class MotifSuffixTreeNode
{
    public:
    MotifSuffixTreeNode* parent;
    MotifSuffixTreeNode* children[4];
    int32_t n_children[4];
    std::string suffix;
    int32_t count;

    /**
     * Constructs a MotifSuffixTreeNode.
     */
    MotifSuffixTreeNode();

    /**
     * Constructs a MotifSuffixTreeNode and initialize it with a parent and suffix.
     */
    MotifSuffixTreeNode(MotifSuffixTreeNode* parent, std::string& suffix);

    /**
     * Clear the suffix tree node.
     */
    void clear();

    /**
     * Constructs an MotifSuffixTreeNode.
     */
    ~MotifSuffixTreeNode();

    private:
};

/**
 * Motif Suffix Tree for selecting candidate motifs.
 */
class MotifSuffixTree
{
    public:
    MotifSuffixTreeNode* root;
            
    /**
     * Constructor.
     */
    MotifSuffixTree();

    /**
     * Destructor.
     */
    ~MotifSuffixTree();

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
        
    private:
        
    /**
     * Adds a suffix of sequence from start to end.
     */
    void add_suffix(char* sequence, int32_t start, int32_t end);
    
    /**
     * Converts base to index.
     */
    int32_t base2index(char base);
};

#undef A
#undef C
#undef G
#undef T

#endif