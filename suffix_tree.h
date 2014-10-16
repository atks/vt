/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H

#include <cstdint>
#include <vector>
#include <iostream>

#define A 0
#define C 1
#define G 2
#define T 3

class CandidateMotif
{
    public:
    std::string motif;
    float score;
};

class SuffixTreeNode
{
    public:
    SuffixTreeNode* parent;
    SuffixTreeNode* children[4];
    SuffixTreeNode* n_children[4];
    std::string suffix;
    int32_t count;

    /**
     * Constructs an SuffixTreeNode and initialize it with an interval.
     */
    SuffixTreeNode();

    /**
     * Constructs an SuffixTreeNode.
     */
    ~SuffixTreeNode();
    
    private:
};

class SuffixTree
{
    public:
    SuffixTreeNode* root;

    /**
     * Constructor.
     */
    SuffixTree();

    /**
     * Destructor.
     */
    ~SuffixTree();

    /**
     * Construct suffix tree based on sequence.
     */
    void set_sequence(std::string& sequence);

    /**
     * Construct suffix tree based on sequence up to max_motif_len.
     */
    void construct_tree(std::string& sequence, int32_t max_motif_len);

    /**
     * Gets candidate motifs up to max_motif_len.
     */
    void get_candidate_motifs(std::vector<CandidateMotif*>& candidate_motifs, int32_t max_motif_len);
    
    private:

    /**
     * Adds a suffix to the tree.
     */
    void add_suffix(std::string seq, int32_t start);   
};

#endif