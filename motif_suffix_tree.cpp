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

#include "motif_suffix_tree.h"

#define A 0
#define C 1
#define G 2
#define T 3
#define N 4

/**
 * Constructs a MotifSuffixTreeNode.
 */
MotifSuffixTreeNode::MotifSuffixTreeNode()
{
    parent = NULL;
    for (size_t b = A; b<=T; ++b)
    {
        children[b] = NULL;
    }
    this->suffix = "";
    count = 0;
};

/**
 * Constructs a MotifSuffixTreeNode and initialize it with a parent and a suffix.
 */
MotifSuffixTreeNode::MotifSuffixTreeNode(MotifSuffixTreeNode* parent, std::string& suffix)
{
    parent = parent;
    for (size_t b = A; b<=T; ++b)
    {
        children[b] = NULL;
    }
    this->suffix = suffix;
    count = 0;
};


/**
 * Clear the suffix tree node.
 */
void MotifSuffixTreeNode::clear()
{
    for (size_t b = A; b<=T; ++b)
    {
        if (children[b])
        {    
            children[b]->clear();
        }
    }
    count = 0;
};

/**
 * Constructs an MotifSuffixTreeNode.
 */
MotifSuffixTreeNode::~MotifSuffixTreeNode()
{
    for (size_t b = A; A<=T; ++b)
    {
        if (children[b])
        {
            delete children[b];
        }

        children[b] = NULL;
    }

    parent = NULL;
};

/**
 * Constructor.
 */
MotifSuffixTree::MotifSuffixTree()
{
    root = new MotifSuffixTreeNode();
};

/**
 * Destructor.
 */
MotifSuffixTree::~MotifSuffixTree()
{
    if (root) delete root;
    root = NULL;
};

/**
 * Construct suffix tree based on sequence.
 */
void MotifSuffixTree::set_sequence(char* sequence)
{
    set_sequence(sequence, strlen(sequence));
};

/**
 * Construct suffix tree based on sequence up to max_motif_len.
 */
void MotifSuffixTree::set_sequence(char* sequence, int32_t max_motif_len)
{
    size_t len = strlen(sequence);

    //i is starting
    for (size_t i=0; i<len; ++i)
    {
        add_suffix(sequence, i, i+max_motif_len-1);
    }
};

/**
 * Gets candidate motifs up to max_motif_len.
 */
void MotifSuffixTree::get_candidate_motifs(std::priority_queue<CandidateMotif, std::vector<CandidateMotif *>, CompareCandidateMotif>& candidate_motifs)
{
    //travel through tree
    MotifSuffixTreeNode* node = root;
     
     
   // std::priority_queue<CandidateMotif, std::vector<CandidateMotif *>, CompareCandidateMotif> pq;
     
     
     
    
};

/**
 * Adds a suffix of sequence from start to end.
 */
void MotifSuffixTree::add_suffix(char* sequence, int32_t start, int32_t end)
{
    MotifSuffixTreeNode* node = root;

    for (size_t i = start; i<=end; ++i)
    {
        int32_t base = base2index(sequence[i]);

        if (base==N) break;

        if (node->children[base]==NULL)
        {
            node->children[base] = new MotifSuffixTreeNode(node, node->suffix.append(1,base));
        }

        node = node->children[base];
        ++node->count;
    }
};

/**
 * Converts base to index.
 */
int32_t MotifSuffixTree::base2index(char base)
{
    switch (base)
    {
        case 'A':
            return A;
            break;
        case 'C':
            return C;
            break;
        case 'G':
            return G;
            break;
        case 'T':
            return T;
            break;
        default:
            return N;
    }
};

#undef A
#undef C
#undef G
#undef T
#undef N
