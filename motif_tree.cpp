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

#include "motif_tree.h"

#define A 0
#define C 1
#define G 2
#define T 3

/**
 * Constructor.
 */
MotifTree::MotifTree(uint32_t max_len)
{   
    mm = new MotifMap(max_len); 

    tree = (node *) malloc(sizeof(node)*mm->max_index+1);

        
    //perform mapping.
    for (uint32_t len=1; len<=max_len; ++len)
    {
        for (uint32_t index=mm->len_count[len-1]; index<mm->len_count[len]; ++index)
        {
            uint32_t seq = mm->index2seq(index);
            tree[index] = {mm->seq2index(seq, len), 0, len, seq};

            if (index>1349515 && index<1349520)
            {   
//                uint32_t seq = mm->index2seq(index);
//                tree[index] = {mm->seq2index(seq, len), 0, len, seq};
                std::cerr << "idx: " << index << "\n";
                print_node(&tree[index]);
                std::cerr << "mlen: " << max_len << "\n";

//               if (index==24517) exit(1);
            }
        }

    //    if (len==7) exit(1);    
    }
    
    
};

/**
 * Destructor.
 */
MotifTree::~MotifTree()
{
    if (tree) delete tree; 
    if (mm) delete mm;
};

/**
 * Construct suffix tree based on sequence.
 */
void MotifTree::set_sequence(char* sequence)
{
    //translate sequence to binary form
    uint32_t len = strlen(sequence);



};

/**
 * Construct suffix tree based on sequence up to max_motif_len.
 */
void MotifTree::set_sequence(char* sequence, int32_t max_motif_len)
{

};

/**
 * Gets candidate motifs up to max_motif_len.
 */
void MotifTree::get_candidate_motifs(std::vector<CandidateMotif>& candidate_motifs)
{

};

/**
 * Gets index of child.
 */
uint32_t MotifTree::get_first_child(uint32_t index)
{
    return 0;
}

/**
 * Adds a suffix of sequence from start to end.
 */
void MotifTree::add_suffix(char* sequence, int32_t start, int32_t end)
{
    
    

};

/**
 * Converts base to index.
 */
int32_t MotifTree::base2index(char base)
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
            return T;
    }
};

/**
 * Print node.
 */
void MotifTree::print_node(node* n)
{
    std::cerr << "seq: " << mm->seq2str(n->seq, n->len) << "\n";
    std::cerr << "cnt: " << n->count << "\n";
    std::cerr << "len: " << n->len << "\n";
    std::cerr << "can: " << mm->seq2str(mm->canonical(n->seq, n->len), n->len) << "\n";
    std::cerr << "\n";    
}

#undef A
#undef C
#undef G
#undef T
#undef N
