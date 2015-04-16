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
#define N 4

/**
 * Constructor.
 */
MotifTree::MotifTree(uint32_t max_len)
{   
    this->max_len = max_len;
    mm = new MotifMap(max_len); 

    tree = (node *) malloc(sizeof(node)*mm->max_index+1);
        
    //perform mapping.
    for (uint32_t len=1; len<=max_len; ++len)
    {
        for (uint32_t index=mm->len_count[len-1]; index<mm->len_count[len]; ++index)
        {
            uint32_t seq = mm->index2seq(index);
            tree[index] = {mm->seq2index(seq, len), 0, len, seq};

//            if (index>1349515)
//            {   
//                uint32_t seq = mm->index2seq(index);
//                tree[index] = {mm->seq2index(seq, len), 0, len, seq};
//                std::cerr << "idx: " << index << "\n";
//                print_node(&tree[index]);
//                std::cerr << "mlen: " << max_len << "\n";

//               if (index==24517) exit(1);
//            }
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
 * Gets subsequence from a C string.
 */
uint32_t MotifTree::get_sub_seq(char* seq, uint32_t len, uint32_t spos0, uint32_t& s)
{
    if (spos0==0)
    {
        s = 0;
        uint32_t slen = std::min(max_len, len);
        for (uint32_t i=0; i<slen; ++i)
        {
            uint32_t j = base2index(seq[i]);
            s = set_seqi(s, i, j);
        }
        
        return slen;
    }   
    else
    {
        uint32_t slen = 0;
         
        //shorter than max_len 
        if (len-spos0<max_len)
        {
            slen = len-spos0;
            s = shift(s, slen);
            s = set_seqi(s, len-spos0, 0);
        }
        else
        {
            slen = max_len;
            s = shift(s, slen);
            //std::cerr << seq[spos0+max_len-1] << " " << spos0+max_len-1 <<" " << base2index(seq[spos0+max_len-1]) << "\n";
            s = set_seqi(s, max_len-1, base2index(seq[spos0+max_len-1]));
        }   
        
        return slen;
    } 
    
}

/**
 * Inserts prefix s into tree.
 */
void MotifTree::insert_prefix(uint32_t s, uint32_t len)
{
    for (uint32_t i=1; i<=len; ++i)
    {
        uint32_t index = mm->seq2index(s, i);
        //std::cerr << index << " " << mm->seq2str(s,i) << "\n";
        ++tree[index].count;
    }
}

/**
 * Construct suffix tree based on sequence up to max_motif_len.
 */
void MotifTree::set_sequence(char* seq, int32_t len)
{
    uint32_t s = 0;
    for (uint32_t i=0; i<len; ++i)
    {
        //update as uint32_t representation
        uint32_t l = get_sub_seq(seq, len, i, s);
        //std::cerr << mm->seq2str(s,l) << "\n";
        
        //insert prefix
        insert_prefix(s, l);
        //insert(s, l);
        

    }
    
    
    print_tree();
    
    //summarise results
};

/**
 * Gets candidate motifs up to max_motif_len.
 */
void MotifTree::detect_candidate_motifs(std::vector<CandidateMotif>& candidate_motifs, char* seq, uint32_t max_motif_len)
{
    std::cerr << "============================================\n";
    std::cerr << "detecting motifs for an str\n";
    std::cerr << "seq: " << seq << "\n";
    
    std::string s(seq);
    if (strchr(seq, 'N'))
    {
        uint32_t len = s.size();
        s.clear();
        for (uint32_t i=0; i<len; ++i)
        {
            if (seq[i]!='N')
            {
                s.append(1,seq[i]); 
            }
        }
    }

    
    //construct tree
    set_sequence(const_cast<char*>(s.c_str()), s.size());
    
    //count occurrences
    
    
    
    //populate candidate_motifs   
    
    
    
    std::cerr << "============================================\n";
    
};

/**
 * Gets index of first child.
 */
uint32_t MotifTree::get_first_child(uint32_t seq, uint32_t len)
{
    //not necessary sequence should be zeroed out after valid length.
    //uint32_t cseq = set_seq(s, len, 0);
    uint32_t index = mm->seq2index(seq, len+1);
    return index;
}

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
            return N;
    }
};

/**
 * Print node.
 */
void MotifTree::print_node(node* n)
{
    std::cerr << std::string(n->len, '\t') << "seq: " << mm->seq2str(n->seq, n->len) << "\n";
    std::cerr << std::string(n->len, '\t') << "cnt: " << n->count << "\n";
    std::cerr << std::string(n->len, '\t') << "len: " << n->len << "\n";
    std::cerr << std::string(n->len, '\t') << "can: " << mm->seq2str(mm->canonical(n->seq, n->len), n->len) << "\n";
    std::cerr << "\n";    
}

/**
 * Print tree.
 */
void MotifTree::print_tree()
{
    print_tree(&tree[0]);
    print_tree(&tree[1]);
    print_tree(&tree[2]);
    print_tree(&tree[3]);
}

/**
 * Print tree.
 */
void MotifTree::print_tree(node* n)
{
    if (n->count)
    {
        print_node(n);
    
        if (n->len != max_len)
        {
            uint32_t index = get_first_child(n->seq, n->len);
        
            print_tree(&tree[index]);
            print_tree(&tree[index+1]);
            print_tree(&tree[index+2]);
            print_tree(&tree[index+3]);
        }
    }
}

#undef A
#undef C
#undef G
#undef T
#undef N
