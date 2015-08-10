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
MotifTree::MotifTree(uint32_t max_len, bool debug)
{
    this->max_len = max_len;
    mm = new MotifMap(max_len);

    this->debug = debug;

    lc.resize(max_len+1,0);

    tree = (node *) malloc(sizeof(node)*mm->max_index+1);

    //perform mapping.
    for (uint32_t len=1; len<=max_len; ++len)
    {
        for (uint32_t index=mm->len_count[len-1]; index<mm->len_count[len]; ++index)
        {
            uint32_t seq = mm->index2seq(index);
            tree[index] = {index, mm->seq2index(mm->canonical(seq, len), len), 0, len, seq};

            uint32_t c = mm->canonical(seq, len);
            //bool a = mm->is_aperiodic(c, len);
            //std::cerr << " " << mm->seq2str(c, len) << " " << mm->is_aperiodic(c, len) << "\n";
        }
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
 * Gets subsequence from a string at the spos0 position of length cmax_len.
 */
uint32_t MotifTree::get_sub_seq(char* seq, uint32_t len, uint32_t spos0, uint32_t& s)
{
    if (spos0==0)
    {
        s = 0;
        for (uint32_t i=0; i<cmax_len; ++i)
        {
            uint32_t j = base2index(seq[i]);
            s = set_seqi(s, i, j);
        }

        return cmax_len;
    }
    else
    {
        uint32_t slen = 0;

        //shorter than max_len
        if (len-spos0<cmax_len)
        {
            slen = len-spos0;
            s = shift(s, slen+1);
            s = set_seqi(s, len-spos0, 0);
        }
        else
        {
            slen = cmax_len;
            s = shift(s, cmax_len);
            //std::cerr << seq[spos0+max_len-1] << " " << spos0+max_len-1 <<" " << base2index(seq[spos0+max_len-1]) << "\n";
            s = set_seqi(s, cmax_len-1, base2index(seq[spos0+cmax_len-1]));
        }

        return slen;
    }
}

/**
 * Construct suffix tree based on sequence of length len.
 */
void MotifTree::set_sequence(char* seq, uint32_t len)
{
    //computes the relevant maximum possible length of motif to check
    cmax_len = (len >> 1) < max_len ? (len >> 1) : max_len;

    if (debug)
    {
        std::cerr << "len : " << len << "\n";
        std::cerr << "cmax_len : " << cmax_len << "\n";
    }

    uint32_t s = 0;
    for (uint32_t i=0; i<len; ++i)
    {
        //update as uint32_t representation
        uint32_t l = get_sub_seq(seq, len, i, s);

        //mm->print_seq(s,l);

        //insert prefix
        insert_prefix(s, l);
    }
};

/**
 * Inserts prefix s into tree.
 */
void MotifTree::insert_prefix(uint32_t s, uint32_t len)
{
    for (uint32_t i=1; i<=len; ++i)
    {
        uint32_t index = mm->seq2index(s, i);
        ++tree[index].count;
    }
}

/**
 * Consolidate motif counts.
 */
void MotifTree::consolidate_motif_counts()
{
    cm.clear();
    std::fill(lc.begin(), lc.end(), 0);
    while (!pcm.empty()) pcm.pop();
    consolidate_motif_counts(&tree[0]);
    consolidate_motif_counts(&tree[1]);
    consolidate_motif_counts(&tree[2]);
    consolidate_motif_counts(&tree[3]);
}

/**
 * Consolidate motif counts.
 */
void MotifTree::consolidate_motif_counts(node* n)
{
    if (n->count)
    {
        if (cm.find(n->cindex)==cm.end())
        {
            cm[n->cindex] = n->count;
        }
        else
        {
            cm[n->cindex] += n->count;
        }

        lc[n->len] += n->count;

        n->count = 0;

        if (n->len != cmax_len)
        {
            uint32_t index = get_first_child(n->seq, n->len);

            consolidate_motif_counts(&tree[index]);
            consolidate_motif_counts(&tree[index+1]);
            consolidate_motif_counts(&tree[index+2]);
            consolidate_motif_counts(&tree[index+3]);
        }
    }
}

/**
 * Shifts a string.
 */
std::string MotifTree::shift_str(std::string& seq, uint32_t i)
{
    std::string sseq = seq;
    if (i)
    {
        sseq = seq.substr(i) + seq.substr(0,i);
    }    
    
    return sseq;
}

/**
 * Checks if two copies of a motif exists in a seq.
 */
bool MotifTree::exist_two_copies(std::string& seq, std::string& motif)
{
    //all phases of motif
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string ru = shift_str(motif, i);
        const char* s = seq.c_str();

        if ((s = strstr(s, ru.c_str())))
        {
            s += ru.size();
            
            if ((s = strstr(s, ru.c_str())))
            {
                return true;
            }
        }
    }
    
    return false;
}

/**
 * Compute fit of expected content of nucleotides.
 */
float MotifTree::compute_fit(uint32_t index, scontent* sc)
{
    uint32_t e[] = {0,0,0,0};
    node* n = &tree[index];
    
//    std::cerr << "compute fit: " << mm->seq2str(n->seq, n->len) << "\n";
    
    //compute expected fit
    for (uint32_t i=0; i<n->len; ++i)
    {
        ++e[get_seqi(n->seq, i)];
    }
    
//    std::cerr << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << "\n";
    
    float fit = 0;
    float t = 0;
    
//    std::cerr << "n " << sc->n << " (" << sc->base[A] << "," << sc->base[C] << "," << sc->base[G] << "," << sc->base[T] << ")\n";
    
    for (uint32_t i=0; i<4; ++i)
    {
        t = (float)sc->base[i]/sc->n - (float)e[i]/n->len;
        fit += t*t;
    }
    
//    std::cerr << "compute fit " << fit << "\n";
    
    return fit;
}

/**
 * Gets candidate motifs up to max_motif_len.
 */
void MotifTree::detect_candidate_motifs(std::string& seq)
{
    detect_candidate_motifs(const_cast<char*>(seq.c_str()), seq.size());
}

/**
 * Detects candidate motifs from seq of length len.
 */
void MotifTree::detect_candidate_motifs(char* seq, uint32_t len)
{
    if (debug)
    {
        std::cerr << "detecting motifs for an str\n";
        std::cerr << "seq: " << seq << "\n";
    }

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

//    print_tree();

    //gather distribution
    scontent sc;
    sc.base[A] = tree[A].count;
    sc.base[C] = tree[C].count;
    sc.base[G] = tree[G].count;
    sc.base[T] = tree[T].count;
    sc.n = sc.base[A] + sc.base[C] + sc.base[G] + sc.base[T];
    
    //count occurrences
    consolidate_motif_counts();

    if (debug)
        std::cerr << "candidate motifs: " << cm.size() << "\n";

    float sthreshold = std::min(0.80, (s.size()-1.0)/s.size());

    for (std::map<uint32_t, uint32_t>::iterator i = cm.begin(); i!=cm.end(); ++i)
    {
        if (mm->is_aperiodic(mm->index2seq(i->first), tree[i->first].len))
        {
            //p - purity proportion
            float p = (float)i->second/(lc[tree[i->first].len]);
            //f - fit based on content
            float f = compute_fit(i->first, &sc);
            //p -= f;
             
            //if ((tree[i->first].len==1 && p+f>sthreshold) || (tree[i->first].len>1 && p>0.2))
            if (len<10 || (tree[i->first].len==1 && p>0.6) || (tree[i->first].len>1))
            {
                std::string motif = mm->seq2str(mm->index2seq(i->first), tree[i->first].len);
                if (exist_two_copies(s, motif))
                {
                    if (debug) std::cerr << motif << " : " << p << " " << tree[i->first].len << " " << f << "\n";
                    pcm.push(CandidateMotif(motif, p, tree[i->first].len, f));
                }
                else
                {
                    if (debug) std::cerr << motif << " : " << p << " " << tree[i->first].len << " " << f << " (< 2 copies)\n";
                }
            }
        }
    }

    //if no pickups
    if (pcm.size()==0)
    {
        std::map<uint32_t, uint32_t>::iterator i = cm.begin();
        
        float p = (float)i->second/(lc[tree[i->first].len]);
        float f = compute_fit(i->first, &sc);
        p -= f;
        if (debug) std::cerr << mm->seq2str(mm->index2seq(i->first), tree[i->first].len) << " : " << p << " " << tree[i->first].len << " " << f << "\n";
        pcm.push(CandidateMotif(mm->seq2str(mm->index2seq(i->first), tree[i->first].len), p, tree[i->first].len, f));
            
    }    

//    print_tree();
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
