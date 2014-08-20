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

#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "variant_manip.h"

#ifndef SV_Tree_H
#define SV_Tree_H

#define VT_SV_DEL    1
#define VT_SV_INS    2
#define VT_SV_DUP    3
#define VT_SV_INV    4
#define VT_SV_CNV    5
#define VT_SV_TRA    6
#define VT_SV_TANDEM 7
#define VT_SV_ME     8
#define VT_SV_MT     9

/**
 * Class for filtering VCF records.
 */
class SVNode
{
    public:

    SVNode* parent;
    std::vector<SVNode*> children;
    int32_t depth, count;
    kstring_t desc;

    /**
     * Constructor.
     */
    SVNode();

    /**
     * Constructor.
     */
    SVNode(const char* desc);

    /**
     * Constructor.
     */
    SVNode(const char* desc, int32_t depth);

    /**
     * Destructor.
     */
    ~SVNode();

    /**
     * Set depth.
     */
    void set_depth(int32_t depth);

    /**
     * Increment count.
     */
    void increment_count();
    
    /**
     * Clear values.
     */
    void clear();
    
    /**
     * Print values.
     */
    void print();
    
    /**
     *  For translating reserved keywords.
     */
    const char* tags2desc();    
};

KHASH_MAP_INIT_STR(xdict, SVNode*);

/**
 * SV Tree for counting SVs.
 */
class SVTree
{
    public:

    SVNode* root;
    int32_t max_depth;    
    std::vector<SVNode*> df_order;
//    bcf_hdr_t *h;
//    bcf1_t *v;
    khash_t(xdict) *m;

    /**
     * Constructor.
     */
    SVTree();

    /**
     * Destructor.
     */
    ~SVTree();

    /**
     * Adds a new tag, returns true if successful.
     */
    bool add(const char* desc);
    
    /**
     * Observes and update the count of a new tag.
     */
    void count(char* desc);

    /**
     * Enumerates the children in a depth first order.
     */
    SVNode* enumerate();

    /**
     * Iterator, returns first node by depth first search.
     */
    SVNode* begin();

    /**
     * Iterator, returns node by depth first search.
     */
    SVNode* next();

    /**
     * Print this tree.
     */
    void print();

    private:
};

#endif