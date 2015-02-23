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

/**
 * Class for filtering VCF records.
 */
class SVNode
{
    public:

    SVNode* parent;
    std::vector<SVNode*> children;
    int32_t depth, count, mcount;
    kstring_t type;

    /**
     * Constructor.
     */
    SVNode();

    /**
     * Constructor.
     */
    SVNode(const char* type);

    /**
     * Constructor.
     */
    SVNode(const char* type, int32_t depth);

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
     * Increment mcount.
     */
    void increment_mcount();

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
    std::string sv_type2string(char* sv_type);
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
    khash_t(xdict) *m;
    int32_t mixed_sv_count;

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
    void count(Variant& variant);

    /**
     * Enumerates the children in a depth first order.
     */
    std::vector<SVNode*> enumerate_dfs();

    /**
     * Helper function for enumerating the children in a depth first order.
     */
    void enumerate_dfs(std::vector<SVNode*>& s, SVNode* node);

    /**
     * Print this tree.
     */
    void print();

    private:
};

#endif