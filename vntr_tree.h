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

#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "variant_manip.h"
#include "vntr_annotator.h"

#ifndef VNTR_TREE_H
#define VNTR_TREE_H

#define SEQUENCE 1 
#define MOTIF    2
#define BASIS    3

/**
 * Class for filtering VCF records.
 */
class VNTRNode
{
    public:
    std::string motif;
    std::string basis;
    int32_t exact_count, fuzzy_count;


    /**
     * Constructor.
     */
    VNTRNode();

    /**
     * Constructor.
     */
    VNTRNode(std::string motif, std::string basis, int32_t exact_count, int32_t fuzzy_count);

    /**
     * Destructor.
     */
    ~VNTRNode();

    /**
     * Clear values.
     */
    void clear();

    /**
     * Print values.
     */
    void print();
};

/**
 * VNTR Tree for counting VNTRs.
 */
class VNTRTree
{
    public:
    
    std::vector<std::list<VNTRNode*> > vntrs[4];
    std::map<std::string, VNTRNode*> motif_map;
         
    /**
     * Constructor.
     */
    VNTRTree();

    /**
     * Destructor.
     */
    ~VNTRTree();

    /**
     * Assumes that the variant is appropriately updated on its VNTR chracteristics.
     * Updates tree on motif, exact/fuzziness of VNTR.
     */
    void count(Variant& variant);

    /**
     * Print this tree.
     */
    void print(int32_t level);

    private:
};

#endif