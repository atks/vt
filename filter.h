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

#ifndef FILTER_H
#define FILTER_H

#include "htslib/vcf.h"
#include "variant_manip.h"

//ordered by precedence level
//binary ops
#define VT_OP_NOT      0
#define VT_OP_AND      1
#define VT_OP_OR       2
#define VT_OP_BIT_AND  3
#define VT_OP_BIT_OR   4

#define VT_OP_EQ       5
#define VT_OP_NE       6
#define VT_OP_LT       7
#define VT_OP_LE       8
#define VT_OP_GT       9
#define VT_OP_GE       10

//binary math ops
#define ADD  11
#define SUB  12
#define MUL  13
#define DIV  14

//unary ops (data getters for vcf)
#define VT_VARIANT_TYPE_OP  65
#define VT_N_ALLELE_OP      66
#define VT_VARIANT_DLEN_OP  67
#define VT_INFO_INT_OP      68
#define VT_INFO_FLT_OP      69
#define VT_INFO_STR_OP      70
#define VT_FILTER_OP        257

#define VT_INT              64
#define VT_FLT              128
#define VT_STR              256




/**
 * Class for filtering VCF records.
 * Filters based on several fields:
 * QUAL
 * FILTER
 * INFO
 * VARIANT (inferred)
 *
 * QUAL>40
 * FILTER==PASS
 * VARIANT==SNP
 * AF>0.5
 * VARIANT==SNP && AF>0.5

 [-F] filter expression

based on QUAL, FILTER, INFO and Variant type


-f QUAL>2&&PASS&&AF>0.05
-f PASS && AF*>0.05
-f PASS && (AF*>0.05 || AC/AN>0.05)
-f PASS && SNP

reserved key words
PASS
QUAL
AF
AC
AN

-f, -g, -h
In the case that a field is not found, it evaluates to false
AF* - intelligent parsing - if AF is not present, estimate from AC/AN
 

Actual working examples

(N_ALLELE==2)&&(VTYPE==INDEL)&&PASS

 *
 */
class Node
{
    public:

    Node* parent;
    Node* left;
    Node* right;

    int32_t type;
    bool value;
    
    kstring_t tag;
    union
    {
        bool b;
        int32_t i;
        float f;
    };
    
    
    Node();
    
    Node(int32_t type);
    
    /**
     * Evaluates the actions for this node.
     */
    void evaluate(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug=false);
};

class Filter
{
    public:

    //encodes the filter expression
    Node* tree;

    //useful pointers
    bcf_hdr_t *h;
    bcf1_t *v;
    Variant *variant;

    Filter();

    Filter(std::string exp);

    /**
     * Applies filter to vcf record.
     */
    bool apply(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug=false);

    /**
     * Recursive call for apply.
     */
    void apply(Node* node, bool debug=false);

    /**
     * Constructs the expression tree.
     */
    void parse(const char* exp, bool debug=false);
    void parse(const char* exp, int32_t len, Node * node, bool debug=false);
    
    /**
     * Parse literals.
     */
    void parse_literal(const char* exp, int32_t len, Node * node, bool debug=false);
        
    bool is_literal(const char* exp, int32_t len);
    
    void update_node(int32_t OP, Node* node)
    {
    }
    
    /**
     * Recursive call for parse.
     */
    
    
};

#endif