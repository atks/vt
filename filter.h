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
#include "htslib/kstring.h"
#include "variant_manip.h"

//TYPES
#define VT_LOGIC_OP  2048
#define VT_MATH_CMP  4096
#define VT_MATH_OP   8192
#define VT_BCF_OP    16384

#define VT_BOOL     64
#define VT_INT      128
#define VT_FLT      256
#define VT_STR      512
#define VT_FLG      1024

//common unary and binary ops (ordered by precedence level)
#define VT_NOT      (0|VT_LOGIC_OP)
#define VT_AND      (1|VT_LOGIC_OP)
#define VT_OR       (2|VT_LOGIC_OP)

#define VT_EQ       (3|VT_MATH_CMP)
#define VT_NE       (4|VT_MATH_CMP)
#define VT_LT       (5|VT_MATH_CMP)
#define VT_LE       (6|VT_MATH_CMP)
#define VT_GT       (7|VT_MATH_CMP)
#define VT_GE       (8|VT_MATH_CMP)

#define VT_ADD      (9|VT_MATH_OP)
#define VT_SUB      (10|VT_MATH_OP)
#define VT_MUL      (11|VT_MATH_OP)
#define VT_DIV      (12|VT_MATH_OP)
#define VT_BIT_AND  (13|VT_INT|VT_MATH_OP)
#define VT_BIT_OR   (14|VT_INT|VT_MATH_OP)

//unary ops (data getters for vcf)
#define VT_FILTER        (33|VT_BCF_OP)
#define VT_INFO          (34|VT_BCF_OP)
#define VT_N_ALLELE      (35|VT_INT|VT_BCF_OP)
#define VT_VARIANT_TYPE  (36|VT_INT|VT_BCF_OP)
#define VT_VARIANT_DLEN  (37|VT_INT|VT_BCF_OP)
#define VT_VARIANT_LEN   (38|VT_INT|VT_BCF_OP)

#define VT_UNKNOWN -1

/**
 * Class for filtering VCF records.
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
    kstring_t s;
    bool b;
    int32_t i;
    float f;

    /**
     * Constructor.
     */
    Node();

    /**
     * Constructor with type initlialization.
     */
    Node(int32_t type);

    /**
     * Evaluates the actions for this node.
     */
    void evaluate(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug=false);
};

/**
 * Filter for VCF records.
 */
class Filter
{
    public:

    //filter expression
    Node* tree;

    //useful pointers for applying the filter to a vcf record
    bcf_hdr_t *h;
    bcf1_t *v;
    Variant *variant;
    bool need_to_classify_variant;

    /**
     * Constructor.
     */
    Filter();

    /**
     * Constructor with expression initialization.
     */
    Filter(std::string exp, VariantManip *vm=NULL);

    /**
     * Parses filter expression.
     */
    void parse(const char* exp, bool debug=false);

    /**
     * Applies filter to vcf record.
     */
    bool apply(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug=false);

    /**
     * Attempts to simplify the expression tree by collapsing nodes that can be precomputed.
     */
    void simplify();

    private:

    /**
     * Recursive call for parse.
     */
    void parse(const char* exp, int32_t len, Node * node, bool debug=false);

    /**
     * Parse literals.
     */
    bool is_literal(const char* exp, int32_t len);

    /**
     * Parse literals.
     */
    void parse_literal(const char* exp, int32_t len, Node * node, bool debug=false);

    /**
     * Trim brackets from an expression.
     */
    void trim_brackets(const char* &exp, int32_t &len);

    /**
     * Moves r to the closing bracket if this expression starts with an open bracket.
     * Returns -1 if end of r else 0.
     */
    int32_t fwd_to_closing_bracket(const char* &r, int32_t &len);

    /**
     * Returns -1 if no operator found. Updates oplen to be the length of the operator observed.
     */
    int32_t peek_op(const char* &r, int32_t len, int32_t &oplen, bool debug);

    /**
     * Recursive call for apply.
     */
    void apply(Node* node, bool debug=false);
    
    /**
     * Help message on filter expressions.
     */
    void print_filter_help();
};

#endif