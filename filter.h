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

#include <algorithm>
#include <cctype>
#include "hts_utils.h"
#include "variant.h"
#include "pregex.h"

//TYPES
//this is defined from the 7th to 13th bit by setting the bit
#define VT_LOGIC_OP  2048   //0x1000
#define VT_MATH_CMP  4096   //0x2000
#define VT_MATH_OP   8192   //0x4000
#define VT_BCF_OP    16384  //0x8000

#define VT_BOOL     64      //0x0040
#define VT_INT      128     //0x0080 
#define VT_FLT      256     //0x0100
#define VT_STR      512     //0x0200
#define VT_FLG      1024    //0x0400

//common unary and binary ops (ordered by precedence level)
//this is identified over the first 6 bits which gives 64 possible operations
#define VT_NOT      (0|VT_LOGIC_OP|VT_BOOL)
#define VT_AND      (1|VT_LOGIC_OP|VT_BOOL)
#define VT_OR       (2|VT_LOGIC_OP|VT_BOOL)

#define VT_EQ       (3|VT_MATH_CMP|VT_BOOL)
#define VT_NE       (4|VT_MATH_CMP|VT_BOOL)
#define VT_LT       (5|VT_MATH_CMP|VT_BOOL)
#define VT_LE       (6|VT_MATH_CMP|VT_BOOL)
#define VT_GT       (7|VT_MATH_CMP|VT_BOOL)
#define VT_GE       (8|VT_MATH_CMP|VT_BOOL)
#define VT_MATCH    (9|VT_MATH_CMP|VT_BOOL)
#define VT_NO_MATCH (10|VT_MATH_CMP|VT_BOOL)
#define VT_ADD      (11|VT_MATH_OP|VT_FLT)
#define VT_SUB      (12|VT_MATH_OP|VT_FLT)
#define VT_MUL      (13|VT_MATH_OP|VT_FLT)
#define VT_DIV      (14|VT_MATH_OP|VT_FLT)
#define VT_BIT_AND  (15|VT_INT|VT_MATH_OP|VT_BOOL)
#define VT_BIT_OR   (16|VT_INT|VT_MATH_OP|VT_BOOL)

//unary ops (data getters for vcf)
#define VT_VARIANT_TYPE        (33|VT_INT|VT_BCF_OP|VT_BOOL)
#define VT_VARIANT_DLEN        (34|VT_INT|VT_BCF_OP)
#define VT_VARIANT_LEN         (35|VT_INT|VT_BCF_OP)
#define VT_VARIANT_CONTAINS_N  (36|VT_BCF_OP|VT_BOOL)
#define VT_N_ALLELE            (37|VT_INT|VT_BCF_OP)
#define VT_QUAL                (38|VT_BCF_OP|VT_FLT)
#define VT_FILTER              (39|VT_BCF_OP|VT_BOOL)
#define VT_N_FILTER            (40|VT_INT|VT_BCF_OP)
#define VT_INFO                (41|VT_BCF_OP)
#define VT_REF_COL             (42|VT_BCF_OP|VT_STR)
#define VT_ALT                 (43|VT_BCF_OP|VT_STR)

//problems will arise once you pass 63.
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

//VT_LOGIC_OP  2048
//VT_MATH_CMP  4096
//VT_MATH_OP   8192
//VT_BCF_OP    16384
//
//VT_BOOL     64
//VT_INT      128
//VT_FLT      256
//VT_STR      512
//#define VT_FLG      1024
    int32_t type;  // data type
//BCF_VL_FIXED 0  - implemented this
//BCF_VL_VAR   1  todo: implement?
//BCF_VL_A     2  todo: implement
//BCF_VL_G     3  todo: implement
//BCF_VL_R     4  todo: implement
    int32_t var_length;     //variable length
    int32_t number; //actual length
    kstring_t tag;  //store the INFO tag of a BCF type
    int32_t index;  //store index value of interest

    bool value_exists; // if value exists

    kstring_t s;   // string value
    bool b;        // boolean value
    int32_t i;     // integer value
    float f;       // float value

    //for cases of fix multiple values
    std::vector<std::string> svec;
    std::vector<bool> bvec;
    std::vector<int32_t> ivec;
    std::vector<float> fvec;

    PERLregex pregex;
    bool regex_set;

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

    /**
     * Converts type to string.
     */
    std::string type2string(int32_t type);
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

    /**
     * Constructor.
     */
    Filter();

    /**
     * Constructor with expression initialization.
     */
    Filter(std::string exp);

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

    /**
     * Resets filter.
     */
    void reset();

    private:

    /**
     * Recursive call for parse.
     */
    void parse(const char* exp, int32_t len, Node * node, bool debug=false);

    /**
     * Checks if exp is a literal.
     */
    bool is_literal(const char* exp, int32_t len, bool debug=false);

    /**
     * Checks if exp is a unary op.
     */
    bool is_unary_op(const char* exp, int32_t len, bool debug=false);

    /**
     * Checks is expression is bracketed.
     */
    bool is_bracketed_expression(const char* exp, int32_t len, bool debug);

    /**
     * Detect index width.
     * e.g. AC[1] => 3
     * e.g. EVIDENCE[12] => 4
     *
     * Populated index with the index queried.
     */
    int32_t get_index_width(const char *exp, int32_t n, int32_t *index);

    /**
     * Parse literals.
     */
    void parse_literal(const char* exp, int32_t len, Node * node, bool debug=false);

    /**
     * Trim brackets from an expression.
     */
    void trim_brackets(const char* &exp, int32_t &len, bool debug=false);

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