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

#include "filter.h"

#define BOOLEAN_OP  0
#define MATH_OP     1
#define VARIANT_TYPE_OP  2 
#define INFO_OP     3
#define FILTER_OP   4

#define NOT  0
#define AND  1
#define OR   2

#define ADD  3
#define SUB  4
#define MUL  5
#define DIV  6

/**
 * Evaluates the actions for this node.
 */
void Node::evaluate(bcf_hdr_t *h, bcf1_t *v, Variant *variant)
{
    //filter operation
    if (type==FILTER_OP)
    {
        if (bcf_has_filter(h, v, tag)!=1)
        {
            value = true;
        }
        else
        {
            value = false;
        }
    }
    else if (type==VARIANT_TYPE_OP)
    {
        if (variant->type==type)
        {
            value = true;
        }
        else
        {
            value = false;
        }
    }
    
//    float *f = 0;
//    int32_t n = 0;
//    int32_t *ac = 0;
//    int32_t *an = 0;
//    if (bcf_get_info_float(h, v, tag.c_str(), &f, &n))
//    {
//        switch (comparison)
//        {
//            case LT : return *f<value;
//            case LE : return *f<=value;
//            case EQ : return *f==value;
//            case GT : return *f>value;
//            case GE : return *f>=value;
//            default : return true;
//        }
//    }
//    else if (bcf_get_info_int(h, v, "AC", &ac, &n) && bcf_get_info_int(h, v, "AN", &an, &n))
//    {
//        *f = (float)(*ac)/(float)(*an);
//
//        switch (comparison)
//        {
//            case LT : return *f<value;
//            case LE : return *f<=value;
//            case EQ : return *f==value;
//            case GT : return *f>value;
//            case GE : return *f>=value;
//            default : return true;
//        }
//    }
    
}

/**
 * Applies filter to vcf record.
 */
bool Filter::apply(bcf_hdr_t *h, bcf1_t *v, Variant *variant) //recursive
{
    if (tree==NULL)
    {
        return true;
    }
    
    this->h = h;
    this->v = v;
    this->variant = variant;    
    
    apply(tree);
    
    return tree->value;
}

/**
 * Recursive call for apply.
 */
void Filter::apply(Node* node)
{
    //evaluate downstream
    if (node->left!=NULL)
    {
        apply(node->left);
    }    
    
    if (node->right!=NULL)
    {
        apply(node->right);
    }
    
    node->evaluate(h, v, variant);
}


/**
 * Constructs the expression tree.
 */
void Filter::parse(const char* filter_expression)
{
    
    //if ()
    
    
//        #define LE 0 <=
//        #define LT 1 <
//        #define EQ 2 ==
//        #define GE 3 >=
//        #define GT 4 >
//
//    kstring_t tag;
//    tag.s=0; tag.l=tag.m=0;
//
//    kstring_t value;
//    value.s=0; value.l=value.m=0;
//
//    comparison = 0;
//    char c = filter.at(0);
//    int32_t mode = 0;
//    for (uint32_t i=0; i<filter.size(); ++i)
//    {
//        c = filter.at(i);
//        if (c=='<' || c=='=' || c=='>')
//        {
//            if (c=='=')
//            {
//                ++comparison;
//            }
//            else if (c=='<')
//            {
//                comparison = LT;
//            }
//            else if (c=='>')
//            {
//                comparison = GE;
//            }
//
//            mode = 1;
//        }
//        else if (mode==0)
//        {
//            kputc(c, &tag);
//        }
//        else
//        {
//            kputc(c, &value);
//        }
//    }
//
//    this->tag = std::string(tag.s);
//    this->value = atof(value.s);
//
//    if (tag.s) free(tag.s);
//    if (value.s) free(value.s);

//        std::cerr << filter << "\n";
//        std::cerr << this->tag << ":";
//        std::cerr << this->comparison << ":";
//        std::cerr << this->value << "\n";
}

/**
 * Recursive call for parse.
 */
void parse(const char* filter, int32_t len)
{
}