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


Node::Node()
{
    parent = NULL;
    left = NULL;
    right = NULL;
    tag = {0,0,0};
};
    
Node::Node(int32_t type)
{
    parent = NULL;
    left = NULL;
    right = NULL;
    this->type = type;
    
};

Filter::Filter()
{
    this->tree = NULL;
};

Filter::Filter(std::string exp) 
{
    this->tree = NULL;
};
    
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
 * Constructs the expression tree.
 */
void Filter::parse(const char* exp)
{
    if (tree!=NULL)
    {
        delete tree;
        tree = NULL;
    }    
    
    if (tree==NULL)
    {
        tree = new Node();
        parse(exp, strlen(exp), tree);
    }
}

/**
 * Constructs the expression tree.
 */
void Filter::parse(const char* exp, int32_t len, Node * node)
{
    std::cerr << "parsing \"";
    for (int32_t i=0; i<len; ++i)
        std::cerr << exp[i] ;
    std::cerr << "\" " << len << "\n";
    
    
    //trim
    while (*exp==' ') ++exp;
    while (exp[len-1]==' ') --len;

    std::cerr << "\tafter trimming \n";
    std::cerr << "\tparsing \"";
    for (int32_t i=0; i<len; ++i)
        std::cerr << exp[i] ;
    std::cerr << "\" " << len << "\n";

    //parse literals
    if (strcmp(exp, "PASS")==0)
    {
        node->type = VT_FILTER_OP;
        kputsn(exp, len, &node->tag);
        std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strcmp(exp, "VTYPE")==0)
    {
        node->type = VT_VARIANT_TYPE_OP;
        std::cerr << "is variant_op\n";
        return;
    }

    const char* p = exp;
    const char* q = exp;

    int32_t type = -1;

    //handle brackets
    if (*p=='(')
        while (*q!=')') ++q;
    
    
    std::cerr << "\tparsing q \"";
    for (int32_t i=0; i<len-(p-exp); ++i)
        std::cerr << q[i] ;
    std::cerr << "\" " << len << "\n";

            
    while (q-p<len)
    {
        if (*q=='&' && (q-exp<len) && *(q+1)=='&')
        {
            type = VT_OP_AND;
            std::cerr << "\tis && operator\n";
            break;
        }
        else if (*q=='&')
        {
            type = VT_OP_BIT_AND;
            std::cerr << "\tis & operator\n";
            break;
        }
        else if (*q=='|' && (q-exp<len) && *(q+1)=='|')
        {
            type = VT_OP_OR;
            std::cerr << "\tis || operator\n";
            break;
        }
        else if (*q=='|')
        {
            type = VT_OP_BIT_OR;
            std::cerr << "\tis | operator\n";
            break;
        }
        else if (*q=='=' && (q-exp<len) && *(q+1)=='=')
        {
            type = VT_OP_OR;
            std::cerr << "\tis == operator\n";
            break;
        }
        
        ++q;
    }
    
    if (q-p==len)
    {
        for (int32_t i=0; i<len; ++i)
            std::cerr << exp[i] ;
        std::cerr << "\" " << len;
        
        std::cerr << " not recognised\n";
        exit(1);
    }    
    
    //std::cerr << "p,q address: " << &p << " " <<  &q << " " << (q-p) << "\n";
    
    node->type = type;
    node->left = new Node();
    parse(p, q-p, node->left);
        
    
    if (false)
    {
        node->right = new Node();
        parse(q+2, len-(q+2-exp), node->right);
    }    
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
 * Evaluates the actions for this node.
 */
void Node::evaluate(bcf_hdr_t *h, bcf1_t *v, Variant *variant)
{
    //filter operation
    if (type==VT_FILTER_OP)
    {
        if (bcf_has_filter(h, v, tag.s)!=1)
        {
            value = false;
        }
        else
        {
            value = true;
        }
    }
    else if (type==VT_VARIANT_TYPE_OP)
    {
        if (variant->type==i)
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
 * Constructs the expression tree.
 */
//void Filter::parse(const char* filter_expression)
//{
//    
//    
//    
//    //if ()
//    
//    
////        #define LE 0 <=
////        #define LT 1 <
////        #define EQ 2 ==
////        #define GE 3 >=
////        #define GT 4 >
////
////    kstring_t tag;
////    tag.s=0; tag.l=tag.m=0;
////
////    kstring_t value;
////    value.s=0; value.l=value.m=0;
////
////    comparison = 0;
////    char c = filter.at(0);
////    int32_t mode = 0;
////    for (uint32_t i=0; i<filter.size(); ++i)
////    {
////        c = filter.at(i);
////        if (c=='<' || c=='=' || c=='>')
////        {
////            if (c=='=')
////            {
////                ++comparison;
////            }
////            else if (c=='<')
////            {
////                comparison = LT;
////            }
////            else if (c=='>')
////            {
////                comparison = GE;
////            }
////
////            mode = 1;
////        }
////        else if (mode==0)
////        {
////            kputc(c, &tag);
////        }
////        else
////        {
////            kputc(c, &value);
////        }
////    }
////
////    this->tag = std::string(tag.s);
////    this->value = atof(value.s);
////
////    if (tag.s) free(tag.s);
////    if (value.s) free(value.s);
//
////        std::cerr << filter << "\n";
////        std::cerr << this->tag << ":";
////        std::cerr << this->comparison << ":";
////        std::cerr << this->value << "\n";
//}

