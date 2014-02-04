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

    //trim brackets
    if (*exp=='(' && exp[len-1]==')')
    {
        ++exp;
        len -=2;
    }

    std::cerr << "\tafter trimming : \"";
    for (int32_t i=0; i<len; ++i)
        std::cerr << exp[i] ;
    std::cerr << "\" " << len << "\n";

    //parse literals
    if (strncmp(exp, "PASS", len)==0)
    {
        node->type = VT_FILTER_OP;
        kputsn(exp, len, &node->tag);
        std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strncmp(exp, "VTYPE", len)==0)
    {
        node->type = VT_VARIANT_TYPE_OP;
        std::cerr << "\tis variant_op\n";
        return;
    }
    else if (strncmp(exp, "N_ALLELE", len)==0)
    {
        node->type = VT_N_ALLELE_OP;
        std::cerr << "\tis nallele_op\n";
        return;
    }
    else if (strncmp(exp, "INDEL", len)==0)
    {
        node->type = VT_VARIANT_TYPE_OP;
        node->i = VT_INDEL;
        std::cerr << "\tis INDEL\n";
        return;
    }
    else if (strncmp(exp, "SNP", len)==0)
    {
        node->type = VT_VARIANT_TYPE_OP;
        node->i = VT_SNP;
        std::cerr << "\tis SNP\n";
        return;
    }
    else if (strncmp(exp, "MNP", len)==0)
    {
        node->type = VT_VARIANT_TYPE_OP;
        node->i = VT_MNP;
        std::cerr << "\tis MNP\n";
        return;
    }
    else if (strncmp(exp, "CLUMPED", len)==0)
    {
        node->type = VT_VARIANT_TYPE_OP;
        node->i = VT_SNP;
        std::cerr << "\tis CLUMPED\n";
        return;
    }
    else
    {
        const char* start = exp;
        char *end = NULL;
        int32_t i = std::strtoul(exp, &end, 10);
        if (end!=start)
        {
            node->i = i;
            std::cerr << "\tis int\n";
            return;
        }
        
        start = exp;
        end = NULL;
        int32_t f = std::strtod(exp, &end);
        if (end!=start)
        {
            node->f = f;
            std::cerr << "\tis float\n";
            return;
        }
        
        //is this a number
    }
        
    const char* p = exp;
    const char* q = exp;

    int32_t type = -1;

    //handle brackets
    if (*p=='(')
        while (*q!=')') ++q;
    
    int32_t oplen = 0;
    while (q-p<len)
    {
        if (*q=='&' && (q-exp<len) && *(q+1)=='&')
        {
            type = VT_OP_AND;
            oplen = 2;
            std::cerr << "\tis && operator\n";
            break;
        }
        else if (*q=='&')
        {
            type = VT_OP_BIT_AND;
            oplen = 1;
            std::cerr << "\tis & operator\n";
            break;
        }
        else if (*q=='|' && (q-exp<len) && *(q+1)=='|')
        {
            type = VT_OP_OR;
            oplen = 2;
            std::cerr << "\tis || operator\n";
            break;
        }
        else if (*q=='|')
        {
            type = VT_OP_BIT_OR;
            oplen = 1;
            std::cerr << "\tis | operator\n";
            break;
        }
        else if (*q=='=' && (q-exp<len) && *(q+1)=='=')
        {
            type = VT_OP_EQ;
            oplen = 2;
            std::cerr << "\tis == operator\n";
            break;
        }
        else
        {
            
        }
        
        ++q;
    }

    //std::cerr << "p,q address: " << &p << " " <<  &q << " " << (q-p) << "\n";
    
    node->type = type;
    node->left = new Node();
    parse(p, q-p, node->left);
    
    if (q-p!=len)
    {
        q += oplen;
        std::cerr << "\tparsing q \"";
        for (int32_t i=0; i<len-(q-p); ++i)
            std::cerr << q[i] ;
        std::cerr << "\" " << (len-(q-p)) << "\n";
            
            
        node->right = new Node();
        parse(q, len-(q-exp), node->right);
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
    if (type==VT_OP_AND)
    {
        value = (left->value==right->value);
    }
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
        i = variant->type;
    }
    else if (type==VT_N_ALLELE_OP)
    {
        i = bcf_get_n_allele(v);
    }
    else if (type==VT_OP_EQ)
    {
        if (left->type==VT_N_ALLELE_OP || left->type==VT_VARIANT_TYPE_OP)
        {
            value = (left->i==right->i);
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

