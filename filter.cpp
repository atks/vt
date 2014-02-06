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
 * Constructs the expression tree.
 */
void Filter::parse(const char* exp, bool debug)
{
    if (tree!=NULL)
    {
        delete tree;
        tree = NULL;
    }

    if (tree==NULL)
    {
        tree = new Node();
        parse(exp, strlen(exp), tree, debug);
    }
}

void trim_brackets(const char* &exp, int32_t &len)
{
    if (*exp=='(' && exp[len-1]==')')
    {
//        std::cerr << "call trim brackets\n";
//        std::cerr << "parsing \"";
//        for (int32_t i=0; i<len; ++i)
//            std::cerr << exp[i] ;
//        std::cerr << "\" " << len << "\n";

        int32_t opened_brackets = 1;
        bool nested = true;
        int32_t j=1;
        while(j<len-1)
        {
            if(exp[j]=='(')
            {
               if (opened_brackets<0)
               {
                    kstring_t s = {0,0,0};
                    kputsn(exp, len, &s);
                    fprintf(stderr, "[%s:%d %s] brackets not correct %s\n", __FILE__, __LINE__, __FUNCTION__, s.s);
                    if (s.m) free(s.s);
                    exit(1);
               }

               ++opened_brackets;
            }
            else if (exp[j]==')')
            {
               if (opened_brackets<=0)
               {
                    kstring_t s = {0,0,0};
                    kputsn(exp, len, &s);
                    fprintf(stderr, "[%s:%d %s] brackets not correct %s\n", __FILE__, __LINE__, __FUNCTION__, s.s);
                    if (s.m) free(s.s);
                    exit(1);
               }

               --opened_brackets;
            }

            if (opened_brackets==0)
            {
                nested = false;
                break;
            }
            ++j;


        }

        if (nested)
        {
            if (opened_brackets == 1)
            {
                ++exp;
                len -=2;

                trim_brackets(exp, len);
            }
            else
            {
                std::cerr << "Illegal expression: brackets not correct\n";
            }
        }
    }
}

/**
 * Moves r to the closing bracket if this expression starts with an open bracket.
 * Returns -1 if end of r else 0.
 */
int32_t Filter::fwd_to_closing_bracket(const char* &r, int32_t &len)
{
    const char* s = r;
    if (*r=='(')
    {
        ++s;
        int32_t opened_brackets = 1;
        bool nested = true;
        while(s-r!=len)
        {
            if(*s=='(')
            {
               if (opened_brackets<0)
               {
                    kstring_t s = {0,0,0};
                    kputsn(r, len, &s);
                    fprintf(stderr, "[%s:%d %s] brackets not correct %s\n", __FILE__, __LINE__, __FUNCTION__, s.s);
                    if (s.m) free(s.s);
                    exit(1);   
               }

               ++opened_brackets;
            }
            else if (*s==')')
            {
               if (opened_brackets<=0)
               {
                    kstring_t s = {0,0,0};
                    kputsn(r, len, &s);
                    fprintf(stderr, "[%s:%d %s] brackets not correct %s\n", __FILE__, __LINE__, __FUNCTION__, s.s);
                    if (s.m) free(s.s);
                    exit(1);
               }

               --opened_brackets;
            }

            if (opened_brackets==0)
            {
                r = s;
                
                if (s-r==len-1)
                {    
                    return -1;
                }
                else
                {
                    return 0;
                }
            }
            
            ++s;
        }
        
        kstring_t s = {0,0,0};
        kputsn(r, len, &s);
        fprintf(stderr, "[%s:%d %s] brackets not correct %s\n", __FILE__, __LINE__, __FUNCTION__, s.s);
        if (s.m) free(s.s);
        exit(1);
    }
    
    return 0;
}

/**
 * Parse literals.
 */
void Filter::parse_literal(const char* exp, int32_t len, Node * node, bool debug)
{
    if (strncmp(exp, "PASS", len)==0)
    {
        node->type = VT_FILTER_OP;
        kputsn(exp, len, &node->tag);
        if (debug) std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strncmp(exp, "VTYPE", len)==0)
    {
        node->type = VT_VARIANT_TYPE_OP;
        if (debug) std::cerr << "\tis variant_op\n";
        return;
    }
    else if (strncmp(exp, "N_ALLELE", len)==0)
    {
        node->type = VT_N_ALLELE_OP;
        if (debug) std::cerr << "\tis nallele_op\n";
        return;
    }
    else if (strncmp(exp, "INDEL", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_INDEL;
        if (debug) std::cerr << "\tis INDEL\n";
        return;
    }
    else if (strncmp(exp, "SNP", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_SNP;
        if (debug) std::cerr << "\tis SNP\n";
        return;
    }
    else if (strncmp(exp, "MNP", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_MNP;
        if (debug) std::cerr << "\tis MNP\n";
        return;
    }
    else if (strncmp(exp, "CLUMPED", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_CLUMPED;
        if (debug) std::cerr << "\tis CLUMPED\n";
        return;
    }
    else if (strncmp(exp, "DLEN", len)==0)
    {
        node->type = VT_VARIANT_DLEN_OP;
        if (debug) std::cerr << "\tis dlen\n";
        return;
    }
    else
    {
        const char* start = exp;
        char *end = NULL;
        int32_t i = std::strtoul(exp, &end, 10);
        if (end!=start)
        {
            node->type = VT_INT;
            node->i = i;
            if (debug) std::cerr << "\tis int\n";
            return;
        }

        start = exp;
        end = NULL;
        int32_t f = std::strtod(exp, &end);
        if (end!=start)
        {
            node->type = VT_FLT;
            node->f = f;
            if (debug) std::cerr << "\tis float\n";
            return;
        }
    }

    return;

}

bool Filter::is_literal(const char* exp, int32_t len)
{
    const char* q = exp;
    while (q-exp<len)
    {
        if(*q=='=' ||
           *q=='&' ||
           *q=='|' ||
           *q=='>' ||
           *q=='<' ||
           *q=='!')
        {
            return false;
        }
        
        ++q;
    }
    
    return true;
}

/**
 * Returns -1 if no operator found. Updates oplen to be the length of the operator observed.
 */
int32_t Filter::peek_op(const char* &r, int32_t len, int32_t &oplen, bool debug)
{
    const char* s = r;
    oplen = -1;
    if (*s=='&' && (s+1-r<len) && *(s+1)=='&')
    {
        if (debug) std::cerr << "\tis && operator\n";
        oplen = 2;
        return VT_OP_AND;
    }
    else if (*s=='|' && (s+1-r<len) && *(s+1)=='|')
    {
        if (debug) std::cerr << "\tis || operator\n";
        oplen = 2;
        return VT_OP_OR;
    }
    else if (*s=='=' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis == operator\n";
        oplen = 2;
        return VT_OP_EQ;
    }
    else if (*s=='>' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis >= operator\n";
        oplen = 2;
        return VT_OP_GE;
    }
    else if (*s=='<' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis <= operator\n";
        oplen = 2;
        return VT_OP_LE;
    }
    else if (*s=='!' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis != operator\n";
        oplen = 2;
        return VT_OP_NE;
    }
    else if (*s=='&')
    {
        if (debug) std::cerr << "\tis & operator\n";
        oplen = 1;
        return VT_OP_BIT_AND;
    }
    else if (*s=='|')
    {
        if (debug) std::cerr << "\tis | operator\n";
        oplen = 1;
        return VT_OP_BIT_OR;
    }
    else if (*s=='>')
    {
        if (debug) std::cerr << "\tis > operator\n";
        oplen = 1;
        return VT_OP_GT;
    }
    else if (*s=='<')
    {
        if (debug) std::cerr << "\tis < operator\n";
        oplen = 1;
        return VT_OP_LT;
    }
    
    return -1;
}

/**
 * Constructs the expression tree.
 */
void Filter::parse(const char* exp, int32_t len, Node *node, bool debug)
{
    if (debug)
    {
        std::cerr << "parsing \"";
        for (int32_t i=0; i<len; ++i)
            std::cerr << exp[i] ;
        std::cerr << "\" " << len << "\n";
    }

    //******************************
    //trim white spaces and brackets
    //******************************
    while (*exp==' ') ++exp;
    while (exp[len-1]==' ') --len;
    trim_brackets(exp, len);

    if (debug)
    {
//        std::cerr << "\tafter trimming white spaces and brackets: \"";
//        for (int32_t i=0; i<len; ++i)
//            std::cerr << exp[i] ;
//        std::cerr << "\" " << len << "\n";
    }

    //this is a literal
    if (is_literal(exp, len))
    {
        //will not recurse any further
        return parse_literal(exp, len, node, debug);
    }
    //this is guaranteed to be decomposed unless there is an error in the expression
    else
    {
        const char* p = exp; //points to end of first part
        const char* q = exp; //points to start of second part
        const char* r = exp; //for iteration

        int32_t type = INT_MAX;
            
        while(r-exp!=len)
        {
            fwd_to_closing_bracket(r, len);
            
            int32_t oplen = -1;
            int32_t ctype = peek_op(r, len, oplen, debug);
            
            if(ctype!=-1)
            {
                if (ctype<type)
                {
                    std::cerr<< "\tupdating type\n";
                    type = ctype;
                    p = r-1;
                    q = r+oplen;
                }   
                
                r += oplen-1;
            }
            
            ++r;
        }

        if (type==-1)
        {
            kstring_t s = {0,0,0};
            kputsn(exp, len, &s);
            fprintf(stderr, "[%s:%d %s] expression not correct %s\n", __FILE__, __LINE__, __FUNCTION__, s.s);
            if (s.m) free(s.s);
            exit(1);
        }
        
        std::cerr << "creating nodes\n";
        
        node->type = type;
        
        node->left = new Node();
        parse(exp, p-exp+1, node->left, debug);
        
        node->right = new Node();
        parse(q, len-(q-exp), node->right, debug);
    }
}

/**
 * Applies filter to vcf record.
 */
bool Filter::apply(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug) //recursive
{
    if (tree==NULL)
    {
        return true;
    }

    this->h = h;
    this->v = v;
    this->variant = variant;

    if (debug) std::cerr << "==========\n";
    apply(tree, debug);
    if (debug) std::cerr << "==========\n";

    return tree->value;
}

/**
 * Recursive call for apply.
 */
void Filter::apply(Node* node, bool debug)
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

    node->evaluate(h, v, variant, debug);
}


/**
 * Evaluates the actions for this node.
 */
void Node::evaluate(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug)
{
    if (debug)
        std::cerr << "evaluation  "  << type << "\n";

    //filter operation
    if (type==VT_OP_AND)
    {
        if (debug)
            std::cerr << "\tVT_OP_AND "   <<  left->value << "&" << right->value    <<  " \n";
        value = (left->value && right->value);
    }
    else if (type==VT_OP_OR)
    {
        value = (left->value || right->value);
    }
    else if (type==VT_OP_BIT_AND)
    {
        if (debug)
            std::cerr << "\tVT_OP_BIT_AND "   <<  left->i << "&" << right->i << "=" << (left->i&right->i)    <<  " \n";
        i = (left->i & right->i);
        value = i;
    }
    else if (type==VT_OP_BIT_OR)
    {
        i = (left->i | right->i);
        value = i;
    }
    else if (type==VT_FILTER_OP)
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
        value = i;
    }
    else if (type==VT_VARIANT_DLEN_OP)
    {
        i = variant->alleles[0].dlen;
        value = i;

      //  std::cerr << "DLEN ASSIGN: " << i << "\n";
    }
    else if (type==VT_N_ALLELE_OP)
    {
        i = bcf_get_n_allele(v);
    }
    else if (type==VT_OP_EQ)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            value = (left->i==right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            value = (left->f==right->f);
        }
        else if ((left->type&VT_STR) && (right->type&VT_STR))
        {
            value = strcmp(left->tag.s, right->tag.s)==0 ? true : false;
        }
        else
        {
            std::cerr << "evaluation not supported, different types to compare\n";
        }
    }
    else if (type==VT_OP_NE)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            value = (left->i!=right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            value = (left->f!=right->f);
        }
        else if ((left->type&VT_STR) && (right->type&VT_STR))
        {
            value = strcmp(left->tag.s, right->tag.s)!=0 ? true : false;
        }
        else
        {
            std::cerr << "evaluation not supported, different types to compare\n";
        }
    }
    else if (type==VT_OP_LE)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            value = (left->i<=right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            value = (left->f<=right->f);
        }
        else if ((left->type&VT_STR) && (right->type&VT_STR))
        {
            value = strcmp(left->tag.s, right->tag.s) ? true : false;
        }
        else
        {
            std::cerr << "evaluation not supported, different types to compare\n";
        }
    }
    else if (type==VT_OP_GE)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            value = (left->i>=right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            value = (left->f>=right->f);
        }
        else if ((left->type&VT_STR) && (right->type&VT_STR))
        {
            value = strcmp(left->tag.s, right->tag.s)>=0 ? true : false;
        }
        else
        {
            std::cerr << "evaluation not supported, different types to compare\n";
        }
    }
    else if (type==VT_OP_GT)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            value = (left->i>right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            value = (left->f>right->f);
        }
        else if ((left->type&VT_STR) && (right->type&VT_STR))
        {
            value = strcmp(left->tag.s, right->tag.s)>0 ? true : false;
        }
        else
        {
            std::cerr << "evaluation not supported, different types to compare\n";
        }
    }
    else if (type==VT_OP_LT)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            value = (left->i<right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            value = (left->f<right->f);
        }
        else if ((left->type&VT_STR) && (right->type&VT_STR))
        {
            value = strcmp(left->tag.s, right->tag.s)<0 ? true : false;
        }
        else
        {
            std::cerr << "evaluation not supported, different types to compare\n";
        }
    }
}
