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
                    std::cerr << "Illegal expression: brackets not correct\n";
               }

               ++opened_brackets;
            }
            else if (exp[j]==')')
            {
               if (opened_brackets<=0)
               {
                    std::cerr << "Illegal expression: brackets not correct\n";
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

    //*****************
    //trim white spaces
    //*****************
    while (*exp==' ') ++exp;
    while (exp[len-1]==' ') --len;

    //*************
    //trim brackets
    //*************
    trim_brackets(exp, len);

    if (debug)
    {
        std::cerr << "\tafter trimming white spacecs and brackets: \"";
        for (int32_t i=0; i<len; ++i)
            std::cerr << exp[i] ;
        std::cerr << "\" " << len << "\n";
    }

    //this is a literal
    if (is_literal(exp, len))
    {
        parse_literal(exp, len, node, debug);
    }
    else
    {
        const char* p = exp;
        const char* q = exp;

        int32_t type = -1;

        //handle brackets
        if (*p=='(')
        {
            int32_t opened_brackets = 1;
            while(q-p!=len)
            {
                if(*q=='(')
                {
                   if (opened_brackets<0)
                   {
                        std::cerr << "Illegal expression: brackets not correct\n";
                   }

                   ++opened_brackets;
                }
                else if (*q==')')
                {
                   if (opened_brackets<=0)
                   {
                        std::cerr << "Illegal expression: brackets not correct\n";
                   }

                   --opened_brackets;
                }

                ++q;

                if (opened_brackets==0)
                {
                    break;
                }
            }
        }
        else
        {
            int32_t no_ops = 0;
            
            //implements left hand associativity
            while(q-p!=len)
            {
                if((*q=='&' && (q-p<len) && *(q+1)=='&') ||
                   (*q=='|' && (q-p<len) && *(q+1)=='|') ||
                   (*q=='=' && (q-p<len) && *(q+1)=='=') ||
                   (*q=='!' && (q-p<len) && *(q+1)=='=') ||
                   (*q=='>' && (q-p<len) && *(q+1)=='=') ||
                   (*q=='<' && (q-p<len) && *(q+1)=='=') ||
                   *q=='&' ||
                   *q=='|' ||
                   *q=='>' ||
                   *q=='<' ||
                   *q=='!')
                {
                    if (no_ops==0)
                    {
                        ++no_ops;
                    }
                    else if (no_ops==1)
                    {
                        ++q;
                        break;    
                    }
                }
                
                ++q;
            }
        }

        int32_t oplen = 0;
        while (q-p<len)
        {
            if (*q=='&' && (q-exp<len) && *(q+1)=='&')
            {
                type = VT_OP_AND;
                oplen = 2;
                if (debug) std::cerr << "\tis && operator\n";
                break;
            }
            else if (*q=='&')
            {
                type = VT_OP_BIT_AND;
                oplen = 1;
                if (debug) std::cerr << "\tis & operator\n";
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
                if (debug) std::cerr << "\tis | operator\n";
                break;
            }
            else if (*q=='=' && (q-exp<len) && *(q+1)=='=')
            {
                type = VT_OP_EQ;
                oplen = 2;
                if (debug) std::cerr << "\tis == operator\n";
                break;
            }
            else if (*q=='>' && (q-exp<len) && *(q+1)=='=')
            {
                type = VT_OP_GE;
                oplen = 2;
                if (debug) std::cerr << "\tis >= operator\n";
                break;
            }
            else if (*q=='<' && (q-exp<len) && *(q+1)=='=')
            {
                type = VT_OP_LE;
                oplen = 2;
                if (debug) std::cerr << "\tis <= operator\n";
                break;
            }
            else if (*q=='!' && (q-exp<len) && *(q+1)=='=')
            {
                type = VT_OP_NE;
                oplen = 2;
                if (debug) std::cerr << "\tis != operator\n";
                break;
            }
            else if (*q=='>')
            {
                type = VT_OP_GT;
                oplen = 1;
                if (debug) std::cerr << "\tis > operator\n";
                break;
            }
            else if (*q=='<')
            {
                type = VT_OP_LT;
                oplen = 1;
                if (debug) std::cerr << "\tis < operator\n";
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
        parse(p, q-p, node->left, debug);

        if (q-p!=len)
        {
            q += oplen;
    //        std::cerr << "\tparsing q \"";
    //        for (int32_t i=0; i<len-(q-p); ++i)
    //            std::cerr << q[i] ;
    //        std::cerr << "\" " << (len-(q-p)) << "\n";

            node->right = new Node();
            parse(q, len-(q-exp), node->right, debug);
        }
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
