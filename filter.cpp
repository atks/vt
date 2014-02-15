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

/**
 * Constructor.
 */
Node::Node()
{
    parent = NULL;
    left = NULL;
    right = NULL;
    tag = {0,0,0};
};

/**
 * Constructor with type initlialization.
 */
Node::Node(int32_t type)
{
    parent = NULL;
    left = NULL;
    right = NULL;
    tag = {0,0,0};
    this->type = type;
};

/**
 * Evaluates the actions for this node.
 */
void Node::evaluate(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug)
{
    if (debug)
        std::cerr << "evaluation  "  << type << "\n";

    //filter operation
    if (type==VT_OP_NOT)
    {
        if (debug)
            std::cerr << "\tVT_OP_NOT "   <<  left->value << " \n";
        value = !(left->value);
    }
    else if (type==VT_OP_AND)
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
        {
            std::cerr << "\tVT_OP_BIT_AND "   <<  left->i << "&" << right->i << "=" << (left->i&right->i)    <<  " \n";
        }
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
    else if (type==VT_INFO_OP)
    {
        int32_t *data = NULL;
        int32_t n=0;

        if (bcf_get_info_int(h, v, tag.s, &data, &n)>0)
        {
            std::cerr << "INFO INT\n";
            type = VT_INFO_INT_OP;
            i = *data;
        }
        else if (bcf_get_info_float(h, v, tag.s, &data, &n)>0)
        {
            type = VT_INFO_FLT_OP;
            f = (float)(*data);
        }
        else if (bcf_get_info_string(h, v, tag.s, &data, &n)>0)
        {
            type = VT_INFO_STR_OP;
            s.l=0;
            for (int32_t i=0; i<n; ++i)
            {
                kputc(data[i], &s);
            }
        }

        if (n) free(data);
    }
    else if (type==VT_INFO_INT_OP)
    {
        int32_t *data = NULL;
        int32_t n=0;

        if (bcf_get_info_int(h, v, tag.s, &data, &n)>0)
        {
            i = *((int*)data);
        }

        if (n) free(data);
    }
    else if (type==VT_INFO_FLT_OP)
    {
        int32_t *data = NULL;
        int32_t n=0;

        if (bcf_get_info_float(h, v, tag.s, &data, &n)>0)
        {
            f = *((float*)data);

        //    std::cerr << "FLOAT " << f << "\n";
        }

        if (n) free(data);
    }
    else if (type==VT_INFO_STR_OP)
    {
        int32_t *data = NULL;
        int32_t n=0;

        if (bcf_get_info_string(h, v, tag.s, &data, &n)>0)
        {
            s.l=0;
            for (int32_t i=0; i<n; ++i)
            {
                kputc(data[i], &s);
            }
        }

        if (n) free(data);
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
    }
    else if (type==VT_VARIANT_LEN_OP)
    {
        i = abs(variant->alleles[0].dlen);
        value = i;
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
            fprintf(stderr, "[%s:%d %s] evaluation not supported : ==\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
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
            fprintf(stderr, "[%s:%d %s] evaluation not supported : !=\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
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
            fprintf(stderr, "[%s:%d %s] evaluation not supported : <=\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
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
            fprintf(stderr, "[%s:%d %s] evaluation not supported : >=\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
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
            fprintf(stderr, "[%s:%d %s] evaluation not supported : >\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
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


           // std::cerr << left->f << "<" << right->f << " " <<value << "\n";

        }
        else if ((left->type&VT_STR) && (right->type&VT_STR))
        {
            value = strcmp(left->tag.s, right->tag.s)<0 ? true : false;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] evaluation not supported : < %d %d\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
    }
    else if (type==VT_OP_ADD)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            f = (left->i+right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            f = (left->f+right->f);
        }
        else if ((left->type&VT_INT) && (right->type&VT_FLT))
        {
            f = (left->i+right->f);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_INT))
        {
            f = (left->f+right->i);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] evaluation not supported : +\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
    }
    else if (type==VT_OP_SUB)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            f = (left->i-right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            f = (left->f-right->f);
        }
        else if ((left->type&VT_INT) && (right->type&VT_FLT))
        {
            f = (left->i-right->f);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_INT))
        {
            f = (left->f-right->i);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] evaluation not supported : -\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
    }
    else if (type==VT_OP_MUL)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            f = (left->i*right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            f = (left->f*right->f);
        }
        else if ((left->type&VT_INT) && (right->type&VT_FLT))
        {
            f = (left->i*right->f);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_INT))
        {
            f = (left->f*right->i);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] evaluation not supported : *\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
    }
    else if (type==VT_OP_DIV)
    {
        if ((left->type&VT_INT) && (right->type&VT_INT))
        {
            f = (left->i/right->i);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_FLT))
        {
            f = (left->f/right->f);
        }
        else if ((left->type&VT_INT) && (right->type&VT_FLT))
        {
            f = (left->i/right->f);
        }
        else if ((left->type&VT_FLT) && (right->type&VT_INT))
        {
            f = (left->f/right->i);
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] evaluation not supported : /\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
    }
}

/**
 * Constructor.
 */
Filter::Filter()
{
    this->tree = NULL;
};

/**
 * Constructor with expression initialization.
 */
Filter::Filter(std::string exp)
{
    parse(exp.c_str(), false);
};

/**
 * Parses filter expression.
 */
void Filter::parse(const char* exp, bool debug)
{
    if (strlen(exp)!=0)
    {
        if (tree!=NULL)
        {
            delete tree;
            tree = NULL;
        }
        else
        {
            tree = new Node();
            parse(exp, strlen(exp), tree, debug);
        }
    }
    else
    {
        tree = NULL;
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
                    if (debug) std::cerr<< "\tupdating type\n";
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

        node->type = type;

        node->left = new Node();
        parse(exp, p-exp+1, node->left, debug);

        node->right = new Node();
        parse(q, len-(q-exp), node->right, debug);
    }
}

/**
 * Parse literals.
 */
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
           *q=='!' ||
           *q=='+' ||
           *q=='*' ||
           *q=='/' ||
           (*q=='-' && exp!=q))
        {
            return false;
        }

        ++q;
    }

    return true;
}

/**
 * Parse literals.
 */
void Filter::parse_literal(const char* exp, int32_t len, Node * node, bool debug)
{
    // if not sign in front of literal
    if (exp[0]=='~')
    {
        node->type = VT_OP_NOT;
        node->left = new Node();
        if (debug) std::cerr << "\tis not_op\n";

        node = node->left;
        ++exp;
        --len;
    }

    if (strncmp(exp, "PASS", len)==0)
    {
        node->type = VT_FILTER_OP;
        kputsn(exp, len, &node->tag);
        if (debug) std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strncmp(exp, "FILTER.", 7)==0)
    {
        node->type = VT_FILTER_OP;
        exp += 7;
        kputsn(exp, len-7, &node->tag);
        if (debug) std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strncmp(exp, "INFO.", 5)==0)
    {
        node->type = VT_INFO_OP;
        exp += 5;
        kputsn(exp, len-5, &node->tag);
        if (debug) std::cerr << "\tis info_op\n";
        return;
    }
    else if (strncmp(exp, "VTYPE", 5)==0)
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
    else if (strncmp(exp, "LEN", len)==0)
    {
        node->type = VT_VARIANT_LEN_OP;
        if (debug) std::cerr << "\tis len\n";
        return;
    }
    else
    {
        const char* start = exp;
        char *end = NULL;
        int32_t i = std::strtoul(exp, &end, 10);
        if (end==exp+len)
        {
            node->type = VT_INT;
            node->i = i;
            if (debug) std::cerr << "\tis int\n";
            return;
        }

        start = exp;
        end = NULL;
        float f = std::strtod(exp, &end);
        if (end==exp+len)
        {
            node->type = VT_FLT;
            node->f = f;
            if (debug) std::cerr << "\tis float: " << f << "\n";
            return;
        }

        kputsn(exp, len, &node->tag);
        if (debug) std::cerr << "\tis string\n";
        return;
    }

    if (debug)
    {
        std::cerr << "\tvalue " << node->value << "\n";
        std::cerr << "\ttag   " << node->tag.s << "\n";
        std::cerr << "\tb     " << node->b << "\n";
        std::cerr << "\ti     " << node->i << "\n";
        std::cerr << "\tf     " << node->f << "\n";
    }

    return;
}

/**
 * Trim brackets from an expression.
 */
void Filter::trim_brackets(const char* &exp, int32_t &len)
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
    else if (*s=='+')
    {
        if (debug) std::cerr << "\tis + operator\n";
        oplen = 1;
        return VT_OP_ADD;
    }
    else if (*s=='-')
    {
        if (debug) std::cerr << "\tis - operator\n";
        oplen = 1;
        return VT_OP_SUB;
    }
    else if (*s=='*')
    {
        if (debug) std::cerr << "\tis * operator\n";
        oplen = 1;
        return VT_OP_MUL;
    }
    else if (*s=='/')
    {
        if (debug) std::cerr << "\tis / operator\n";
        oplen = 1;
        return VT_OP_DIV;
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
