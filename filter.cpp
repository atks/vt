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

    if (type&VT_LOGIC_OP)
    {
        if (type==VT_NOT)
        {
            if (debug)
                std::cerr << "\tVT_NOT "   <<  left->value << " \n";
            value = !(left->value);
        }
        else if (type==VT_AND)
        {
            if (debug)
                std::cerr << "\tVT_AND "   <<  left->value << "&" << right->value    <<  " \n";
            value = (left->value && right->value);
        }
        else if (type==VT_OR)
        {
            value = (left->value || right->value);
        }
    }
    else if (type&VT_MATH_CMP)   
    {
        if (type==VT_EQ)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->i << "&" << right->i    <<  " \n";
                    value = (left->i==right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->i << "&" << right->f    <<  " \n";
                    value = (left->i==right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->f << "&" << right->i    <<  " \n";
                    value = (left->f==right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->f << "&" << right->f    <<  " \n";
                    value = (left->f==right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->tag.s << "&" << right->tag.s    <<  " \n";
                value = strcmp(left->tag.s, right->tag.s)==0 ? true : false;
                return;
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported : == %d %d\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }    
        else if (type==VT_NE)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->i!=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->i!=right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->f!=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->f!=right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                value = strcmp(left->tag.s, right->tag.s)==0 ? false : true;
                return;
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported: %d %d: !=\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
        else if (type==VT_LE)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->i<=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->i<=right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    value = (left->f<=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->f<=right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                value = strcmp(left->tag.s, right->tag.s)<=0 ? true : false;
                return;
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported: %d %d: <=\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
        else if (type==VT_GE)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->i>=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->i>=right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->f>=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->f>=right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                value = strcmp(left->tag.s, right->tag.s)>=0 ? true : false;
                return;
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported: %d %d: >=\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
        else if (type==VT_GT)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->i>right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->i>right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->f>right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->f>right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                value = strcmp(left->tag.s, right->tag.s)>0 ? true : false;
                return;
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported: %d %d: >\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
        else if (type==VT_LT)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->i<right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->i<right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    value = (left->f<right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    value = (left->f<right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                value = strcmp(left->tag.s, right->tag.s)<0 ? true : false;
                return;
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported: %d %d: <\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
    }
    else if (type&VT_BCF_OP)   
    {
        if (type==VT_FILTER)
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
        else if (type==VT_INFO)
        {
            int32_t *data = NULL;
            int32_t n=0;
    
            if (bcf_get_info_int32(h, v, tag.s, &data, &n)>0)
            {
                type |= VT_INT;
                i = *data;
                f = (float)i;
            }
            else if (bcf_get_info_float(h, v, tag.s, &data, &n)>0)
            {
                type |= VT_FLT;
                f = (float)(*data);
            }
            else if (bcf_get_info_string(h, v, tag.s, &data, &n)>0)
            {
                type |= VT_STR;
                s.l=0;
                for (int32_t i=0; i<n; ++i)
                {
                    kputc(data[i], &s);
                }
            }
            else if (bcf_get_info_flag(h, v, tag.s, 0, 0)>0)
            {
                type |= VT_FLG;
                i = 1;
                f = 1;
                b = true;
                value = true;
                s.l=0; 
            }
            else
            {
                i = 0;
                f = 0;
                b = false;
                value = false;
                s.l=0;
            }
    
            if (n) free(data);
        }
        else if (type==(VT_INFO|VT_INT))
        {
            int32_t *data = NULL;
            int32_t n=0;
    
            if (bcf_get_info_int32(h, v, tag.s, &data, &n)>0)
            {
                i = *((int*)data);
            }
    
            if (n) free(data);
        }
        else if (type==(VT_INFO|VT_FLT))
        {
            int32_t *data = NULL;
            int32_t n=0;
    
            if (bcf_get_info_float(h, v, tag.s, &data, &n)>0)
            {
                f = *((float*)data);
            }
    
            if (n) free(data);
        }
        else if (type==(VT_INFO|VT_STR))
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
        else if (type==(VT_INFO|VT_FLG))
        {
            if (bcf_get_info_flag(h, v, tag.s, 0, 0)>0)
            {
                i = 1;
                f = 1;
                b = true;
                value = true;
                //s.l=0; kputc('1', &s);
            }
            else
            {
                i = 0;
                f = 0;
                b = false;
                value = false;
                s.l=0;
            }
            
            if (debug)
                std::cerr << "\tVT_INFO|VT_FLG "   << i << " " << f << " " << b << " " << value << " " << s.s <<  " \n";
        }
        else if (type==VT_VARIANT_TYPE)
        {
            if (debug)
                std::cerr << "\tVTYPE "   <<  variant->vtype2string(variant->type) <<  " \n";
            i = variant->type;
            value = i;
        }
        else if (type==VT_VARIANT_DLEN)
        {
            if (debug)
                std::cerr << "\tDLEN "   <<  variant->alleles[0].dlen <<  " \n";
            i = variant->alleles[0].dlen;
            value = i;
        }
        else if (type==VT_VARIANT_LEN)
        {
            if (debug)
                std::cerr << "\tLEN "   <<  abs(variant->alleles[0].dlen) <<  " \n";
            i = abs(variant->alleles[0].dlen);
            value = i;
        }
        else if (type==VT_N_ALLELE)
        {
            if (debug)
                std::cerr << "\tN_ALLELE "   <<  bcf_get_n_allele(v) <<  " \n";
            i = bcf_get_n_allele(v);
        }
    }
    else if (type&VT_MATH_OP)
    {   
        if ((type&8207)==VT_ADD)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    i = (left->i+right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    type |= VT_FLT;
                    f = (left->i+right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_FLT;
                    f = (left->f+right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    type |= VT_FLT;
                    f = (left->f+right->f);
                    return;
                }
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported : +\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
        else if ((type&8207)==VT_SUB)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    i = (left->i-right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    type |= VT_FLT;
                    f = (left->i-right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_FLT;
                    f = (left->f-right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    type |= VT_FLT;
                    f = (left->f-right->f);
                    return;
                }
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported : -\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
        else if ((type&8207)==VT_MUL)
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    i = (left->i*right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    type |= VT_FLT;
                    f = (left->i*right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_FLT;
                    f = (left->f*right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    type |= VT_FLT;
                    f = (left->f*right->f);
                    return;
                }
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported : *\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
        else if ((type&8207)==VT_DIV)
        {
            if (left->type&VT_INT)
            {
                if (right->type&VT_INT)
                {
                    type |= VT_FLT;
                    f = ((float)left->i/right->i);
                    return;
                }
                else if (right->type&VT_FLT)
                {
                    type |= VT_FLT;
                    f = (left->i/right->f);
                    return;
                }
            }
            else if (left->type&VT_FLT)
            {
                if (right->type&VT_INT)
                {
                    type |= VT_FLT;
                    f = (left->f/right->i);
                    return;
                }
                else if (right->type&VT_FLT)
                {
                    type |= VT_FLT;
                    f = (left->f/right->f);
                    return;
                }
            }
    
            fprintf(stderr, "[%s:%d %s] evaluation not supported : /\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
        else if (type==VT_BIT_AND)
        {
            if ((left->type&VT_INT) && (right->type&VT_INT))
            {
                i = (left->i & right->i);
                value = i;
                return;
            }
            
            fprintf(stderr, "[%s:%d %s] evaluation not supported for & :  %d %d\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
        else if (type==VT_BIT_OR)
        {
            if ((left->type&VT_INT) && (right->type&VT_INT))
            {
                i = (left->i | right->i);
                value = i;
                return;
            }
            
            fprintf(stderr, "[%s:%d %s] evaluation not supported for | : %d %d\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }    
        else
        {
            fprintf(stderr, "[%s:%d %s] math op not supported : %d\n", __FILE__, __LINE__, __FUNCTION__, (type&15));
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
Filter::Filter(std::string exp, VariantManip *vm)
{
    this->tree = NULL;
    parse(exp.c_str(), false);
};

/**
 * Parses filter expression.
 */
void Filter::parse(const char* exp, bool debug)
{
    need_to_classify_variant = false;
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
    node->type = VT_UNKNOWN;
    
    // if not sign in front of literal
    if (exp[0]=='~')
    {
        node->type = VT_NOT;
        node->left = new Node();
        if (debug) std::cerr << "\tis not_op\n";

        node = node->left;
        ++exp;
        --len;
    }

    if (strncmp(exp, "PASS", len)==0)
    {
        node->type = VT_FILTER;
        kputsn(exp, len, &node->tag);
        if (debug) std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strncmp(exp, "FILTER.", 7)==0)
    {
        node->type = VT_FILTER;
        exp += 7;
        kputsn(exp, len-7, &node->tag);
        if (debug) std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strncmp(exp, "INFO.", 5)==0)
    {
        node->type = VT_INFO;
        exp += 5;
        kputsn(exp, len-5, &node->tag);
        if (debug) std::cerr << "\tis info_op\n";
        return;
    }
    else if (strncmp(exp, "VTYPE", 5)==0)
    {
        need_to_classify_variant = true;
        node->type = VT_VARIANT_TYPE;
        if (debug) std::cerr << "\tis variant_op\n";
        return;
    }
    else if (strncmp(exp, "N_ALLELE", len)==0)
    {
        node->type = VT_N_ALLELE;
        if (debug) std::cerr << "\tis nallele_op\n";
        return;
    }
    else if (strncmp(exp, "INDEL", len)==0)
    {
        need_to_classify_variant = true;
        node->type = VT_INT;
        node->i = VT_INDEL;
        if (debug) std::cerr << "\tis INDEL\n";
        return;
    }
    else if (strncmp(exp, "SNP", len)==0)
    {
        need_to_classify_variant = true;
        node->type = VT_INT;
        node->i = VT_SNP;
        if (debug) std::cerr << "\tis SNP\n";
        return;
    }
    else if (strncmp(exp, "MNP", len)==0)
    {
        need_to_classify_variant = true;
        node->type = VT_INT;
        node->i = VT_MNP;
        if (debug) std::cerr << "\tis MNP\n";
        return;
    }
    else if (strncmp(exp, "CLUMPED", len)==0)
    {
        need_to_classify_variant = true;
        node->type = VT_INT;
        node->i = VT_CLUMPED;
        if (debug) std::cerr << "\tis CLUMPED\n";
        return;
    }
    else if (strncmp(exp, "DLEN", len)==0)
    {
        need_to_classify_variant = true;
        node->type = VT_VARIANT_DLEN;
        if (debug) std::cerr << "\tis dlen\n";
        return;
    }
    else if (strncmp(exp, "LEN", len)==0)
    {
        need_to_classify_variant = true;
        node->type = VT_VARIANT_LEN;
        if (debug) std::cerr << "\tis len\n";
        return;
    }
    else
    {
        //integer type
        const char* start = exp;
        char *end = NULL;
        int32_t i = std::strtoul(exp, &end, 10);
        if (end==exp+len)
        {
            node->type = VT_INT;
            node->i = i;
            node->f = (float)i;
            if (debug) std::cerr << "\tis int\n";
            return;
        }

        //float type
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
    
        //string type
        if (exp[0]=='"' && exp[len-1]=='"')
        {
            node->type = VT_STR;
            kputsn(exp, len, &node->tag);
            if (debug) std::cerr << "\tis string\n";
            return;
        }
        
        if (node->type==VT_UNKNOWN)
        {
            kstring_t s = {0,0,0}; 
            kputsn(exp, len, &s);
            fprintf(stderr, "[E:%s] %s is not recognised\n", __FUNCTION__, s.s);
            print_filter_help();
            exit(1);
        }
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
 * Help message on filter expressions.
 */
void Filter::print_filter_help()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "  Variant characteristics\n");
    fprintf(stderr, "    VTYPE,N_ALLELE,DLEN,LEN\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Variant value types\n");
    fprintf(stderr, "    SNP,MNP,INDEL,CLUMPED\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Biallelic SNPs only                       : VTYPE==SNP&&N_ALLELE==2\n");
    fprintf(stderr, "  Biallelic Indels with embedded SNP        : VTYPE==(SNP|INDEL)&&N_ALLELE==2\n");
    fprintf(stderr, "  Biallelic variants involving insertions   : VTYPE&INDEL&&DLEN>0&&N_ALLELE==2\n");
    fprintf(stderr, "  Biallelic variants involving 1bp variants : LEN==1&&N_ALLELE==2\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  FILTER fields\n");
    fprintf(stderr, "    PASS, FILTER.<tag>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  INFO fields\n");
    fprintf(stderr, "    INFO.<tag>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Passed biallelic SNPs only                  : PASS&&VTYPE==SNP&&N_ALLELE==2\n");
    fprintf(stderr, "  Passed Common biallelic SNPs only           : PASS&&VTYPE==SNP&&N_ALLELE==2&&INFO.AF>0.005\n");
    fprintf(stderr, "  Passed Common biallelic SNPs or rare indels : (PASS&&VTYPE==SNP&&N_ALLELE==2&&INFO.AF>0.005)||(VTYPE&INDEL&&INFO.AF<=0.005)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Operations\n");
    fprintf(stderr, "    ==,~,&&,||,&,|,+,-,*,/\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Failed rare variants : ~PASS&&(INFO.AC/INFO.AN<0.005)\n");
    fprintf(stderr, "\n");
    
}

/**
 * Trim brackets from an expression.
 */
void Filter::trim_brackets(const char* &exp, int32_t &len)
{
    if (*exp=='(' && exp[len-1]==')')
    {
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
        return VT_AND;
    }
    else if (*s=='|' && (s+1-r<len) && *(s+1)=='|')
    {
        if (debug) std::cerr << "\tis || operator\n";
        oplen = 2;
        return VT_OR;
    }
    else if (*s=='=' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis == operator\n";
        oplen = 2;
        return VT_EQ;
    }
    else if (*s=='>' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis >= operator\n";
        oplen = 2;
        return VT_GE;
    }
    else if (*s=='<' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis <= operator\n";
        oplen = 2;
        return VT_LE;
    }
    else if (*s=='!' && (s+1-r<len) && *(s+1)=='=')
    {
        if (debug) std::cerr << "\tis != operator\n";
        oplen = 2;
        return VT_NE;
    }
    else if (*s=='+')
    {
        if (debug) std::cerr << "\tis + operator\n";
        oplen = 1;
        return VT_ADD;
    }
    else if (*s=='-')
    {
        if (debug) std::cerr << "\tis - operator\n";
        oplen = 1;
        return VT_SUB;
    }
    else if (*s=='*')
    {
        if (debug) std::cerr << "\tis * operator\n";
        oplen = 1;
        return VT_MUL;
    }
    else if (*s=='/')
    {
        if (debug) std::cerr << "\tis / operator\n";
        oplen = 1;
        return VT_DIV;
    }
    else if (*s=='&')
    {
        if (debug) std::cerr << "\tis & operator\n";
        oplen = 1;
        return VT_BIT_AND;
    }
    else if (*s=='|')
    {
        if (debug) std::cerr << "\tis | operator\n";
        oplen = 1;
        return VT_BIT_OR;
    }
    else if (*s=='>')
    {
        if (debug) std::cerr << "\tis > operator\n";
        oplen = 1;
        return VT_GT;
    }
    else if (*s=='<')
    {
        if (debug) std::cerr << "\tis < operator\n";
        oplen = 1;
        return VT_LT;
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
        apply(node->left, debug);
    }

    //can do some lazy evaluation here for && and || types.
    if (node->type&VT_LOGIC_OP)
    {
        if (node->type==VT_AND)
        {
            if (!node->left->value)
            {
                node->value = false;
                return;
            }
            else
            {
                 apply(node->right, debug);
            }
        }
        else if (node->type==VT_OR)
        {
            if (node->left->value)
            {
                node->value = true;
                return;
            }
            else
            {
                 apply(node->right, debug);
            }
        }
    }
    else
    {
        if (node->right!=NULL)
        {
            apply(node->right, debug);
        }
    }
   
    node->evaluate(h, v, variant, debug);
}
