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
    s = {0,0,0};
    regex_set = false;
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
    s = {0,0,0};
    regex_set = false;
    this->type = type;
};

/**
 * Evaluates the actions for this node.
 */
void Node::evaluate(bcf_hdr_t *h, bcf1_t *v, Variant *variant, bool debug)
{
    if (debug)
        std::cerr << "evaluation  "  << type2string(type) << "\n";

    //by default
    value_exists = true;

    if (type&VT_LOGIC_OP)
    {
        if (!left->value_exists && type==VT_NOT)
        {
            if (debug)
                std::cerr << "\tVT_NOT "   <<  left->b << " \n";

            b = true;
            return;
        }

        //think about this, what happens when INFO values are not present, do you treat as a flag?
        if (!left->value_exists)
        {
            b = false;
            value_exists = false;
            return;
        }

        if (type==VT_NOT)
        {
            if (debug)
                std::cerr << "\tVT_NOT "   <<  left->b << " \n";

            b = !(left->b);
            return;
        }

        if (!right->value_exists)
        {
            b = false;
            value_exists = false;
            return;
        }

        if (type==VT_AND)
        {
            if (debug)
                std::cerr << "\tVT_AND "   <<  left->b << "&" << right->b    <<  " \n";
            b = (left->b && right->b);
            return;
        }
        else if (type==VT_OR)
        {
            b = (left->b || right->b);
            return;
        }
    }
    else if (type&VT_MATH_CMP)
    {
        if (!left->value_exists && !right->value_exists)
        {
            value_exists = false;
            return;
        }

        if (type==VT_EQ)
        {
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }
            
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->i << "&" << right->i    <<  " \n";
                    b = (left->i==right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->i << "&" << right->f    <<  " \n";
                    b = (left->i==right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->f << "&" << right->i    <<  " \n";
                    b = (left->f==right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    if (debug)
                        std::cerr << "\tVT_EQ "   <<  left->f << "&" << right->f    <<  " \n";
                    b = (left->f==right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                if (debug)
                    std::cerr << "\tVT_EQ "   <<  left->s.s << "&" << right->s.s    <<  " \n";
                b = strcmp(left->s.s, right->s.s)==0 ? true : false;
                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported : == %s %s\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
        else if (type==VT_MATCH)
        {
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }
            
            if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                if (debug)
                    std::cerr << "\tVT_MATCH "   <<  left->s.s << "&" << right->s.s    <<  " \n";

                if (!regex_set)
                {
                    pregex.set(right->s.s);
                    regex_set = true;
                }

                b = pregex.match(left->s.s);

                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported : =~ %s %s\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
        else if (type==VT_NO_MATCH)
        {
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }
            
            if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                if (debug)
                    std::cerr << "\tVT_NO_MATCH "   <<  left->s.s << "&" << right->s.s    <<  " \n";

                if (!regex_set)
                {
                    pregex.set(right->s.s);
                    regex_set = true;
                }

                b = !pregex.match(left->s.s);

                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported : ~~ %s %s\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
        else if (type==VT_NE)
        {
            //fprintf(stderr, "[%s:%d %s] check: %s %s: !=\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }

            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->i!=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->i!=right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->f!=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->f!=right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                b = strcmp(left->s.s, right->s.s)==0 ? false : true;
                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported: %s %s: !=\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
        else if (type==VT_LE)
        {
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }
            
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->i<=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->i<=right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    b = (left->f<=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->f<=right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                b = strcmp(left->s.s, right->s.s)<=0 ? true : false;
                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported: %s %s: <=\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
        else if (type==VT_GE)
        {
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }
            
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->i>=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->i>=right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->f>=right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->f>=right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                b = strcmp(left->s.s, right->s.s)>=0 ? true : false;
                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported: %s %s: >=\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
        else if (type==VT_GT)
        {
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }
            
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->i>right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->i>right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->f>right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->f>right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                b = strcmp(left->s.s, right->s.s)>0 ? true : false;
                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported: %s %s: >\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
        else if (type==VT_LT)
        {
            //if an INFO field is involved and does not exist, then we evaluate as false
            if (!left->value_exists || !right->value_exists)
            {
                value_exists = true;
                b = false;
                return;
            }
            
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->i<right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->i<right->f);
                    return;
                }
            }
            else if ((left->type&VT_FLT))
            {
                if ((right->type&VT_INT))
                {
                    b = (left->f<right->i);
                    return;
                }
                else if ((right->type&VT_FLT))
                {
                    b = (left->f<right->f);
                    return;
                }
            }
            else if ((left->type&VT_STR) && (right->type&VT_STR))
            {
                b = strcmp(left->s.s, right->s.s)<0 ? true : false;
                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported: %s %s: <\n", __FILE__, __LINE__, __FUNCTION__, type2string(left->type).c_str(), type2string(right->type).c_str());
            exit(1);
        }
    }
    else if (type&VT_BCF_OP)
    {
        value_exists = true;

        if (type==VT_QUAL)
        {
            if (bcf_float_is_missing(bcf_get_qual(v)))
            {
                f = -1;
                value_exists = false;
            }
            else
            {
                f = bcf_get_qual(v);
                value_exists = true;
            }
        }
        else if (type==VT_FILTER)
        {
            if (bcf_has_filter(h, v, tag.s)!=1)
            {
                b = false;
            }
            else
            {
                b = true;
            }
        }
        else if (type==VT_N_FILTER)
        {
            i = bcf_get_n_filter(v);
            f = i;
            b = true;
            value_exists = true;
        }
        else if (type==VT_INFO) //figure out type
        {
            if (!bcf_hdr_idinfo_exists(h, BCF_HL_INFO, bcf_hdr_id2int(h, BCF_DT_ID, tag.s)))
            {
                fprintf(stderr, "[%s:%d %s] INFO tag %s does not exist in header of VCF file.\n", __FILE__, __LINE__, __FUNCTION__, tag.s);
                exit(1);
            }

            int32_t info_type = bcf_hdr_id2type(h, BCF_HL_INFO, bcf_hdr_id2int(h, BCF_DT_ID, tag.s));

            if (info_type==BCF_HT_FLAG)
            {
                type |= VT_FLG;
                if (bcf_get_info_flag(h, v, tag.s, 0, 0)>0)
                {
                    i = 1;
                    f = 1;
                    b = true;
                    s.l=0;
                }
                else
                {
                    b = false;
                    value_exists = false;
                }
            }
            else if (info_type==BCF_HT_INT)
            {
                type |= VT_INT;
                int32_t ns = 0;
                int32_t *is = NULL;
                if (bcf_get_info_int32(h, v, tag.s, &is, &ns)>0)
                {
                    i = is[0];
                    f = (float)is[0];
                }
                else
                {
                    b = false;
                    value_exists = false;
                }

                if (ns) free(is);
            }
            else if (info_type==BCF_HT_REAL)
            {
                type |= VT_FLT;
                int32_t ns = 0;
                float *fs = NULL;
                if (bcf_get_info_float(h, v, tag.s, &fs, &ns)>0)
                {
                    b = true;
                    f = (float)fs[0];
                }
                else
                {
                    b = false;
                    value_exists = false;
                }

                if (ns) free(fs);
            }
            else if (info_type==BCF_HT_STR)
            {
                type |= VT_STR;
                //todo: how do you handle a vector of strings?

                int32_t n = s.m;

                if (bcf_get_info_string(h, v, tag.s, &s.s, &n)>0)
                {
                    b = true;
                    s.m = n;
                }
                else
                {
                    b = false;
                    value_exists = false;
                }
            }
        }
        else if (type==(VT_INFO|VT_INT))
        {
            int32_t *data = NULL;
            int32_t n=0;

            if (bcf_get_info_int32(h, v, tag.s, &data, &n)>0)
            {
                b = true;
                i = *((int*)data);
                f = (float) i;
            }
            else
            {
                b = false;
                value_exists = false;
            }

            if (n) free(data);
        }
        else if (type==(VT_INFO|VT_FLT))
        {
            int32_t *data = NULL;
            int32_t n=0;

            if (bcf_get_info_float(h, v, tag.s, &data, &n)>0)
            {
                b = true;
                f = *((float*)data);
            }
            else
            {
                b = false;
                value_exists = false;
            }

            if (n) free(data);
        }
        else if (type==(VT_INFO|VT_STR))
        {
            int32_t n=s.m;

            if (bcf_get_info_string(h, v, tag.s, &s.s, &n)>0)
            {
                b = true;
                s.m=n;
            }
            else
            {
                b = false;
                value_exists = false;
            }
        }
        else if (type==(VT_INFO|VT_FLG))
        {
            if (bcf_get_info_flag(h, v, tag.s, 0, 0)>0)
            {
                b = true;
            }
            else
            {
                b = false;
            }

            if (debug)
                std::cerr << "\tVT_INFO|VT_FLG "   << i << " " << f << " " << b << " " << s.s <<  " \n";
        }
        else if (type==VT_VARIANT_TYPE)
        {
            if (debug)
                std::cerr << "\tVTYPE "   <<  Variant::vtype2string(variant->type) <<  " \n";
            i = variant->type;
            b = i;
        }
        else if (type==VT_VARIANT_DLEN)
        {
            if (debug)
                std::cerr << "\tDLEN "   <<  variant->alleles[0].dlen <<  " \n";
            i = variant->alleles[0].dlen;
            f = i;
        }
        else if (type==VT_VARIANT_LEN)
        {
            if (debug)
                std::cerr << "\tLEN "   <<  abs(variant->alleles[0].dlen) <<  " \n";
            i = abs(variant->alleles[0].dlen);
            f = i;
        }
        else if (type==VT_N_ALLELE)
        {
            if (debug)
                std::cerr << "\tN_ALLELE "   <<  bcf_get_n_allele(v) <<  " \n";
            i = bcf_get_n_allele(v);
            f = i;
        }
    }
    else if (type&VT_MATH_OP)
    {
        if (!left->value_exists && !right->value_exists)
        {
            value_exists = false;
            return;
        }

        if ((type&15)==11) //ADD
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    i = (left->i+right->i);
                    f = i;
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
        else if ((type&15)==12) //SUB
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    i = (left->i-right->i);
                    f = i;
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
        else if ((type&15)==13) //MUL
        {
            if ((left->type&VT_INT))
            {
                if ((right->type&VT_INT))
                {
                    type |= VT_INT;
                    i = (left->i*right->i);
                    f = i;
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
        else if ((type&15)==14) //DIV
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
        else if ((type&15)==15) //BIT AND
        {
            if ((left->type&VT_INT) && (right->type&VT_INT))
            {
                i = (left->i & right->i);
                b = i;
                return;
            }

            fprintf(stderr, "[%s:%d %s] evaluation not supported for & :  %d %d\n", __FILE__, __LINE__, __FUNCTION__, left->type, right->type);
            exit(1);
        }
        else if ((type&15)==16) //BIT OR
        {
            if ((left->type&VT_INT) && (right->type&VT_INT))
            {
                i = (left->i | right->i);
                b = i;
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
 * Converts type to string.
 */
std::string Node::type2string(int32_t type)
{
    std::string s = "";
    if (type&VT_BOOL)
    {
        s += (s==""? "" : "|");
        s += "BOOL";
    }

    if (type&VT_INT)
    {
        s += (s==""? "" : "|");
        s += "INT";
    }

    if (type&VT_FLT)
    {
        s += (s==""? "" : "|");
        s += "FLT";
    }

    if (type&VT_STR)
    {
        s += (s==""? "" : "|");
        s += "STR";
    }

    if (type&VT_FLG)
    {
        s += (s==""? "" : "|");
        s += "FLG";
    }

    if (type&VT_LOGIC_OP)
    {
        s += (s==""? "" : "|");
        s += "LOGIC";
    }

    if (type&VT_MATH_CMP)
    {
        s += (s==""? "" : "|");
        s += "MATH_CMP";
    }

    if (type&VT_MATH_OP)
    {
        s += (s==""? "" : "|");
        s += "MATH_OP";
    }

    if (type&VT_BCF_OP)
    {
        s += (s==""? "" : "|");
        s += "BCF_OP";
    }

    if (type==VT_QUAL)
    {
        s += (s==""? "" : "|");
        s += "QUAL";
    }

    if (type==VT_FILTER)
    {
        s += (s==""? "" : "|");
        s += "FILTER";
    }

    if (type==VT_N_FILTER)
    {
        s += (s==""? "" : "|");
        s += "FILTER";
    }

    if (type==VT_INFO)
    {
        s += (s==""? "" : "|");
        s += "INFO";
    }

    if (type==VT_N_ALLELE)
    {
        s += (s==""? "" : "|");
        s += "N_ALLELE";
    }

    if (type==VT_VARIANT_TYPE)
    {
        s += (s==""? "" : "|");
        s += "VARIANT_TYPE";
    }

    if (type==VT_VARIANT_DLEN)
    {
        s += (s==""? "" : "|");
        s += "VARIANT_DLEN";
    }

    if (type==VT_VARIANT_LEN)
    {
        s += (s==""? "" : "|");
        s += "VARIANT_LEN";
    }

    return s;
};

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
    this->tree = NULL;
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

        if (!tree->type&VT_BOOL)
        {
            fprintf(stderr, "[%s:%d %s] filter expression not boolean %s\n", __FILE__, __LINE__, __FUNCTION__, exp);
            exit(1);
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

    if (tree->value_exists)
    {
        return tree->b;
    }
    else
    {
        return false;
    }
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
    trim_brackets(exp, len, debug);

    //this is a literal
    if (is_literal(exp, len, debug))
    {
        //will not recurse any further
        return parse_literal(exp, len, node, debug);
    }
    //unary operation
    else if (is_unary_op(exp, len, debug))
    {
        node->type = VT_NOT;
        if (debug) std::cerr << "\tis not_op\n";

        node->left = new Node();
        parse(exp+1, len-1, node->left, debug);
    }
    //binary operation
    else
    {
        if (debug) std::cerr << "\tis binary op\n";

        const char* p = exp; //points to end of first part
        const char* q = exp; //points to start of second part
        const char* r = exp; //for iteration

        int32_t type = INT_MAX;

        while(r-exp!=len)
        {
            //bypasses bracketed expressions
            fwd_to_closing_bracket(r, len);

            int32_t oplen = -1;
            int32_t ctype = peek_op(r, len, oplen, debug);

            if(ctype!=-1)
            {
                if (debug) std::cerr<< "\ttype : " << ctype << " \n";
                //this implements order of operations
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

        //valid binary operator not found
        if (type==INT_MAX)
        {
            kstring_t s = {0,0,0};
            kputsn(exp, len, &s);
            fprintf(stderr, "[%s:%d %s] binary operator not found in \"%s\". Valid operators are  ==,!=,=~,~~,&&,||,&,|,+,-,*,/\n", __FILE__, __LINE__, __FUNCTION__, s.s);
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
 * Checks if exp is a unary op.
 */
bool Filter::is_unary_op(const char* exp, int32_t len, bool debug)
{
    //check to make sure not a binary operator

    //NOT operator
    if (exp[0]=='~')
    {
        if (is_literal(exp+1, len-1, debug)||is_bracketed_expression(exp+1, len-1, debug))
        {
            if (debug) std::cerr << "\tis unary op\n";
            return true;
        }
    }

    if (debug) std::cerr << "\tis not unary op\n";
    return false;
}

/**
 * Checks is expression is bracketed.
 */
bool Filter::is_bracketed_expression(const char* exp, int32_t len, bool debug)
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

            ++j;
        }

        return opened_brackets==1;
    }

    return false;
}

/**
 * Check if exp is a literal (no binary operations, no unary operation).
 */
bool Filter::is_literal(const char* exp, int32_t len, bool debug)
{
    //checks if string literal
    if (exp==strchr(exp,'\'') && exp+len-1==strchr(exp+1,'\''))
    {
        return true;
    }

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
           *q=='~' ||
           (*q=='-' && exp!=q))
        {
            if (debug) std::cerr << "\tis not literal\n";
            return false;
        }

        ++q;
    }

    if (debug) std::cerr << "\tis literal\n";
    return true;
}

/**
 * Parse literals.
 */
void Filter::parse_literal(const char* exp, int32_t len, Node * node, bool debug)
{
    node->type = VT_UNKNOWN;

    if (strncmp(exp, "PASS", len)==0)
    {
        node->type = VT_FILTER;
        kputsn(exp, len, &node->tag);
        if (debug) std::cerr << "\tis filter_op\n";
        return;
    }
    else if (strncmp(exp, "QUAL", 4)==0)
    {
        node->type = VT_QUAL;
        exp += 4;
        if (debug) std::cerr << "\tis qual_op\n";
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
    else if (strncmp(exp, "N_FILTER", len)==0)
    {
        node->type = VT_N_FILTER;
        if (debug) std::cerr << "\tis n_filter_op\n";
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
        node->type = VT_INT;
        node->i = VT_INDEL;
        node->value_exists = true;
        if (debug) std::cerr << "\tis INDEL\n";
        return;
    }
    else if (strncmp(exp, "SNP", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_SNP;
        node->value_exists = true;
        if (debug) std::cerr << "\tis SNP\n";
        return;
    }
    else if (strncmp(exp, "MNP", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_MNP;
        node->value_exists = true;
        if (debug) std::cerr << "\tis MNP\n";
        return;
    }
    else if (strncmp(exp, "CLUMPED", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_CLUMPED;
        node->value_exists = true;
        if (debug) std::cerr << "\tis CLUMPED\n";
        return;
    }
    else if (strncmp(exp, "VNTR", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_VNTR;
        if (debug) std::cerr << "\tis VNTR\n";
        return;
    }
    else if (strncmp(exp, "SV", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_SV;
        if (debug) std::cerr << "\tis SV\n";
        return;
    }
    else if (strncmp(exp, "REF", len)==0)
    {
        node->type = VT_INT;
        node->i = VT_REF;
        node->value_exists = true;
        if (debug) std::cerr << "\tis REF\n";
        return;
    }
    else if (strncmp(exp, "DLEN", len)==0)
    {
        node->type = VT_VARIANT_DLEN;
        node->value_exists = false;
        if (debug) std::cerr << "\tis dlen\n";
        return;
    }
    else if (strncmp(exp, "LEN", len)==0)
    {
        node->type = VT_VARIANT_LEN;
        node->value_exists = false;
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
            node->value_exists = true;
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
            node->value_exists = true;
            if (debug) std::cerr << "\tis float: " << f << "\n";
            return;
        }

        //string type
        if (exp[0]=='\'' && exp[len-1]=='\'')
        {
            node->type = VT_STR;
            kputsn(exp+1, len-2, &node->tag);
            kputc(0, &node->tag);
            kputsn(exp+1, len-2, &node->s);
            kputc(0, &node->s);
            node->value_exists = true;
            if (debug) std::cerr << "\tis string\n";
            return;
        }

        if (node->type==VT_UNKNOWN)
        {
            kstring_t s = {0,0,0};
            kputsn(exp, len, &s);
            fprintf(stderr, "[E:%s] %s is not recognized. If it is a string, you should place it in single inverted commas.\n", __FUNCTION__, s.s);
            print_filter_help();
            exit(1);
        }
    }

    if (debug)
    {
        std::cerr << "\tvalue_exists " << node->value_exists << "\n";
        std::cerr << "\ttag          " << node->tag.s << "\n";
        std::cerr << "\tb            " << node->b << "\n";
        std::cerr << "\ti            " << node->i << "\n";
        std::cerr << "\tf            " << node->f << "\n";
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
    fprintf(stderr, "  QUAL field\n");
    fprintf(stderr, "    QUAL\n");
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
    fprintf(stderr, "  Regular expressions for string fields using pcre2\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Passed variants in intergenic regions or UTR : PASS&&INFO.ANNO=~'Intergenic|UTR'\n");
    fprintf(stderr, "  Passed variants in intergenic regions or UTR : PASS&&INFO.ANNO=~'(?i)Intergenic|UTR'\n");
    fprintf(stderr, "  ignoring case\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Operations\n");
    fprintf(stderr, "     ==,!=,=~,~~,~,&&,||,&,|,+,-,*,/\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Failed rare variants : ~PASS&&(INFO.AC/INFO.AN<0.005)\n");
    fprintf(stderr, "\n");
}

/**
 * Trim brackets from an expression.
 */
void Filter::trim_brackets(const char* &exp, int32_t &len, bool debug)
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
                if (debug) std::cerr << "\t\ttrimmed brackets\n";
                trim_brackets(exp, len, debug);
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
    else if (*s=='=' && (s+1-r<len) && *(s+1)=='~')
    {
        if (debug) std::cerr << "\tis =~ operator\n";
        oplen = 2;
        return VT_MATCH;
    }
    else if (*s=='~' && (s+1-r<len) && *(s+1)=='~')
    {
        if (debug) std::cerr << "\tis ~~ operator\n";
        oplen = 2;
        return VT_NO_MATCH;
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
    node->value_exists = false;

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
            if (node->left->value_exists && !node->left->b)
            {
                node->b = false;
                node->value_exists = true;
                return;
            }
            else
            {
                 apply(node->right, debug);
            }
        }
        else if (node->type==VT_OR)
        {
            if (node->left->value_exists && node->left->b)
            {
                node->b = true;
                node->value_exists = true;
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
};