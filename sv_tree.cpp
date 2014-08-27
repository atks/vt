/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#include "sv_tree.h"

/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#include "sv_tree.h"

/**
 * Constructor.
 */
SVNode::SVNode()
{
    this->parent = NULL;
    this->type = {0,0,0};
    this->depth = 0;
    this->count = 0;
};

/**
 * Constructor.
 */
SVNode::SVNode(const char* type)
{
    this->parent = NULL;
    this->type = {0,0,0};
    this->depth = 0;
    count = 0;
    kputs(type, &this->type);
};

/**
 * Constructor.
 */
SVNode::SVNode(const char* type, int32_t depth)
{
    this->parent = NULL;
    this->type = {0,0,0};
    this->depth = depth;
    count = 0;
    kputs(type, &this->type);
};

/**
 * Destructor.
 */
SVNode::~SVNode()
{
    for (int32_t i=0; i<children.size(); ++i)
    {
        delete children[i];
    }
    children.clear();
    this->parent = NULL;

    if (type.m) free(type.s);
};

/**
 * Set depth.
 */
void SVNode::set_depth(int32_t depth)
{
    this->depth = depth;
};

/**
 * Increment count.
 */
void SVNode::increment_count()
{
    std::cerr << "\tincrement " << this->type.s << "\n";
    
    ++count;

    if (parent!=NULL)
    {
        parent->increment_count();
    }
};

/**
 * Clear values.
 */
void SVNode::clear()
{
    if (type.m) free(type.s);
    count = 0;
};

/**
 * Print values.
 */
void SVNode::print()
{
    for (int32_t i=0; i<depth; ++i) std::cerr << "\t";
        
    std::cerr << sv_type2string(type.s) << " (" << count << ")\n";

    for (int32_t i=0; i<children.size(); ++i)
    {
        children[i]->print();
    }
};

/**
 *  For translating reserved keywords.
 */
std::string SVNode::sv_type2string(char* sv_type)
{
    if (!strcmp(sv_type,"TRA"))
    {
        return "translocation";
    }
    else if (!strcmp(sv_type,"DEL"))
    {
        return "deletion";
    }
    else if (!strcmp(sv_type,"INS"))
    {
        return "insertion";
    }
    else if (!strcmp(sv_type,"INV"))
    {
        return "inversion";
    }
    else if (!strcmp(sv_type,"ME"))
    {
        return "mobile element";
    }
    else if (!strcmp(sv_type,"MT"))
    {
        return "nuclear mitochondrial DNA";
    }
    else if (!strcmp(sv_type,"DUP"))
    {
        return "duplication";
    }
    else if (!strcmp(sv_type,"TANDEM"))
    {
        return "tandem repeats";
    }
    else if (!strcmp(sv_type,"CNV"))
    {
        return "copy number variation";
    }
    else
    {
        return sv_type;
    }
};

/**
 * Constructor.
 */
SVTree::SVTree()
{
    root = new SVNode("root", 0);
    m = kh_init(xdict);
    
    this->add("<TRA>");
    this->add("<DEL>");
    this->add("<INS>");
    this->add("<DUP>");
    this->add("<INV>");
    this->add("<CNV>");
    this->add("<DUP:TANDEM>");
    this->add("<DEL:ME>");
    this->add("<INS:ME>");
    this->add("<INS:MT>");
    this->add("<INS:ME:ALU>");
    this->add("<INS:ME:LINE1>");
    this->add("<INS:ME:SVA>");

    print();
};

/**
 * Destructor.
 */
SVTree::~SVTree()
{
    if (root)
    {
        delete root;
    }
    root = NULL;

    m = kh_init(xdict);
};

/**
 * Adds a new tag, returns true if successful.
 */
bool SVTree::add(const char* sv_type)
{    
    //update hash
    khiter_t k;
    int32_t ret = 0;
    if ((k=kh_get(xdict, m, sv_type))==kh_end(m))
    {
        k = kh_put(xdict, m, sv_type, &ret);
        std::vector<std::string> vec;
        split(vec, "<:>", sv_type);
        
        SVNode* cnode = root;
        for (size_t i=0; i<vec.size(); ++i)
        {
            bool found_type = false;
            for (size_t j=0; j<cnode->children.size(); ++j)
            {
                if (!strcmp(cnode->children[j]->type.s, vec[i].c_str()))
                {
                    cnode = cnode->children[j];
                    found_type = true;
                    break;
                }
            }

            if (!found_type)
            {
                max_depth = i+1>max_depth?i+1:max_depth;
                std::cerr << "adding " << vec[i] << "\n";
                SVNode* newnode = new SVNode(vec[i].c_str(), i+1);
                kh_value(m, k) = newnode;
                cnode->children.push_back(newnode);
                newnode->parent = cnode;
                cnode = newnode;
            }
        }

        return true;
    }

    return false;
}

/**
 * Observes and update the count of a new tag.
 */
void SVTree::count(Variant& variant)
{
    khiter_t k;
    int32_t ret = 0;

    for (size_t i=0; i<variant.alleles.size(); ++i)
    {
        variant.alleles[i].print();
        const char* sv_type = variant.alleles[i].sv_type.c_str();
        std::cerr  << "counting " << sv_type << "\n";
        if ((k=kh_get(xdict, m, sv_type))==kh_end(m))
        {
           std::cerr  << "\tadding for the first time\n";
           this->add(sv_type);
           k=kh_get(xdict, m, sv_type);
        }
    }
    
    kh_value(m, k)->increment_count();
};

/**
 * Enumerates the children in a depth first order.
 */
SVNode* SVTree::enumerate()
{
    return NULL;
};

/**
 * Iterator, returns first node by depth first search.
 */
SVNode* SVTree::begin()
{
    return NULL;
};

/**
 * Iterator, returns node by depth first search.
 */
SVNode* SVTree::next()
{
    return NULL;
};

/**
 * Print this tree.
 */
void SVTree::print()
{
    root->print();
};
