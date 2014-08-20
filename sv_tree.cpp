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
    this->desc = {0,0,0};
    this->depth = 0;
    this->count = 0;
};

/**
 * Constructor.
 */
SVNode::SVNode(const char* desc)
{
    this->parent = NULL;
    this->desc = {0,0,0};
    this->depth = 0;
    count = 0;
    kputs(desc, &this->desc);
};

/**
 * Constructor.
 */
SVNode::SVNode(const char* desc, int32_t depth)
{
    this->parent = NULL;
    this->desc = {0,0,0};
    this->depth = depth;
    count = 0;
    kputs(desc, &this->desc);
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

    if (desc.m) free(desc.s);
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
    ++count;
};

/**
 * Clear values.
 */
void SVNode::clear()
{
    if (desc.m) free(desc.s);
    count = 0;
};

/**
 * Print values.
 */
void SVNode::print()
{
    for (int32_t i=0; i<depth; ++i)
        std::cerr << "\t";
    std::cerr << tags2desc() << " (" << count << ")\n";

    for (int32_t i=0; i<children.size(); ++i)
    {
        children[i]->print();
    }
};

/**
 *  For translating reserved keywords.
 */
const char* SVNode::tags2desc()
{
    if (!strcmp(desc.s,"TRA"))
    {
        return "translocation";
    }
    else if (!strcmp(desc.s,"DEL"))
    {
        return "deletion";
    }
    else if (!strcmp(desc.s,"INS"))
    {
        return "insertion";
    }
    else if (!strcmp(desc.s,"INV"))
    {
        return "inversion";
    }
    else if (!strcmp(desc.s,"ME"))
    {
        return "mobile element";
    }
    else if (!strcmp(desc.s,"MT"))
    {
        return "nuclear mitochondrial";
    }
    else if (!strcmp(desc.s,"DUP"))
    {
        return "duplication";
    }
    else if (!strcmp(desc.s,"TANDEM"))
    {
        return "tandem repeats";
    }
    else if (!strcmp(desc.s,"CNV"))
    {
        return "copy number variation";
    }
    else
    {
        return desc.s;
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
bool SVTree::add(const char* desc)
{
    //update hash
    khiter_t k;
    int32_t ret = 0;
    if ((k=kh_get(xdict, m, desc))==kh_end(m))
    {
        k = kh_put(xdict, m, desc, &ret);
        SVNode* svnode = NULL;
        if (ret)
        {
            svnode = new SVNode(desc);
            kh_value(m, k) = svnode;
        }
        else
        {
            kh_value(m, k)->clear();
        }

        //update tree
        std::vector<std::string> vec;
        split(vec, "<:>", desc);

        SVNode* cnode = root;
        for (int32_t i=0; i<vec.size(); ++i)
        {
            bool found_type = false;
            for (int32_t j=0; j<cnode->children.size(); ++j)
            {
                if (!strcmp(cnode->children[j]->desc.s, vec[i].c_str()))
                {
                    cnode = cnode->children[j];
                    found_type = true;
                    break;
                }
            }

            if (!found_type)
            {
                max_depth = i+1>max_depth?i+1:max_depth;
                SVNode* newnode = new SVNode(vec[i].c_str(), i+1);
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
void SVTree::count(char* desc)
{
    khiter_t k;
    int32_t ret = 0;
    if ((k=kh_get(xdict, m, desc))==kh_end(m))
    {
        this->add(desc);
        k=kh_get(xdict, m, desc);
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
