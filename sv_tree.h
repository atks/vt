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

#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "variant_manip.h"

#ifndef SV_Tree_H
#define SV_Tree_H

#define VT_SV_DEL    1
#define VT_SV_INS    2
#define VT_SV_DUP    3
#define VT_SV_INV    4
#define VT_SV_CNV    5
#define VT_SV_TRA    6
#define VT_SV_TANDEM 7
#define VT_SV_ME     8
#define VT_SV_MT     9

/**
 * Class for filtering VCF records.
 */
class SVNode
{
    public:

    SVNode* parent;
    std::vector<SVNode*> children;

    int32_t depth, count;
    kstring_t desc;

    /**
     * Constructor.
     */
    SVNode()
    {
        desc = {0,0,0};
        depth = 0;
        count = 0;
    };

    /**
     * Constructor.
     */
    SVNode(const char* desc)
    {
        this->desc = {0,0,0};
        this->depth = 0;
        count = 0;
        kputs(desc, &this->desc);
    };

    /**
     * Constructor.
     */
    SVNode(const char* desc, int32_t depth)
    {
        this->desc = {0,0,0};
        this->depth = depth;
        count = 0;
        kputs(desc, &this->desc);
    };

    /**
     * Destructor.
     */
    ~SVNode()
    {
        for (int32_t i=0; i<children.size(); ++i)
        {
            delete children[i];
        }
        children.clear();

        if (desc.m) free(desc.s);
    };

    /**
     * Set depth.
     */
    void set_depth(int32_t depth)
    {
        this->depth = depth;
    };

    /**
     * Clear values.
     */
    void clear()
    {
        if (desc.m) free(desc.s);
        count = 0;
    };
};

KHASH_MAP_INIT_STR(xdict, SVNode*);

/**
 * SV Tree for counting SVs.
 */
class SVTree
{
    public:

    //filter expression
    SVNode* tree;

    //useful pointers for applying the filter to a vcf record
    bcf_hdr_t *h;
    bcf1_t *v;

    khash_t(xdict) *m;

    /**
     * Constructor.
     */
    SVTree()
    {
        tree = new SVNode("root", 0);
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
    };

    /**
     * Destructor.
     */
    ~SVTree()
    {
        if (tree)
        {
            delete tree;
        }
        tree = NULL;

        m = kh_init(xdict);
    };

    /**
     * Adds a new tag, returns true if successful.
     */
    bool add(const char* desc)
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

            SVNode* cnode = tree;
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
                    SVNode* newnode = new SVNode(vec[i].c_str(), i+1);
                    cnode->children.push_back(newnode);
                    cnode = newnode;
                }    
            }

            return true;
        }

        return false; // already exists
    }

    /**
     * Observes and update the count of a new tag.
     */
    void count(char* desc)
    {
         khiter_t k;
        int32_t ret = 0;
        if ((k=kh_get(xdict, m, desc))!=kh_end(m))
        {


        }
        else
        {
            this->add(desc);
            if (ret)
            {
                kh_value(m, k) = new SVNode();
            }
            else
            {
                kh_value(m, k)->clear();
            }

            kh_value(m, k)->count = 1;
        }
    };

    /**
     * Iterator, returns first node by depth first search.
     */
    SVNode* begin()
    {
        return NULL;
    };

    /**
     * Iterator, returns node by depth first search.
     */
    SVNode* next()
    {
        return NULL;
    };

    private:
};

#endif