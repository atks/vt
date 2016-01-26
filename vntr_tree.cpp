/* The MIT License
   Copyright (c) 2015 Adrian Tan <atks@umich.edu>
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

#include "vntr_tree.h"

/**
 * Constructor.
 */
VNTRNode::VNTRNode()
{

};

/**
 * Constructor.
 */
VNTRNode::VNTRNode(std::string motif, std::string basis, int32_t exact_count, int32_t fuzzy_count)
{
    this->motif = motif;
    this->basis = basis;
    this->exact_count = exact_count;
    this->fuzzy_count = fuzzy_count;
};

/**
 * Destructor.
 */
VNTRNode::~VNTRNode()
{
};

/**
 * Clear values.
 */
void VNTRNode::clear()
{
};

/**
 * Print values.
 */
void VNTRNode::print()
{
};

/**
 * Constructor.
 */
VNTRTree::VNTRTree()
{
    vntrs[0].resize(1);
    vntrs[1].resize(1);
    vntrs[2].resize(1);
    vntrs[3].resize(1);
    
    //for all subsets of ACGT
    std::string motifs[15] = {"A", "C", "G", "T",
                              "AC", "AG", "AT", "CG", "CT", "GT",
                              "ACG", "ACT", "AGT", "CGT",
                              "ACGT"};
    
    for (uint32_t i=0; i<15; ++i)
    {
        VNTRNode* node = new VNTRNode(motifs[i], motifs[i], 0, 0);
        motif_map[motifs[i]] = node;
        int32_t basis_len = motifs[i].size();
        int32_t motif_len = motifs[i].size();
        vntrs[basis_len-1][motif_len-basis_len].push_back(node);
    }
};

/**
 * Destructor.
 */
VNTRTree::~VNTRTree()
{
};

/**
 * Assumes that the variant is appropriately updated on its VNTR chracteristics.
 * Updates tree on motif, exact/fuzziness of VNTR.
 */
void VNTRTree::count(Variant& variant)
{
    if (variant.type == VT_VNTR)
    {
        VNTR& vntr = variant.vntr;

        VNTRNode* node = NULL; 
        if (motif_map.find(vntr.motif)==motif_map.end())
        {
            node = new VNTRNode(vntr.motif, vntr.basis, 0, 0);
            motif_map[vntr.motif] = node;
            int32_t basis_len = vntr.basis.size();
            int32_t motif_len = vntr.motif.size();
            
            if (vntrs[basis_len-1].size()<motif_len-basis_len+1)
            {        
                vntrs[basis_len-1].resize(motif_len-basis_len+1);
            }
            vntrs[basis_len-1][motif_len-basis_len].push_back(node);
            
        }
        else
        {
            node = motif_map[vntr.motif];
        }
        
        
        float concordance = -1;
        float *score = NULL;
        int32_t n = 0;
        if (bcf_get_info_float(variant.h, variant.v, "SCORE", &score, &n)>0)
        {
            concordance = score[0];
            free(score);
        }
        else if (bcf_get_info_float(variant.h, variant.v, "CONCORDANCE", &score, &n)>0)
        {
            concordance = score[0];
            free(score);
        }
        
        if (concordance==1 ||
           (vntr.fuzzy_rbeg1!=0 && vntr.exact_rbeg1==vntr.fuzzy_rbeg1 &&
            vntr.exact_rend1==vntr.fuzzy_rend1))
        {
            ++node->exact_count;
        }
        else
        {
            ++node->fuzzy_count;
        }
    }
};

/**
 * Print this tree.
 */
void VNTRTree::print(int32_t level)
{   
    if (level==SEQUENCE)
    {
        for (uint32_t basis_len=1; basis_len<=4; ++basis_len)
        {
            std::cerr << "basis length " << basis_len << "\n";
            for (uint32_t motif_len=basis_len; motif_len-basis_len<vntrs[basis_len-1].size(); ++motif_len)
            {
                std::cerr << "\tmotif length " << motif_len << "\n";
                std::list<VNTRNode*>::iterator i;
                uint32_t exact_count = 0;
                uint32_t fuzzy_count = 0;
                for (i=vntrs[basis_len-1][motif_len-basis_len].begin(); i!=vntrs[basis_len-1][motif_len-basis_len].end(); ++i)
                {
                    VNTRNode& vntr_node = **i;
                    std::cerr << "\t\t" << vntr_node.motif << " " << (vntr_node.exact_count+vntr_node.fuzzy_count) << " (" << vntr_node.exact_count << "/" << vntr_node.fuzzy_count << ")\n";
                    exact_count += vntr_node.exact_count;
                    fuzzy_count += vntr_node.fuzzy_count;
                }
                std::cerr << "\t\ttotal\t" << (exact_count+fuzzy_count) << " (" << exact_count << "/" << fuzzy_count << ")\n";
            }
        }
    }   
    else if (level==MOTIF)
    {
        for (uint32_t basis_len=1; basis_len<=4; ++basis_len)
        {
            std::cerr << "basis length " << basis_len << "\n";
            std::cerr << "\tmotif length\n";
            for (uint32_t motif_len=basis_len; motif_len-basis_len<vntrs[basis_len-1].size(); ++motif_len)
            {
                std::list<VNTRNode*>::iterator i;
                uint32_t exact_count = 0;
                uint32_t fuzzy_count = 0;
                for (i=vntrs[basis_len-1][motif_len-basis_len].begin(); i!=vntrs[basis_len-1][motif_len-basis_len].end(); ++i)
                {
                    VNTRNode& vntr_node = **i;
                    exact_count += vntr_node.exact_count;
                    fuzzy_count += vntr_node.fuzzy_count;
                }
                std::cerr << "\t\t" << motif_len << "\t" << (exact_count+fuzzy_count) << " (" << exact_count << "/" << fuzzy_count << ")\n";
                
            }
        }
    } 
    else if (level==BASIS)
    {
        for (uint32_t basis_len=1; basis_len<=4; ++basis_len)
        {
            uint32_t basis_exact_count = 0;
            uint32_t basis_fuzzy_count = 0;
            for (uint32_t motif_len=basis_len; motif_len-basis_len<vntrs[basis_len-1].size(); ++motif_len)
            {
                std::list<VNTRNode*>::iterator i;
                uint32_t exact_count = 0;
                uint32_t fuzzy_count = 0;
                for (i=vntrs[basis_len-1][motif_len-basis_len].begin(); i!=vntrs[basis_len-1][motif_len-basis_len].end(); ++i)
                {
                    VNTRNode& vntr_node = **i;
                    exact_count += vntr_node.exact_count;
                    fuzzy_count += vntr_node.fuzzy_count;
                }
                basis_exact_count += exact_count;
                basis_fuzzy_count += fuzzy_count;                
            }
            
            std::cerr << "basis length " << basis_len << "\t" << (basis_exact_count+basis_fuzzy_count) << " (" << basis_exact_count << "/" << basis_fuzzy_count << ")\n";
        }
    } 
};