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

#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <sstream>
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/tbx.h"
#include "rb_tree.h"
#include "hts_utils.h"

class IntervalTreeNode
{
    public:
    std::vector<Interval*> intervals;
    IntervalTreeNode* parent;
    IntervalTreeNode* left;
    IntervalTreeNode* right;
    int32_t start;
    int32_t max;
    int32_t min;
    uint32_t color;
        
    IntervalTreeNode(Interval* item);
    ~IntervalTreeNode();
    
    void print();
    
    void insert(Interval* item);
    
    private:
};
      
class IntervalTree
{
    public:
    IntervalTreeNode* root;   
    uint32_t noElements;
    uint32_t height;
    
    IntervalTree();
    ~IntervalTree();
    IntervalTreeNode* simple_insert(Interval* x); 
    void insert(Interval* x); 
    void search(int32_t low, int32_t high, std::vector<Interval*>& intervals);
    void search_iter(int32_t low, int32_t high, std::vector<Interval*>& intervals, IntervalTreeNode* x);
    void searchBrute(int32_t low, int32_t high, std::vector<Interval*>& intervals);
    void search_iter_brute(int32_t low, int32_t high, std::vector<Interval*>& intervals, IntervalTreeNode* x);
    void print();
    void print_iter(IntervalTreeNode* x);
    void validate();
    void validate_iter(IntervalTreeNode* x, uint32_t height);
    uint32_t size()
    {
        return noElements;
    };
                    
    private:
    void leftRotate(IntervalTreeNode* x);    
    void rightRotate(IntervalTreeNode* y);        
};   
    
#endif