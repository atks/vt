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

#ifndef RB_TREE_H
#define RB_TREE_H

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <list>

#define RED 1
#define BLACK 2

class Interval
{
    public:
    int32_t start;
    int32_t end;

};

class RBTreeNode
{
    public:
    std::vector<Interval*> intervals;
    RBTreeNode* parent;
    RBTreeNode* left;
    RBTreeNode* right;
    int32_t start;
    int32_t max;
    int32_t min;
    uint32_t color;

    RBTreeNode(Interval* item);

    void print();

    void insert(Interval* item);

    private:
};

class RBTree
{
    public:
    RBTreeNode* root;
    uint32_t noElements;
    uint32_t height;

    RBTree();
    ~RBTree();
    RBTreeNode* insert(Interval* x);
    void RBinsert(Interval* x);
    void RBSearch(int32_t low, int32_t high, std::vector<Interval*>& intervals);
    void RBSearch_iter(int32_t low, int32_t high, std::vector<Interval*>& intervals, RBTreeNode* x);
    void RBSearchBrute(int32_t low, int32_t high, std::vector<Interval*>& intervals);
    void RBSearch_iter_brute(int32_t low, int32_t high, std::vector<Interval*>& intervals, RBTreeNode* x);
    void print();
    void print_iter(RBTreeNode* x);
    void validate();
    void validate_iter(RBTreeNode* x, uint32_t height);
    uint32_t size()
    {
        return noElements;
    };

    private:
    void leftRotate(RBTreeNode* x);
    void rightRotate(RBTreeNode* y);
};

#endif