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

#include "utils.h"
#include "interval.h"

#define RED 1
#define BLACK 2

class IntervalTreeNode
{
    public:
    IntervalTreeNode* parent;
    IntervalTreeNode* left;
    IntervalTreeNode* right;
    int32_t start;
    int32_t max;
    int32_t min;
    std::vector<Interval*> intervals;
    uint32_t color;

    /**
     * Constructs an IntervalTreeNode and initialize it with an interval.
     */
    IntervalTreeNode(Interval* interval);

    /**
     * Constructs an IntervalTreeNode.
     */
    ~IntervalTreeNode();

    /**
     * Insert an interval.
     */
    void insert(Interval* interval);

    /**
     * Prints the node.
     */
    void print();

    private:
};

class IntervalTree
{
    public:
    IntervalTreeNode* root;
    uint32_t no_elements;
    uint32_t height;

    /**
     * Constructor.
     */
    IntervalTree();

    /**
     * Destructor.
     */
    ~IntervalTree();

    /**
     * Returns the number of intervals in the tree.
     */
    uint32_t size();

    /**
     * Insert an interval, returns a node if the insertion causes violation of the red black tree.
     */
    void insert(Interval* interval);

    /**
     * Gets overlapping intervals with [start,end]. Returns true on success and false otherwise.
     */
    void search(int32_t start, int32_t end, std::vector<Interval*>& intervals);

    /**
     * Brute force recursive search for overlap for sanity checks.
     */
    void search_brute(int32_t start, int32_t end, std::vector<Interval*>& intervals);

    /**
     * Prints the tree.
     */
    void print();

    /**
     * Validates red black tree property.
     */
    void validate();

    private:

    /**
     * Insert an interval, returns a node if the insertion causes violation of the red black tree.
     */
    IntervalTreeNode* simple_insert(Interval* interval);

    /**
     * Iterative method for search_brute.
     */
    void search_iter_brute(int32_t start, int32_t end, std::vector<Interval*>& intervals, IntervalTreeNode* x);

    /**
     * Iterative method for search.
     */
    void search_iter(int32_t start, int32_t end, std::vector<Interval*>& intervals, IntervalTreeNode* x);

    /**
     * Iterative method for print.
     */
    void print_iter(IntervalTreeNode* x);

    /**
     * Iterative method for validate.
     */
    void validate_iter(IntervalTreeNode* x, uint32_t depth);

    /**
     * Left rotates a node.
     */
    void left_rotate(IntervalTreeNode* x);

    /**
     * Righ rotates a node.
     */
    void right_rotate(IntervalTreeNode* y);
};

#endif