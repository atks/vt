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

#include "interval_tree.h"

/**
 * Constructs an IntervalTreeNode and initialize it with an interval.
 */
IntervalTreeNode::IntervalTreeNode(Interval* interval)
{
    insert(interval);
    color = BLACK;
    parent = NULL;
    left = NULL;
    right = NULL;
};

/**
 * Constructs an IntervalTreeNode.
 */
IntervalTreeNode::~IntervalTreeNode()
{
    parent = NULL;

    if (left!=NULL)
        delete left;

    if (right!=NULL)
        delete right;

    left = NULL;
    right = NULL;
};

/**
 * Insert an interval.
 */
void IntervalTreeNode::insert(Interval* interval)
{
    intervals.push_back(interval);
    max = std::max(interval->end, max);
    min = std::min(interval->start, min);
};

/**
 * Prints the node.
 */
void IntervalTreeNode::print()
{
    std::cerr << "start   : " << start << "\n";
    std::cerr << "address : " << this << "\n";
    std::cerr << "max     : " << max << "\n";
    std::cerr << "min     : " << min << "\n";
    std::cerr << "color   : " << color << "\n";
    std::cerr << "child   : (" << left << "," << right << ")\n";
};

/**
 * Constructor.
 */
IntervalTree::IntervalTree()
{
    root = NULL;

    height = 0;
    no_elements = 0;
};

/**
 * Destructor.
 */
IntervalTree::~IntervalTree()
{
    delete root;
    root = NULL;
};

/**
 * Returns the number of intervals in the tree.
 */
uint32_t IntervalTree::size()
{
    return no_elements;
};

/**
 * Insert an interval, returns a node if the insertion causes violation of the red black tree.
 */
void IntervalTree::insert(Interval* interval)
{
    IntervalTreeNode* x = simple_insert(interval);

    //if no new node inserted
    if (x==NULL)
    {
        //no need to correct tree
        return;
    }

    x->color = RED;

    //percolate up
    while (x!=root && x->parent->color==RED) //violation of property 3
    {
        //x->parent->max = std::max(x->max,x->parent->max);
       // x->parent->parent->max = std::max(x->parent->parent->max,x->parent->max);

        //if parent is a left child
        if (x->parent==x->parent->parent->left)
        {
            IntervalTreeNode* y = x->parent->parent->right;
            //parent and avuncular have same height, no need to rotate
            //color grandparent red and "recurse"
            if (y!=NULL && y->color==RED)
            {
                x->parent->color = BLACK;
                y->color = BLACK;
                x->parent->parent->color = RED; //move up to grandparent
                x = x->parent->parent;
            }
            //parent and avuncular have differing height, rotate
            else
            {
                //shift weight
                if (x==x->parent->right)
                {
                    x = x->parent;
                    left_rotate(x);
                }
                else
                {
                    x->parent->color = BLACK;
                    x->parent->parent->color = RED;
                    right_rotate(x->parent->parent);
                }
            }
        }
        else //parent is a right child
        {
            IntervalTreeNode* y = x->parent->parent->left;
            if (y!=NULL && y->color==RED)
            {
                x->parent->color = BLACK;
                y->color = BLACK;
                x->parent->parent->color = RED;
                x = x->parent->parent;
            }
            else
            {
                if (x==x->parent->left)
                {
                    x = x->parent;
                    right_rotate(x);
                }
                else
                {
                    x->parent->color = BLACK;
                    x->parent->parent->color = RED;
                    left_rotate(x->parent->parent);
                }
            }
        }
    }

    root->color = BLACK;
};

/**
 * Gets overlapping intervals with [start,end].
 */
void IntervalTree::search(int32_t start, int32_t end, std::vector<Interval*>& intervals)
{
    intervals.clear();

    if (root==NULL)
        return;

    search_iter(start, end, intervals, root);
};

/**
 * Brute force recursive search for overlap for sanity checks.
 */
void IntervalTree::search_brute(int32_t start, int32_t end, std::vector<Interval*>& intervals)
{
    intervals.clear();

    if (root==NULL)
        return;

    search_iter_brute(start, end, intervals, root);
};

/**
 * Prints the tree.
 */
void IntervalTree::print()
{
    print_iter(root);
    std::cerr << "\n";
};

/**
 * Validates red black tree property.
 */
void IntervalTree::validate()
{
    height = 0;
    validate_iter(root, 0);
};

/**
 * Insert an interval, returns a node if the insertion causes violation of the red black tree.
 */
IntervalTreeNode* IntervalTree::simple_insert(Interval* interval)
{
    ++no_elements;

    IntervalTreeNode* y = NULL;
    IntervalTreeNode* x = root;

    int32_t start = interval->start;

    //search for leaf node to append z
    while (x!=NULL)
    {
        y = x;

        x->max = std::max(x->max, interval->end);
        x->min = std::min(x->min, interval->start);

        if (start < x->start)
        {
            x = x->left;
        }
        else if (start > x->start)
        {
            x = x->right;
        }
        else //equal
        {
            break;
        }
    }

    //identical key found
    if (y!=NULL && start == y->start)
    {
        y->insert(interval);

        return NULL;
    }
    else
    {
        IntervalTreeNode* z = new IntervalTreeNode(interval);
        z->start = interval->start;
        z->max = interval->end;
        z->min = interval->start;

        z->parent = y;

        //if first element in tree
        if (y == NULL)
        {
            root = z;
            z->color = BLACK;
        }
        //attach to a parent
        else
        {
            if (z->start < y->start)
            {
                y->left = z;
            }
            else
            {
                y->right = z;
            }
        }

        return z;
    }
};

/**
 * Iterative method for search_brute.
 */
void IntervalTree::search_iter_brute(int32_t start, int32_t end, std::vector<Interval*>& intervals, IntervalTreeNode* x)
{
    //check through list
    if (x->start <= end)
    {
        for (uint32_t i=0; i<x->intervals.size(); ++i)
        {
            //overlap
            if (x->intervals[i]->end >= start)
            {
                intervals.push_back(x->intervals[i]);
            }
        }
    }

    if (x->left!=NULL)
    {
        search_iter_brute(start, end, intervals, x->left);
    }

    if (x->right!=NULL)
    {
        search_iter_brute(start, end, intervals, x->right);
    }
};

/**
 * Iterative method for search.
 */
void IntervalTree::search_iter(int32_t start, int32_t end, std::vector<Interval*>& intervals, IntervalTreeNode* x)
{
    //check through list
    if (x->start <= end)
    {
        for (uint32_t i=0; i<x->intervals.size(); ++i)
        {
            //overlap
            if (x->intervals[i]->end >= start)
            {
                intervals.push_back(x->intervals[i]);
            }
        }
    }

    if (x->left!=NULL && x->left->max >= start)
    {
        search_iter(start, end, intervals, x->left);
    }

    if (x->right!=NULL && x->right->min <= end)
    {
        search_iter(start, end, intervals, x->right);
    }
};

/**
 * Iterative method for print.
 */
void IntervalTree::print_iter(IntervalTreeNode* x)
{
    //intervals before target interval && after target interval
    if (x==NULL)
    {
        std::cerr << "NULL";
    }
    else
    {
        std::cerr << "[";
        print_iter(x->left);
        std::cerr << ",(" << x->start << ",";
        for (uint32_t i=0; i<x->intervals.size(); ++i)
        {
            std::cerr << x->intervals[i]->end << (i==x->intervals.size()-1 ? "" : ":");
        }
        std::cerr << "," << x->max << "," << x->color << "),";
        print_iter(x->right);
        std::cerr << "]";
    }
};

/**
 * Iterative method for validate.
 */
void IntervalTree::validate_iter(IntervalTreeNode* x, uint32_t depth)
{
    //intervals before target interval && after target interval
    if (x!=NULL)
    {
        //inspects color
        if (x->color==RED)
        {
            if (x->left!=NULL && x->left->color==RED)
            {
                std::cerr << "LEFT Color violation\n";
            }

            if (x->right!=NULL && x->right->color==RED)
            {
                std::cerr << "RIGHT Color violation\n";
            }
        }

        //inspects order
        if (x->left!=NULL && x->left->start > x->start)
        {
            std::cerr << "LEFT Order violation\n";
            std::cerr << "left->start: " <<  x->left->start << "\n";
            std::cerr << "this->start: " <<  x->start << "\n";
        }

        if (x->right!=NULL && x->right->start < x->start)
        {
            std::cerr << "RIGHT Order violation\n";
            std::cerr << "right->start: " <<  x->right->start << "\n";
            std::cerr << "this->start: " <<  x->start << "\n";
        }

        //inspects max
        for (uint32_t i=0; i<x->intervals.size(); ++i)
        {
            //overlap
            if (x->intervals[i]->end > x->max)
            {
                std::cerr << "CENTER Max violation\n";
                std::cerr << "intervals[" << i << "]->end : " <<  x->intervals[i]->end<< "\n";
                std::cerr << "this->max: " <<  x->max << "\n";
            }

            //overlap
            if (x->intervals[i]->start < x->min)
            {
                std::cerr << "CENTER Min violation\n";
                std::cerr << "intervals[" << i << "]->start : " <<  x->intervals[i]->start<< "\n";
                std::cerr << "this->min: " <<  x->min << "\n";
            }

            //incorrect value
            if (x->intervals[i]->start != x->start)
            {
                std::cerr << "VALUE violation\n";
                std::cerr << "intervals[" << i << "]->start : " <<  x->intervals[i]->start<< "\n";
                std::cerr << "this->start: " <<  x->start << "\n";
            }

        }

        if (x->left!=NULL && x->left->max > x->max)
        {
            std::cerr << "LEFT Max violation\n";
            std::cerr << "left->max: " <<  x->left->max << "\n";
            std::cerr << "this->max: " <<  x->max << "\n";
        }

        if (x->right!=NULL && x->right->max > x->max)
        {
            std::cerr << "RIGHT Max violation\n";
            std::cerr << "right->max: " <<  x->right->max << "\n";
            std::cerr << "this->max: " <<  x->max << "\n";
        }

        //overlap
        if (x->left!=NULL &&  x->left->min < x->min)
        {
            std::cerr << "LEFT Min violation\n";
            std::cerr << "left->min: " <<  x->left->min << "\n";
            std::cerr << "this->min: " <<  x->min << "\n";
        }

        //overlap
        if (x->right!=NULL &&  x->right->min < x->min)
        {
            std::cerr << "RIGHT Min violation\n";
            std::cerr << "right->min: " <<  x->right->min << "\n";
            std::cerr << "this->min: " <<  x->min << "\n";
        }

        validate_iter(x->left, depth+1);
        validate_iter(x->right, depth+1);
    }
    else
    {
        height = height<depth? depth : height;
    }
};

/**
 * Left rotates a node.
 */
void IntervalTree::left_rotate(IntervalTreeNode* x)
{
    //do nothing
    if (x->right==NULL)
        return;

    //move left child of y
    IntervalTreeNode* y = x->right;
    x->right = y->left;
    if (y->left != NULL)
    {
        (y->left)->parent = x;
    }

    //move y
    y->parent = x->parent;

    if (x->parent == NULL)
    {
        root = y;
    }
    else
    {
        if (x==(x->parent)->left)
        {
            (x->parent)->left = y;
        }
        else
        {
            (x->parent)->right = y;
        }
    }

    //move x
    y->left = x;
    x->parent = y;

    //update max for x
    if (x->right!=NULL)
    {
        x->max = std::max(x->right->max, x->max);
        x->min = std::min(x->right->min, x->min);
    }

    //update max for y
    y->max = std::max(y->left->max, y->max);
    y->min = std::min(y->left->min, y->min);
};

/**
 * Righ rotates a node.
 */
void IntervalTree::right_rotate(IntervalTreeNode* y)
{
    if (y->left==NULL)
        return;

    //move right child of x
    IntervalTreeNode* x = y->left;
    y->left = x->right;
    if (x->right != NULL)
    {
        (x->right)->parent = y;
    }

    //move x
    x->parent = y->parent;
    if (y->parent == NULL)
    {
        root = x;
    }
    else
    {
        if (y==(y->parent)->left)
        {
            (y->parent)->left = x;
        }
        else
        {
            (y->parent)->right = x;
        }
    }

    //move y
    x->right = y;
    y->parent = x;

    //update max for y
    if (y->left!=NULL)
    {
        y->max = std::max(y->left->max, y->max);
        y->min = std::min(y->left->min, y->min);
    }

    //update max for x
    x->max = std::max(x->right->max, x->max);
    x->min = std::min(x->right->min, x->min);
};