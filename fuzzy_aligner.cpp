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

#include "fuzzy_aligner.h"

/**
 * Align and compute genotype likelihood.
 */
void FuzzyAligner::align(std::string& sequence, std::string& repeat_unit)
{
    //use shifting idea
//    
//    AAGGAAGGAAGGAAGGAAGGAAGG
//    AAGG----AAGGAAGGAAGGAAGG
//    
//    
//    1. find an occurrence of a motif
//    2. delete it
//    3. fuzzy shift it
//    
//    
//    1. shift without any idea of motif
//    2. shift with idea of motif
//    3. 
//    
//    A. fuzzy left align.
//    while ends the same
//        shift
//        if stuck
//             
//             track - score
//                     exact units
//                     inexact units
//                     position of sequence
//                     position of deleted sequence
//             
//             how many to keep track
//             how many to store?
//             how to prune
//             
//                
//     AAGGAAGGAAGGAAGGAAGGAAGG
//     AAGGAAGGG----AGGAAGGAAGG
//    
//         AAGGAAGGAAGG AAGG AAGG AAGG 1 sub
//         AAGGAAGG---- GAGG AAGG AAGG
//         
//         AAGGAAGGAAGG -AAGG AAGG AAGG  1 ins
//         AAGGAAG----  GGAGG AAGG AAGG
//
//                AAGGAAGGAAG G -AAGG AAGG AAGG  1 ins
//                AAGGAA----  G  GGAGG AAGG AAGG
//
//
//         AAGGAAGGAAGG AAGG AAGG AAGG  1 del
//         AAGGAAGG---- -AGG  AAGG  AAGG
    
}