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

#ifndef VNTR_H
#define VNTR_H

#include "utils.h"

//VNTR update modes
#define FINAL  1 // update final attributes
#define EXACT  2 // update exact attributes
#define FUZZY  4 // update fuzzy attributes

/**
 * Class for representing a VNTR.
 *
 * Includes functions specific to tandem repeats.
 */
class VNTR
{
    public:

    //chromosome
    int32_t rid;  //rid, redundant data with Variant. todo: something about this.

    std::string definition_support; //either exact or fuzzy

    ///////////////////////////
    //decided upon repeat tract
    ///////////////////////////
    std::string motif;                //motif                                                                                  
    std::string ru;                   //repeat unit on the reference                                                           
    int32_t mlen;                     //length of motif                                                                        
    std::string basis;                //unique bases found in motif                                                            
    int32_t blen;                     //basis length                                                                           
    std::string repeat_tract;         //repeat tract                                                                           
    int32_t comp[4];                  //composition of bases in repeat tract                                                   
    float entropy ;                   //sequence entropy of repeat tract                                                       
    float entropy2;                   //dinucleotide sequence entropy of repeat tract                                          
    float kl_divergence;              //Kullback-Leibler divergence of repeat tract                                            
    float kl_divergence2;             //dinucleotide Kullback-Leibler divergence of repeat tract                               
    int32_t beg1;                     //beginning of repeat tract                                                              
    int32_t end1;                     //end of repeat tract                                                                    
    float ref;                        //number of repeat units in reference repeat tract                                       
    float lref;                       //number of repeat units in reference repeat tract including longest allele              
    int32_t rl;                       //length of repeat tract in base pairs                                                 
    int32_t ll;                       //length of repeat tract (including longest alternate allele) in base pairs            
    int32_t no_exact_ru;              //number exact repeat units from hmm                                                     
    int32_t total_no_ru;              //total no of repeat units from hmm                                                      
    float score;                      //motif concordance from hmm                                                             
    int32_t trf_score;                //TRF score of exact repeat tract                                                        
    std::string lflank;               //left flank                                                                             
    std::string rflank;               //right flank                                                                            

    ////////////////////
    //exact repeat tract
    ////////////////////
    std::string exact_motif;          //motif
    std::string exact_ru;             //repeat unit on the reference
    int32_t exact_mlen;               //length of motif
    std::string exact_basis;          //unique bases found in motif
    int32_t exact_blen;               //basis length
    std::string exact_repeat_tract;   //repeat tract
    int32_t exact_comp[4];            //composition of bases in repeat tract
    float exact_entropy ;             //sequence entropy of repeat tract
    float exact_entropy2 ;            //dinucleotide sequence entropy of repeat tract
    float exact_kl_divergence;        //Kullback-Leibler divergence of repeat tract
    float exact_kl_divergence2;       //dinucleotide Kullback-Leibler divergence of repeat tract
    int32_t exact_beg1;               //beginning of repeat tract
    int32_t exact_end1;               //end of repeat tract
    float exact_ref;                  //number of repeat units in exact reference repeat tract
    float exact_lref;                 //number of repeat units in exact reference repeat tract including longest allele
    int32_t exact_rl;                 //length of repeat tract in base pairs
    int32_t exact_ll;                 //length of repeat tract (including longest alternate allele) in base pairs
    int32_t exact_no_exact_ru;        //number exact repeat units from hmm
    int32_t exact_total_no_ru;        //total no of repeat units from hmm
    float exact_score;                //motif concordance from hmm
    int32_t exact_trf_score;          //TRF score of exact repeat tract
    std::string exact_lflank;         //left flank
    std::string exact_rflank;         //right flank

    bool exact_ru_ambiguous;

    /////////////////
    //fuzzy alignment
    /////////////////
    std::string fuzzy_motif;          //motif
    std::string fuzzy_ru;             //repeat unit on the reference
    int32_t fuzzy_mlen;               //length of motif
    std::string fuzzy_basis;          //unique bases found in motif
    int32_t fuzzy_blen;               //basis length
    std::string fuzzy_repeat_tract;   //repeat tract
    int32_t fuzzy_comp[4];            //composition of bases in repeat tract
    float fuzzy_entropy ;             //sequence entropy of repeat tract
    float fuzzy_entropy2;             //dinucleotide sequence entropy of repeat tract
    float fuzzy_kl_divergence;        //Kullback-Leibler divergence of repeat tract
    float fuzzy_kl_divergence2;       //dinucleotide Kullback-Leibler divergence of repeat tract    
    int32_t fuzzy_beg1;               //beginning of repeat tract
    int32_t fuzzy_end1;               //end of repeat tract
    float fuzzy_ref;                  //number of repeat units in fuzzy reference repeat tract
    float fuzzy_lref;                 //number of repeat units in fuzzy reference repeat tract including longest allele
    int32_t fuzzy_rl;                 //length of fuzzy repeat tract in base pairs
    int32_t fuzzy_ll;                 //length of fuzzy repeat tract (including longest alternate allele) in base pairs
    int32_t fuzzy_no_exact_ru;        //number exact repeat units from fuzzy alignment
    int32_t fuzzy_total_no_ru;        //total no of repeat units from fuzzy alignment
    float fuzzy_score;                //motif concordance from hmm
    int32_t fuzzy_trf_score;          //TRF score of fuzzy repeat tract
    std::string fuzzy_lflank;         //left flank
    std::string fuzzy_rflank;         //right flank

    //types of repeat
    bool is_large_repeat_tract;
    bool is_interspersed_repeat_tract;
    
    //associated indels
    std::map<std::string, int32_t> associated_indels;        

    /**
     * Constructor.
     */
    VNTR();

    /**
     * Clear object.
     */
    void clear();

    /**
     * Adds an associated indel.
     */
    void add_associated_indel(std::string& indel);
            
    /**
     * Get associated indels.
     */
    std::string get_associated_indels();
    
    /**
     * Checks for equality.
     */
    bool equals(VNTR& vntr);
    
    #define complement(b) ("TGNCNNNNNA"[((b)-65)>>1])
    
    /**
     * Reverse complement a sequence.
     */
    static std::string reverse_complement(std::string& seq);

    /**
     * Return the canonical representation of a motif.
     * Considers reverse complements too.
     */
    static std::string canonicalize2(std::string& motif);
            
    /**
     * Return the canonical representation of a motif.
     */
    static std::string canonicalize(std::string& motif);
    
    /**
     * Checks if a string is periodic.
     *
     * Returns the length of the sub motif.
     * and returns 0 if the motif is periodic.
     */
    static int32_t is_periodic(std::string& motif);
    
    /**
     * Checks if a string is aperiodic.
     */
    bool is_aperiodic(std::string& motif);
            
    /**
     * Return the string of unique bases in a motif.
     */
    static std::string get_basis(std::string& motif);

    /**
     * Return the string of unique bases in a motif.
     */
    static std::string get_basis(char* motif, uint32_t n);

    /**
     * Shifts a string.
     */
    static std::string shift_str(std::string& seq, uint32_t i);

    /**
     * Print object.
     */
    void print();
};
#endif


//sub canonicalize
//{
//    my $motif = shift;
//
//    my $cmotif = $motif;
//
//    for my $i (0 .. (length($motif)-1))
//    {
//        my $smotif = shift_str($motif, $i);
//        my $rc_smotif = reverse_complement($smotif);
//
//        if (($smotif cmp  $cmotif) < 0)
//        {
//            $cmotif = $smotif;
//        }
//
//        if (($rc_smotif cmp  $cmotif) < 0)
//        {
//            $cmotif = $rc_smotif;
//        }
//    }
//
//    return $cmotif;
//}
//
//sub shift_str
//{
//    my ($motif, $i, $rest) = @_;
//
//    return substr($motif, $i) .  substr($motif, 0, $i);
//}
//
//sub is_periodic
//{
//    my $motif = shift;
//
//    for my $i (1 .. (length($motif)-1))
//    {
//        my $smotif = shift_str($motif, $i);
//
//        if ($smotif eq $motif)
//        {
//#            print "$motif => $smotif => $i\n";
//
//            return $i;
//        }
//    }
//
//#    print "$motif => 0\n";
//
//    return 0;
//}
//
//sub get_aperiodic_submotif
//{
//    my $motif = shift;
//
//    for my $i (1 .. (length($motif)-1))
//    {
//        my $smotif = shift_str($motif, $i);
//
//        if ($smotif eq  $motif)
//        {
//            return substr($motif, 0, $i);
//        }
//    }
//}
//
//sub get_aperiodic_submotif_by_index
//{
//    my ($motif, $i) = @_;
//
//#    print "REURN BY INDEX : $i\n";
//
//    return substr($motif, 0, $i);
//}
//
//sub complement
//{
//    my $b = shift;
//
//    if ($b eq "A")
//    {
//        return "T";
//    }
//    elsif ($b eq "C")
//    {
//        return "G";
//    }
//    elsif ($b eq "G")
//    {
//        return "C";
//    }
//    elsif ($b eq "T")
//    {
//        return "A";
//    }
//}
//
//sub reverse_complement
//{
//    my $seq = shift;
//
//    my $rc = "";
//
//    my $i = length($seq) - 1;
//    while ($i>=0)
//    {
//        my $b = substr($seq, $i, 1);
//        $rc .= complement($b);
//
//        --$i;
//    }
//
//    return $rc;
//}