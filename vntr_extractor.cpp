/* The MIT License

   Copyright (c) 2016 Adrian Tan <atks@umich.edu>

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

#include "vntr_extractor.h"

VNTRExtractor::VNTRExtractor(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& ref_fasta_file)
{
 

    ////////////////////////
    //stats initialization//
    ////////////////////////
    
    ////////////////////////
    //tools initialization//
    ////////////////////////
    
}

/**
 * Inserts a VNTR record.
 * Returns true if successful.
 */
void VNTRExtractor::insert_variant_record_into_buffer(Variant* variant)
{
//    std::list<VNTR>::iterator i = vntr_buffer.begin();
//    while(i!=vntr_buffer.end())
//    {
//        VNTR& cvntr = *i;
//
//        if (vntr.rid > cvntr.rid)
//        {
//            vntr_buffer.insert(i, vntr);
//            return true;
//        }
//        else if (vntr.rid == cvntr.rid)
//        {
//            if (vntr.exact_beg1 > cvntr.exact_beg1)
//            {
//                vntr_buffer.insert(i, vntr);
//                return true;
//            }
//            else if (vntr.exact_beg1 == cvntr.exact_beg1)
//            {
//                if (vntr.exact_end1 > cvntr.exact_end1)
//                {
//                    vntr_buffer.insert(i, vntr);
//                    return true;
//                }
//                else if (cvntr.exact_end1 == vntr.exact_end1)
//                {
//                    if (cvntr.motif > vntr.motif)
//                    {
//                        vntr_buffer.insert(i, vntr);
//                        return true;
//                    }
//                    else if (cvntr.motif == vntr.motif)
//                    {
//                        //do not insert
//                        return false;
//                    }
//                    else // cvntr.motif > vntr.motif
//                    {
//                        ++i;
//                    }
//                }
//                else // cvntr.rend1 > vntr.rend1
//                {
//                    ++i;
//                }
//            }
//            else //vntr.rbeg1 < cvntr.rbeg1
//            {
//                ++i;
//            }
//        }
//        else //vntr.rid < cvntr.rid is impossible if input file is ordered.
//        {
//            fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
//            exit(1);
//        }
//    }
    
  
//    vntr_buffer.push_back(vntr);
  
}

/**
 * Flush variant buffer.
 */
void VNTRExtractor::flush_variant_buffer(Variant* var)
{
//    if (vntr_buffer.empty())
//    {
//        return;
//    }
//
//    int32_t rid = bcf_get_rid(v);
//    int32_t pos1 = bcf_get_pos1(v);
//
//    //search for vntr to start deleting from.
//    std::list<VNTR>::iterator i = vntr_buffer.begin();
//    while(i!=vntr_buffer.end())
//    {
//        VNTR& vntr = *i;
//
//        if (vntr.rid < rid)
//        {
//            break;
//        }
//        else if (vntr.rid == rid)
//        {
//            if (method=="e")
//            {
//                if (vntr.exact_end1 < pos1-2000)
//                {
//                    break;
//                }
//            }
//            else if (method=="f")
//            {
//                if (vntr.fuzzy_end1 < pos1-2000)
//                {
//                    break;
//                }
//            }
//        }
//        else //rid < vntr.rid is impossible
//        {
//            fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
//            exit(1);
//        }
//
//        ++i;
//    }
//
//    while (i!=vntr_buffer.end())
//    {
//        i = vntr_buffer.erase(i);
//    }
}

/**
 * Creates a VNTR record.
 */
void VNTRExtractor::create_vntr_record(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
{
//    VNTR& vntr = variant.vntr;
//
//    //shared fields
//    bcf_set_rid(v, variant.rid);
//    bcf_update_info_string(h, v, MOTIF.c_str(), vntr.motif.c_str());
//    bcf_update_info_string(h, v, BASIS.c_str(), vntr.basis.c_str());
//    bcf_update_info_string(h, v, RU.c_str(), vntr.ru.c_str());
//
//    //exact characteristics
//    bcf_update_info_float(h, v, EX_SCORE.c_str(), &vntr.exact_motif_concordance, 1);
//    bcf_update_info_float(h, v, EX_RL.c_str(), &vntr.exact_rl, 1);
//    bcf_update_info_float(h, v, EX_LL.c_str(), &vntr.exact_ll, 1);
//    int32_t ru_count[2] = {vntr.exact_no_exact_ru, vntr.exact_total_no_ru};
//    bcf_update_info_int32(h, v, EX_RU_COUNTS.c_str(), &ru_count, 2);
//    int32_t flank_pos1[2] = {variant.vntr.exact_beg1-1, variant.vntr.exact_end1+1};
//    bcf_update_info_int32(h, v, EX_FLANKS.c_str(), &flank_pos1, 2);
//
//    //fuzzy characteristics
//    bcf_update_info_float(h, v, FZ_SCORE.c_str(), &vntr.fuzzy_motif_concordance, 1);
//    bcf_update_info_float(h, v, FZ_RL.c_str(), &vntr.fuzzy_rl, 1);
//    bcf_update_info_float(h, v, FZ_LL.c_str(), &vntr.fuzzy_ll, 1);
//    int32_t fuzzy_ru_count[2] = {vntr.fuzzy_no_exact_ru, vntr.fuzzy_total_no_ru};
//    bcf_update_info_int32(h, v, FZ_RU_COUNTS.c_str(), &fuzzy_ru_count, 2);
//    int32_t fuzzy_flank_pos1[2] = {variant.vntr.fuzzy_beg1-1, variant.vntr.fuzzy_end1+1};
//    bcf_update_info_int32(h, v, FZ_FLANKS.c_str(), &fuzzy_flank_pos1, 2);
//
//    if (vntr.is_large_repeat_tract) bcf_update_info_flag(h, v, "LARGE_REPEAT_REGION", NULL, 1);
//
//    if (variant.vntr.definition_support=="e")
//    {
//        bcf_update_info_flag(h, v, "EXACT", NULL, 1);
//
//        //VNTR position and sequences
//        bcf_set_pos1(v, vntr.exact_beg1);
//        s.l = 0;
//        kputs(vntr.exact_repeat_tract.c_str(), &s);
//        kputc(',', &s);
//        kputs("<VNTR>", &s);
//        bcf_update_alleles_str(h, v, s.s);
//
//        //update flanking sequences
//        if (add_flank_annotation)
//        {
//            update_flankseq(h, v, variant.chrom.c_str(),
//                            variant.vntr.exact_beg1-10, variant.vntr.exact_beg1-1,
//                            variant.vntr.exact_end1+1, variant.vntr.exact_end1+10);
//        }
//    }
//    else if (variant.vntr.definition_support=="f")
//    {
//        bcf_update_info_flag(h, v, "FUZZY", NULL, 1);
//
//        //VNTR position and sequences
//        bcf_set_pos1(v, vntr.fuzzy_beg1);
//        s.l = 0;
//        kputs(vntr.fuzzy_repeat_tract.c_str(), &s);
//        kputc(',', &s);
//        kputs("<VNTR>", &s);
//        bcf_update_alleles_str(h, v, s.s);
//
//        //update flanking sequences
//        if (add_flank_annotation)
//        {
//            update_flankseq(h, v, variant.chrom.c_str(),
//                            variant.vntr.fuzzy_beg1-10, variant.vntr.fuzzy_beg1-1,
//                            variant.vntr.fuzzy_end1+1, variant.vntr.fuzzy_end1+10);
//        }
//    }
//
//    //individual fields - just set GT
//    bcf_update_genotypes(h, v, gts, no_samples);
}



