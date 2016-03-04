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
    //////////////////////
    //i/o initialization//
    //////////////////////
    buffer_window_allowance = 5000;

    this->input_vcf_file = input_vcf_file;
    this->output_vcf_file = output_vcf_file;
    odr = new BCFOrderedReader(input_vcf_file, intervals);
    odw = new BCFOrderedWriter(output_vcf_file, buffer_window_allowance);
    odw->link_hdr(odr->hdr);
    odw->write_hdr();

    //for adding empty genotype fields for a VCF file with individual information
    int32_t no_samples = bcf_hdr_nsamples(odw->hdr);
    gts = NULL;
    if (no_samples)
    {    
        gts = (int32_t*) malloc(no_samples*sizeof(int32_t));
        for (uint32_t i=0; i<no_samples; ++i)
        {
            gts[i] = 0;
        }
    }

    ////////////////////////////////////////////
    //add relevant field for adding VNTR records
    ////////////////////////////////////////////    
    bcf_hdr_append(odw->hdr, "##ALT=<ID=VNTR,Description=\"Variable Number of Tandem Repeats.\">");
    bcf_hdr_append(odw->hdr, "##INFO=<ID=ASSOCIATED_INDEL,Number=.,Type=String,Description=\"Indels that were annotated as this VNTR.\">");
    bool rename = false;
    MOTIF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "MOTIF", "1", "String", "Canonical motif in a VNTR.", rename);
    RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU", "1", "String", "Repeat unit in the reference sequence.", rename);
    BASIS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "BASIS", "1", "String", "Basis nucleotides in the motif.", rename);
    MLEN = bcf_hdr_append_info_with_backup_naming(odw->hdr, "MLEN", "1", "Integer", "Motif length.", rename);
    BLEN = bcf_hdr_append_info_with_backup_naming(odw->hdr, "BLEN", "1", "Integer", "Basis length.", rename);
    REPEAT_TRACT = bcf_hdr_append_info_with_backup_naming(odw->hdr, "REPEAT_TRACT", "2", "Integer", "Boundary of the repeat tract detected by exact alignment.", rename);
    COMP = bcf_hdr_append_info_with_backup_naming(odw->hdr, "COMP", "4", "Integer", "Composition(%) of bases in an exact repeat tract.", rename);
    ENTROPY = bcf_hdr_append_info_with_backup_naming(odw->hdr, "ENTROPY", "1", "Float", "Entropy measure of an exact repeat tract [0,2].", rename);
    ENTROPY2 = bcf_hdr_append_info_with_backup_naming(odw->hdr, "ENTROPY2", "1", "Float", "Dinucleotide entropy measure of an exact repeat tract [0,4].", rename);
    KL_DIVERGENCE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "KL_DIVERGENCE", "1", "Float", "Kullback-Leibler Divergence of an exact repeat tract.", rename);
    KL_DIVERGENCE2 = bcf_hdr_append_info_with_backup_naming(odw->hdr, "KL_DIVERGENCE2", "1", "Float", "Dinucleotide Kullback-Leibler Divergence of an exact repeat tract.", rename);
    RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RL", "1", "Integer", "Reference exact repeat tract length in bases.", rename);
    LL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "LL", "1", "Integer", "Longest exact repeat tract length in bases.", rename);
    RU_COUNTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU_COUNTS", "2", "Integer", "Number of exact repeat units and total number of repeat units in exact repeat tract.", rename);
    SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "SCORE", "1", "Float", "Score of repeat unit in exact repeat tract.", rename);
    TRF_SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "TRF_SCORE", "1", "Integer", "TRF Score for M/I/D as 2/-7/-7 in exact repeat tract.", rename);

   

//used in classification
//        std::string TR = bcf_hdr_append_info_with_backup_naming(odw->hdr, "TR", "1", "String", \"Tandem repeat associated with this indel.", rename);
//
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=LARGE_REPEAT_REGION,Number=0,Type=Flag,Description=\"Very large repeat region, vt only detects up to 1000bp long regions.\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=FLANKSEQ,Number=1,Type=String,Description=\"Flanking sequence 10bp on either side of detected repeat region.\">");
//

    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_variants = 0;
    no_added_vntrs = 0;

    ////////////////////////
    //tools initialization//
    ////////////////////////
    refseq = new ReferenceSequence(ref_fasta_file);

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
    Variant& new_variant = *variant;
    
    //create VNTR record
    if (new_variant.type&VT_INDEL)
    {
//        create_variant_record(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
    }    
    
    std::list<Variant*>::iterator i = variant_buffer.begin();
    while(i != variant_buffer.end())
    {
        Variant& current_variant = **i;

        if (new_variant.rid > current_variant.rid)
        {
            variant_buffer.insert(i, &new_variant);
        }
        else if (new_variant.rid == current_variant.rid)
        {
            if (new_variant.beg1 > current_variant.beg1)
            {
                variant_buffer.insert(i, &new_variant);
            }
            else if (new_variant.beg1 == current_variant.beg1)
            {
                if (new_variant.end1 > current_variant.end1)
                {
                    variant_buffer.insert(i, &new_variant);
                }
                else if (current_variant.end1 == new_variant.end1)
                {
                    if (new_variant.type==VT_VNTR && current_variant.type==VT_VNTR)
                    {
                        if (current_variant.vntr.motif > new_variant.vntr.motif)
                        {
                            variant_buffer.insert(i, &new_variant);
                        }
                        else if (current_variant.vntr.motif == new_variant.vntr.motif)
                        {
                            bcf1_t *v = new_variant.vs[0];
                            current_variant.vs.push_back(v);
                            current_variant.indel_vs.push_back(v);
                            
                            //update current_variant
                            ++ current_variant.no_overlapping_indels;
                            
                            //do not insert
//                            return false;
                        }
                        else // current_variant.motif > new_variant.motif
                        {
                            ++i;
                        }
                    }
                }
                else // current_variant.rend1 > new_variant.rend1
                {
                    ++i;
                }
            }
            else //new_variant.rbeg1 < current_variant.rbeg1
            {
                ++i;
            }
        }
        else //new_variant.rid < current_variant.rid is impossible if input file is ordered.
        {
            fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
            exit(1);
        }
    }


   variant_buffer.push_back(variant);

}

/**
 * Flush variant buffer.
 */
void VNTRExtractor::flush_variant_buffer(Variant* var)
{
//    if (variant_buffer.empty())
//    {
//        return;
//    }
//
//    int32_t rid = bcf_get_rid(v);
//    int32_t pos1 = bcf_get_pos1(v);
//
//    //search for variant to start deleting from.
//    std::list<VNTR>::iterator i =variant_buffer.begin();
//    while(i!=variant_buffer.end())
//    {
//        VNTR& variant = *i;
//
//        if (variant.rid < rid)
//        {
//            break;
//        }
//        else if (variant.rid == rid)
//        {
//            if (method=="e")
//            {
//                if (variant.exact_end1 < pos1-2000)
//                {
//                    break;
//                }
//            }
//            else if (method=="f")
//            {
//                if (variant.fuzzy_end1 < pos1-2000)
//                {
//                    break;
//                }
//            }
//        }
//        else //rid < variant.rid is impossible
//        {
//            fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
//            exit(1);
//        }
//
//        ++i;
//    }
//
//    while (i!=variant_buffer.end())
//    {
//        i =variant_buffer.erase(i);
//    }
}

/**
 * Creates a VNTR record.
 */
Variant* VNTRExtractor::create_vntr_record(bcf_hdr_t* h, bcf1_t *v)
{
//    if ()
//    {
//        VNTR& variant = variant.variant;
//    
//        bcf1_t* v = bcf_init();
//    
//    
//        //VNTR position and sequences
//        bcf_set_pos1(v, variant.exact_beg1);
//        s.l = 0;
//        kputs(variant.exact_repeat_tract.c_str(), &s);
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
//    
//        //shared fields
//        
//        bcf_set_rid(v, variant.rid);
//        bcf_update_info_string(h, v, MOTIF.c_str(), variant.motif.c_str());
//        bcf_update_info_string(h, v, BASIS.c_str(), variant.basis.c_str());
//        bcf_update_info_string(h, v, RU.c_str(), variant.ru.c_str());
//    
//        bcf_update_info_float(h, v, SCORE.c_str(), &variant.exact_motif_concordance, 1);
//        bcf_update_info_float(h, v, RL.c_str(), &variant.exact_rl, 1);
//        bcf_update_info_float(h, v, LL.c_str(), &variant.exact_ll, 1);
//        int32_t ru_count[2] = {variant.exact_no_exact_ru, variant.exact_total_no_ru};
//        bcf_update_info_int32(h, v, RU_COUNTS.c_str(), &ru_count, 2);
//        int32_t flank_pos1[2] = {variant.variant.exact_beg1-1, variant.variant.exact_end1+1};
//        bcf_update_info_int32(h, v, FLANKS.c_str(), &flank_pos1, 2);
//    
//    
////        std::string MOTIF;
////        std::string RU;
////        std::string BASIS;
////        std::string MLEN;
////        std::string BLEN;
////        std::string REPEAT_TRACT;
////        std::string COMP;
////        std::string ENTROPY;
////        std::string ENTROPY2;
////        std::string KL_DIVERGENCE;
////        std::string KL_DIVERGENCE2;
////        std::string RL;
////        std::string LL;
////        std::string RU_COUNTS;
////        std::string SCORE;
////        std::string TRF_SCORE;
//    
//    
//            
//    
//        if (no_samples) bcf_update_genotypes(h, v, gts, no_samples);
//    }
}



