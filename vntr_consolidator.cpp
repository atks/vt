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

#include "vntr_consolidator.h"

/**
 * Constructor.
 */
VNTRConsolidator::VNTRConsolidator(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& ref_fasta_file)
{
    //////////////////////
    //i/o initialization//
    //////////////////////
    this->input_vcf_file = input_vcf_file;
    this->output_vcf_file = output_vcf_file;
    odr = new BCFOrderedReader(input_vcf_file, intervals);
    odw = new BCFOrderedWriter(output_vcf_file, 3000);
    odw->link_hdr(odr->hdr);
    bcf_hdr_append(odw->hdr, "##FILTER=<ID=shorter_vntr,Description=\"Another VNTR overlaps with this VNTR.\">");
    //to be removed later
    bcf_hdr_append(odw->hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant.\">");
    bcf_hdr_append(odw->hdr, "##INFO=<ID=TRF_SCORE,Number=1,Type=Integer,Description=\"TRF score if the tandem repeat.\">");
    odw->write_hdr();

    overlap_vntr = const_cast<char*>("overlap_vntr");
    overlap_vntr_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr");

    buffer_window_allowance = 5000;

    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_total_variants = 0;
    no_overlap_vntrs = 0;
    no_dropped_vntrs = 0;

    no_snps = 0;
    no_indels = 0;
    no_vntrs = 0;
    no_other_variants = 0;

    //VNTR types
    no_isolated_exact_vntrs = 0;

    no_perfect_concordance_isolated_exact_vntrs = 0;
    no_imperfect_concordance_isolated_exact_vntrs = 0;

    no_perfect_concordance_isolated_inexact_vntrs = 0;
    no_imperfect_concordance_isolated_inexact_vntrs = 0;

    no_isolated_inexact_vntrs = 0;

    no_clustered_exact_vntrs = 0;
    no_clustered_inexact_vntrs = 0;

    no_isolated_complete_overlap_vntrs = 0;
    no_isolated_incomplete_overlap_vntrs = 0;

    no_isolated_partial_overlap_vntrs = 0;
    no_isolated_no_overlap_vntrs = 0;

    ////////////////////////
    //tools initialization//
    ////////////////////////
    refseq = new ReferenceSequence(ref_fasta_file);
    cmp = new CandidateMotifPicker(debug);
    fd = new FlankDetector(ref_fasta_file, debug);
}

/**
 * Update distribution of overlapping VNTRs
 */
void VNTRConsolidator::update_overlapping_vntr_hist(int32_t no_overlapping_vntrs)
{
    if (overlapping_vntr_hist.size()<no_overlapping_vntrs+1)
    {
        for (uint32_t i=overlapping_vntr_hist.size(); i<no_overlapping_vntrs+1; ++i)
        {
            overlapping_vntr_hist.push_back(0);
        }
    }

    ++overlapping_vntr_hist[no_overlapping_vntrs];
}


/**
 * Inserts a Variant record.
 *
 * VNTR overlapping variants are handled by adding to an existing VNTR record in the vntr_vs vector for consolidation purposes.
 * All other variants are simply added.
 */
void VNTRConsolidator::insert_variant_record_into_buffer(Variant* variant)
{
    std::list<Variant *>::iterator i = variant_buffer.begin();

    if (variant->type==VT_SNP)
    {
        ++no_snps;
    }
    else if (variant->type==VT_INDEL)
    {
        ++no_indels;
    }
    else if (variant->type==VT_VNTR)
    {
        ++no_vntrs;
    }
    else
    {
        ++no_other_variants;
    }

    while(i!=variant_buffer.end())
    {
        Variant *cvariant = *i;

        if (variant->rid > cvariant->rid)
        {
            break;
        }
        else if (variant->rid == cvariant->rid)
        {
            if (variant->end1 < cvariant->beg1) //not possible
            {
                fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
                exit(1);
            }
            //after most recent variant, we need to have the buffer window allowance, because the variants are roughly
            //ordered by start.  It is possible to have the start positions changed when merging VNTRs
            //resulting in unordered variants.
            else if (variant->beg1 > cvariant->end1 + buffer_window_allowance)
            {
                break;
            }
            else if (variant->end1 >= cvariant->beg1 && variant->beg1 <= cvariant->end1) //overlaps
            {
                if (variant->type==VT_VNTR && cvariant->type==VT_VNTR)
                {
                    bcf1_t* v = variant->v;
                    //this will induce order change thus the 1000bp buffer in the prior conditional
                    cvariant->beg1 = std::min(cvariant->beg1, variant->beg1);
                    cvariant->end1 = std::max(cvariant->end1, variant->end1);
                    cvariant->vs.push_back(v);
                    cvariant->vntr_vs.push_back(v);
                    ++cvariant->no_overlapping_vntrs;

                    return;
                }

                ++i;
            }
            else
            {
                ++i;
            }
        }
        else //variant.rid < cvariant.rid is impossible if input file is ordered.
        {
            fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
            exit(1);
        }
    }

    //push all variants
    variant_buffer.push_front(variant);
}

/**
 * Flush variant buffer.
 */
void VNTRConsolidator::flush_variant_buffer(Variant* var)
{
    if (variant_buffer.empty())
    {
        return;
    }

    int32_t rid = var->rid;
    int32_t beg1 = var->beg1;

    while (!variant_buffer.empty())
    {
        Variant* variant = variant_buffer.back();

        if (variant->rid < rid)
        {
            if (variant->type==VT_VNTR)
            {
                if (consolidate_multiple_overlapping_vntrs(variant))
                {
                    odw->write(variant->v);
                    variant->v = NULL;
                    delete variant;
                    variant_buffer.pop_back();
                }
            }
            else
            {
                odw->write(variant->v);
                variant->v = NULL;
                delete variant;
                variant_buffer.pop_back();
            }
        }
        else if (variant->rid == rid)
        {
            if (variant->beg1 < beg1-buffer_window_allowance)
            {
                if (variant->type==VT_VNTR)
                {
                    if (consolidate_multiple_overlapping_vntrs(variant))
                    {


                    }

                    delete variant;
                    variant_buffer.pop_back();

                }
                else
                {
                    odw->write(variant->v);
                    variant->v = NULL;
                    delete variant;
                    variant_buffer.pop_back();
                }
            }
            else
            {
                break;
            }
        }
    }
}



/**
 * Consolidate multiallelic variant based on associated biallelic records
 * stored in vs.  Updates v which is to be the consolidated multiallelic
 * variant.
 *
 * returns true if the multiallelic variant is good to go.
 */
bool VNTRConsolidator::consolidate_multiple_overlapping_vntrs(Variant* variant)
{
    update_overlapping_vntr_hist(variant->no_overlapping_vntrs);

    if (variant->no_overlapping_vntrs==0)
    {
        if ( debug)
        {
            std::cerr  << "################\n";
            std::cerr  << "#1 isolated VNTR\n";
            std::cerr  << "################\n";
            std::cerr << "no overlapping SNPs   " << variant->no_overlapping_snps << "\n";
            std::cerr << "no overlapping Indels " << variant->no_overlapping_indels << "\n";
            std::cerr << "no overlapping VNTRs  " << variant->no_overlapping_vntrs << "\n";
            std::cerr << "consolidating: " << variant->vs.size() << " alleles\n";

//                bcf_print(odw->hdr, variant->v);
        }

        bcf1_t* vntr_v = variant->v;


    }
    else if (variant->no_overlapping_vntrs >= 1)
    {
        if ((true && variant->vntr_vs.size()<6) || debug)
        {
            std::cerr  << "###################################\n";
            std::cerr  << "#2 or more VNTR and multiple Indels\n";
            std::cerr  << "###################################\n";
            std::cerr << "no overlapping SNPs   " << variant->no_overlapping_snps << "\n";
            std::cerr << "no overlapping Indels " << variant->no_overlapping_indels << "\n";
            std::cerr << "no overlapping VNTRs  " << variant->no_overlapping_vntrs << "\n";
            std::cerr << "consolidating: " << (variant->vs.size()+1) << " alleles\n";

            detect_consistent_motifs(variant);


        }

        return false;
    }

    return false;
}

/**
 * Detects a a consistent basis motif in a chain of overlapping VNTRs.
 */
bool VNTRConsolidator::detect_consistent_motifs(Variant* variant)
{
    if (variant->vntr_vs.size()>1)
    {
        VNTR& vntr = variant->vntr;

        std::map<std::string, int32_t> basis_map;
        std::map<std::string, int32_t>::iterator it;

        int32_t merged_beg1 = vntr.beg1;
        int32_t merged_end1 = vntr.end1;

        std::map<std::string, int32_t> motifs;

        ///////////////////
        //merge the regions
        ///////////////////
        for (uint32_t i=0; i<variant->vntr_vs.size(); ++i)
        {
            bcf1_t* vntr_v = variant->vntr_vs[i];
            std::string cbasis = vntr.basis;

            std::string motif = bcf_get_info_str(odr->hdr, vntr_v, "MOTIF");;
            float concordance = bcf_get_info_flt(odr->hdr, vntr_v, "CONCORDANCE");
            std::vector<int32_t> repeat_tract = bcf_get_info_int_vec(odr->hdr, vntr_v, "REPEAT_TRACT");

            cbasis = vntr.get_basis(motif);
            std::cerr << (i+1) << ") " << motif << " (" << cbasis << ")\t" << concordance << "\t" << repeat_tract[0] << "," << repeat_tract[1] << "\t" << "\n";
            std::cerr << "\t" << bcf_get_ref(vntr_v) << "\n";
            bcf_print(odw->hdr, vntr_v);

            merged_beg1 = std::min(merged_beg1, repeat_tract[0]);
            merged_end1 = std::max(merged_end1, repeat_tract[1]);

            std::string current_motif(motif);
            motifs[current_motif] = 1;

            if ((it = basis_map.find(cbasis)) != basis_map.end())
            {
                ++it->second;
            }
            else
            {
                basis_map[cbasis] = 1;
            }
        }

        //clear priority queue
        ///////////////////
        //merge the regions
        ///////////////////
        while (!ordered_basis.empty()) ordered_basis.pop();

        int32_t n = variant->vntr_vs.size();
        it = basis_map.begin();
        while (it!=basis_map.end())
        {
//            std::cerr << "pqueue " << it->first << ", "  << it->second << "(" << ((float) it->second/n) << ")\n";
            basis_proportion bp = {it->first, (float) it->second/n};
            ordered_basis.push(bp);
            ++it;
        }

//        std::cerr << "size of ordered_bp " << ordered_basis.size() << "\n";

        basis_proportion top_bp = ordered_basis.top();
        if (top_bp.proportion>=1)
        {
            bcf1_t *new_v = bcf_dup(variant->v);

            std::string repeat_tract;
            refseq->fetch_seq(variant->chrom.c_str(), merged_beg1, merged_end1, repeat_tract);

//            std::cerr << "\n";
//            std::cerr << "\n";
//            std::cerr << "\n";
//            std::cerr << "merged VNTR\n";
//
//            bcf_print(odw->hdr, new_v);
//            std::cerr << "newly merged ref " << variant->chrom << ":" << merged_beg1 << "-" << merged_end1 << " " << repeat_tract << "\n";

            //collect motifs from overlpaping records
            std::string best_motif = "";
            float best_motif_concordance = 0;
            std::map<std::string, int32_t>::iterator it;
            for (it = motifs.begin(); it!=motifs.end(); ++it)
            {
                std::string motif = it->first;
                fd->compute_purity_score(repeat_tract, motif);

                if (fd->score > best_motif_concordance)
                {
                    best_motif_concordance = fd->score;
                    best_motif = it->first;
                }
            }

//            fd->polish_repeat_tract_ends(repeat_tract, best_motif, true);


            //recompute conconcordance with polished tract
//            fd->compute_purity_score(fd->polished_repeat_tract, best_motif);

            //VNTR position and sequences
            bcf_set_pos1(new_v, merged_beg1+fd->min_beg0);
            kstring_t s = {0,0,0};
            s.l = 0;
            kputs(fd->polished_repeat_tract.c_str(), &s);
            kputc(',', &s);
            kputs("<VNTR>", &s);
            bcf_update_alleles_str(odr->hdr, new_v, s.s);
            if (s.m) free(s.s);

            //VNTR motif
            bcf_update_info_string(odw->hdr, new_v, "MOTIF", best_motif.c_str());
            bcf_update_info_string(odw->hdr, new_v, "RU", fd->ru.c_str());

            //VNTR characteristics
            bcf_update_info_float(odw->hdr, new_v, "CONCORDANCE", &fd->score , 1);
            bcf_update_info_float(odw->hdr, new_v, "RL", &fd->rl, 1);
            int32_t new_fuzzy_flanks[2] = {merged_beg1+fd->min_beg0-1, merged_beg1+fd->min_beg0+(int32_t)fd->polished_repeat_tract.size()};
            bcf_update_info_int32(odw->hdr, new_v, "FLANKS", &new_fuzzy_flanks, 2);

//            //flank positions
//            int32_t fuzzy_flank_pos1[2] = {vntr.fuzzy_beg1-1, vntr.fuzzy_end1+1};
//            bcf_update_info_int32(h, v, "FZ_FLANKS", &fuzzy_flank_pos1, 2);
//            int32_t ru_count[2] = {vntr.fuzzy_no_exact_ru, vntr.fuzzy_total_no_ru};
//            bcf_update_info_int32(h, v, "FZ_RU_COUNTS", &ru_count, 2);

//              bcf_print(odw->hdr, new_v);


            return true;
        }
        else
        {
            return false;
        }
    }


    return false;
}

/**
 * Flush variant buffer.
 */
void VNTRConsolidator::flush_variant_buffer()
{
    while (!variant_buffer.empty())
    {
        Variant* variant = variant_buffer.back();

        if (variant->type==VT_VNTR)
        {
            if (consolidate_multiple_overlapping_vntrs(variant))
            {
                odw->write(variant->v);
                variant->v = NULL;
            }

            delete variant;
            variant_buffer.pop_back();
        }
        else
        {
            odw->write(variant->v);
            variant->v = NULL;
            delete variant;
            variant_buffer.pop_back();
        }
    }
}

/**
 * Close.
 */
void VNTRConsolidator::close()
{
    odw->close();
    odr->close();
}