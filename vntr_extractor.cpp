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

VNTRExtractor::VNTRExtractor(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& fexp, int32_t vntr_classification_code, std::string& ref_fasta_file)
{
    this->vntr_classification = vntr_classification_code;

    //////////////////////
    //i/o initialization//
    //////////////////////
    buffer_window_allowance = 5000;

    this->input_vcf_file = input_vcf_file;
    this->output_vcf_file = output_vcf_file;
    odr = new BCFOrderedReader(input_vcf_file, intervals);
    odw = new BCFOrderedWriter(output_vcf_file, 2*buffer_window_allowance);
    odw->link_hdr(odr->hdr);

    //for adding empty genotype fields for a VCF file with individual information
    no_samples = bcf_hdr_nsamples(odw->hdr);
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
    END = bcf_hdr_append_info_with_backup_naming(odw->hdr, "END", "1", "Integer", "End position of the variant.", rename);
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
    ASSOCIATED_INDEL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "ASSOCIATED_INDEL", ".", "String", "Indels that induced this VNTR.", rename);
    odw->write_hdr();

//    used in classification
//    std::string TR = bcf_hdr_append_info_with_backup_naming(odw->hdr, "TR", "1", "String", \"Tandem repeat associated with this indel.", rename);
//    bcf_hdr_append(odw->hdr, "##INFO=<ID=LARGE_REPEAT_REGION,Number=0,Type=Flag,Description=\"Very large repeat region, vt only detects up to 1000bp long regions.\">");
//    bcf_hdr_append(odw->hdr, "##INFO=<ID=FLANKSEQ,Number=1,Type=String,Description=\"Flanking sequence 10bp on either side of detected repeat region.\">");

    /////////////////////////
    //filter initialization//
    /////////////////////////
    filter.parse(fexp.c_str(), false);
    filter_exists = fexp=="" ? false : true;

    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_variants = 0;
    no_indels = 0;
    no_added_vntrs = 0;
    no_duplicate_vntrs = 0;

    ////////////////////////
    //tools initialization//
    ////////////////////////
    refseq = new ReferenceSequence(ref_fasta_file);
}

/**
 * Inserts a VNTR record.
 * Returns true if successful.
 */
void VNTRExtractor::insert(Variant* var)
{
//    std::cerr << "buffer size : " << vbuffer.size() << "\n";

    flush(var);

    Variant& nvar = *var;

    std::list<Variant*>::iterator i =vbuffer.begin();
    while(i != vbuffer.end())
    {
        Variant& cvar = **i;

        if (nvar.rid > cvar.rid)
        {
            vbuffer.insert(i, &nvar);
            return;
        }
        else if (nvar.rid == cvar.rid)
        {
            if (nvar.beg1 > cvar.beg1)
            {
               vbuffer.insert(i, &nvar);
               return;
            }
            else if (nvar.beg1 == cvar.beg1)
            {
                if (nvar.end1 > cvar.end1)
                {
                    vbuffer.insert(i, &nvar);
                    return;
                }
                else if (cvar.end1 == nvar.end1)
                {
                    if (nvar.type==VT_VNTR && cvar.type==VT_VNTR)
                    {
                        //duplicate, do not print
                        if (cvar.vntr.motif == nvar.vntr.motif)
                        {
                            std::string nvar_associated_indel = bcf_get_info_str(nvar.h, nvar.v, "ASSOCIATED_INDEL", "");
                            cvar.vntr.add_associated_indel(nvar_associated_indel);

                            ++no_duplicate_vntrs;
                            bcf_destroy(var->v);
                            delete var;
                            return;
                        }
                        else
                        {
                            vbuffer.insert(i, &nvar);
                            return;
                        }
                    }
                    else
                    {
                         vbuffer.insert(i, &nvar);
                         return;
                    }
                }
                else //nvar.rend1 < cvar.end1
                {
                    ++i;
                }
            }
            else //nvar.beg1 < cvar.beg1
            {
                ++i;
            }
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
            exit(1);
        }
    }

    vbuffer.push_back(var);

//    std::cerr << "exit insert\n";
}

/**
 * Flush variant buffer.
 */
void VNTRExtractor::flush(Variant* var)
{
    if (vbuffer.empty())
    {
        return;
    }

    if (var)
    {
        Variant& nvar = *var;

        //search for point to start deleting from.
        std::list<Variant*>::iterator i = vbuffer.begin();
        while(i!=vbuffer.end())
        {
            Variant& cvar = **i;

            if (nvar.rid > cvar.rid)
            {
                break;
            }
            else if (nvar.rid == cvar.rid)
            {
                if ((cvar.end1+buffer_window_allowance) < nvar.beg1)
                {
                    break;
                }
            }
            else
            {
                fprintf(stderr, "[%s:%d %s] Buffer is unordered\n", __FILE__, __LINE__, __FUNCTION__);
                exit(1);
            }

            ++i;
        }

        while (i!=vbuffer.end())
        {
            process_exit(*i);
            i = vbuffer.erase(i);
        }

    }
    else
    {
        std::list<Variant*>::iterator i = vbuffer.begin();
        while (i!=vbuffer.end())
        {
            var = *i;
            process_exit(*i);
            i = vbuffer.erase(i);
        }
    }
}

/**
 * Flush variant buffer.
 */
void VNTRExtractor::process()
{
    bcf_hdr_t *h = odr->hdr;
    bcf1_t *v = odw->get_bcf1_from_pool();

    while (odr->read(v))
    {
        Variant* var = new Variant(h, v);

        if (filter_exists)
        {
            if (!filter.apply(h, v, var, false))
            {
                delete(var);
                continue;
            }
        }

//        bcf_print(h, v);
        create_and_insert_vntr(*var);
        insert(var);

        if (var->type&VT_INDEL) ++no_indels;
        ++no_variants;

        v = odw->get_bcf1_from_pool();
    }

    flush();
    close();
};

/**
 * Process exiting variant.run
 */
void VNTRExtractor::process_exit(Variant* var)
{
    if (var->type==VT_VNTR)
    {
        std::string indels = var->vntr.get_associated_indels();
        if (indels!="")
        {
            bcf_update_info_string(var->h, var->v, ASSOCIATED_INDEL.c_str(), indels.c_str());
        }

        ++no_added_vntrs;
    }

//    bcf_print(var->h, var->v);
    odw->write(var->v);
}

/**.
 * Close files.
 */
void VNTRExtractor::close()
{
    odw->close();
    odr->close();
}

/**
 * Copies exact VNTR features to finalized features.
 */
void VNTRExtractor::copy_exact_vntr_features_to_final_vntr_features(VNTR& vntr)
{
    vntr.repeat_tract = vntr.exact_repeat_tract;
    vntr.motif =  vntr.exact_motif;
    vntr.basis =  vntr.exact_basis;
    vntr.ru =  vntr.exact_ru;
    vntr.mlen =  vntr.exact_mlen;
    vntr.blen =  vntr.exact_blen;
    vntr.beg1 =  vntr.exact_beg1;
    vntr.end1 =  vntr.exact_end1;
    vntr.comp[0] =  vntr.exact_comp[0];
    vntr.comp[1] =  vntr.exact_comp[1];
    vntr.comp[2] =  vntr.exact_comp[2];
    vntr.comp[3] =  vntr.exact_comp[3];
    vntr.entropy =  vntr.exact_entropy;
    vntr.entropy2 =  vntr.exact_entropy2;
    vntr.kl_divergence =  vntr.exact_kl_divergence;
    vntr.kl_divergence2 =  vntr.exact_kl_divergence2;
    vntr.rl =  vntr.exact_rl;
    vntr.ll =  vntr.exact_ll;
    vntr.no_perfect_ru =  vntr.exact_no_perfect_ru;
    vntr.no_ru =  vntr.exact_no_ru;
    vntr.score =  vntr.exact_score;
    vntr.trf_score =  vntr.exact_trf_score;
}

/**
 * Copies fuzzy VNTR features to finalized features.
 */
void VNTRExtractor::copy_fuzzy_vntr_features_to_final_vntr_features(VNTR& vntr)
{
    vntr.repeat_tract = vntr.fuzzy_repeat_tract;
    vntr.motif =  vntr.fuzzy_motif;
    vntr.basis =  vntr.fuzzy_basis;
    vntr.ru =  vntr.fuzzy_ru;
    vntr.mlen =  vntr.fuzzy_mlen;
    vntr.blen =  vntr.fuzzy_blen;
    vntr.beg1 =  vntr.fuzzy_beg1;
    vntr.end1 =  vntr.fuzzy_end1;
    vntr.comp[0] =  vntr.fuzzy_comp[0];
    vntr.comp[1] =  vntr.fuzzy_comp[1];
    vntr.comp[2] =  vntr.fuzzy_comp[2];
    vntr.comp[3] =  vntr.fuzzy_comp[3];
    vntr.entropy =  vntr.fuzzy_entropy;
    vntr.entropy2 =  vntr.fuzzy_entropy2;
    vntr.kl_divergence =  vntr.fuzzy_kl_divergence;
    vntr.kl_divergence2 =  vntr.fuzzy_kl_divergence2;
    vntr.rl =  vntr.fuzzy_rl;
    vntr.ll =  vntr.fuzzy_ll;
    vntr.no_perfect_ru =  vntr.fuzzy_no_perfect_ru;
    vntr.no_ru =  vntr.fuzzy_no_ru;
    vntr.score =  vntr.fuzzy_score;
    vntr.trf_score =  vntr.fuzzy_trf_score;
}

/**
 * Creates a VNTR record based on classification schema.
 */
void VNTRExtractor::create_and_insert_vntr(Variant& nvar)
{
    if (!(nvar.type&VT_INDEL))
    {
        return;
    }
    
    bool insert_vntr = false;
    VNTR& vntr = nvar.vntr;
    nvar.update_vntr_from_info_fields();

//    bcf_print(nvar.h, nvar.v);

    //convert all indels into their corresponding VNTR representation using exact parameters
    if (vntr_classification==EXACT_VNTR)
    {
        if (vntr.exact_repeat_tract == "")
        {
            refseq->fetch_seq(nvar.chrom, vntr.exact_beg1, vntr.exact_end1, vntr.exact_repeat_tract);
        }

        copy_exact_vntr_features_to_final_vntr_features(vntr);

        insert_vntr = true;
    }
    else if (vntr_classification==FUZZY_VNTR)
    {
        if (vntr.fuzzy_repeat_tract == "")
        {
            refseq->fetch_seq(nvar.chrom, vntr.fuzzy_beg1, vntr.fuzzy_end1, vntr.fuzzy_repeat_tract);
        }

        copy_fuzzy_vntr_features_to_final_vntr_features(vntr);

        insert_vntr = true;
    }
    else if (vntr_classification==TAN_KANG2015)
    {
        if ((vntr.fuzzy_rl - vntr.fuzzy_mlen) >= 6 && vntr.fuzzy_no_perfect_ru>=2)
        {
            if (vntr.fuzzy_mlen==1 && vntr.fuzzy_score>0.9)
            {
                insert_vntr = true;
            }
            else if (vntr.fuzzy_mlen>1 && vntr.fuzzy_score>0.75)
            {
                insert_vntr = true;
            }
        }

        if (vntr.fuzzy_repeat_tract == "")
        {
            refseq->fetch_seq(nvar.chrom, vntr.fuzzy_beg1, vntr.fuzzy_end1, vntr.fuzzy_repeat_tract);
        }

        copy_fuzzy_vntr_features_to_final_vntr_features(vntr);
    }
    else if (vntr_classification==WILLEMS2014)
    {
        if ((vntr.mlen==1 && vntr.rl>=6) ||
                (vntr.mlen==2 && vntr.rl>=11) ||
                (vntr.mlen==3 && vntr.rl>=14) ||
                (vntr.mlen==4 && vntr.rl>=14) ||
                (vntr.mlen==5 && vntr.rl>=16) ||
                (vntr.mlen==6 && vntr.rl>=17) ||
                (vntr.mlen>=7 && vntr.rl>=vntr.mlen*2))
        {
            insert_vntr = true;
        }
    }
    else if (vntr_classification==ANANDA2013)
    {
        if ((vntr.mlen==1 && vntr.rl>=2) ||
                (vntr.mlen==2 && vntr.rl>=4) ||
                (vntr.mlen==3 && vntr.rl>=6) ||
                (vntr.mlen==4 && vntr.rl>=8) ||
                (vntr.mlen==5 && vntr.rl>=10) ||
                (vntr.mlen==6 && vntr.rl>=12) ||
                (vntr.mlen>=7 && vntr.rl>=vntr.mlen*2))
        {
            insert_vntr = true;
        }
    }
    else if (vntr_classification==FONDON2012)
    {
        if ((vntr.mlen==1 && vntr.rl>=6) ||
                (vntr.mlen==2 && vntr.rl>=13) ||
                (vntr.mlen==3 && vntr.rl>=20) ||
                (vntr.mlen==4 && vntr.rl>=23) ||
                (vntr.mlen==5 && vntr.rl>=27) ||
                (vntr.mlen==6 && vntr.rl>=27))
        {
            insert_vntr = true;
        }
    }
    else if (vntr_classification==KELKAR2008)
    {
        if ((vntr.mlen==1 && vntr.rl>=6) ||
                (vntr.mlen==2 && vntr.rl>=10) ||
                (vntr.mlen==3 && vntr.rl>=6) ||
                (vntr.mlen==4 && vntr.rl>=8) ||
                (vntr.mlen==5 && vntr.rl>=10) ||
                (vntr.mlen==6 && vntr.rl>=12) ||
                (vntr.mlen>=7 && vntr.rl>=vntr.mlen*2))
        {
            insert_vntr = true;
        }
    }
    else if (vntr_classification==LAI2003)
    {
        if ((vntr.mlen==1 && vntr.rl>=6) ||
                (vntr.mlen==2 && vntr.rl>=8) ||
                (vntr.mlen==3 && vntr.rl>=12) ||
                (vntr.mlen==4 && vntr.rl>=16) ||
                (vntr.mlen==5 && vntr.rl>=20) ||
                (vntr.mlen==6 && vntr.rl>=24) ||
                (vntr.mlen>=7 && vntr.rl>=vntr.mlen*2))
        {
            insert_vntr = true;
        }
    }
    else
    {
        fprintf(stderr, "[%s:%d %s] VNTR classification code not recognized :  %d\n", __FILE__, __LINE__, __FUNCTION__, vntr_classification);
        exit(1);
    }

    if (insert_vntr)
    {
        bcf_hdr_t* h = nvar.h;
        bcf1_t* nv = bcf_init1();
        bcf_clear(nv);

        bcf_set_rid(nv, nvar.rid);
        bcf_set_pos1(nv, vntr.beg1);
        kstring_t s = {0,0,0};
        kputs(vntr.repeat_tract.c_str(), &s);
        kputc(',', &s);
        kputs("<VNTR>", &s);
        bcf_update_alleles_str(h, nv, s.s);
        if (s.m) free(s.s);

        if (no_samples) bcf_update_genotypes(h, nv, gts, no_samples);

        bcf_update_info_int32(h, nv, END.c_str(), &vntr.end1, 1);
        bcf_update_info_string(h, nv, MOTIF.c_str(), vntr.motif.c_str());
        bcf_update_info_string(h, nv, BASIS.c_str(), vntr.basis.c_str());
        bcf_update_info_string(h, nv, RU.c_str(), vntr.ru.c_str());
        bcf_update_info_int32(h, nv, MLEN.c_str(), &vntr.mlen, 1);
        bcf_update_info_int32(h, nv, BLEN.c_str(), &vntr.blen, 1);
        int32_t repeat_tract[2] = {vntr.beg1, vntr.end1};
        bcf_update_info_int32(h, nv, REPEAT_TRACT.c_str(), &repeat_tract, 2);
        bcf_update_info_int32(h, nv, COMP.c_str(), &vntr.comp, 4);
        bcf_update_info_float(h, nv, ENTROPY.c_str(), &vntr.entropy, 1);
        bcf_update_info_float(h, nv, ENTROPY2.c_str(), &vntr.entropy2, 1);
        bcf_update_info_float(h, nv, KL_DIVERGENCE.c_str(), &vntr.kl_divergence, 1);
        bcf_update_info_float(h, nv, KL_DIVERGENCE2.c_str(), &vntr.kl_divergence2, 1);
        bcf_update_info_int32(h, nv, RL.c_str(), &vntr.rl, 1);
        bcf_update_info_int32(h, nv, LL.c_str(), &vntr.ll, 1);
        int32_t ru_count[2] = {vntr.no_perfect_ru, vntr.no_ru};
        bcf_update_info_int32(h, nv, RU_COUNTS.c_str(), &ru_count, 2);
        bcf_update_info_float(h, nv, SCORE.c_str(), &vntr.score, 1);
        bcf_update_info_int32(h, nv, TRF_SCORE.c_str(), &vntr.trf_score, 1);

        Variant *nvntr = new Variant(h, nv);
//        bcf_print(h, nv);
        std::string indel = bcf_variant2string(nvar.h, nvar.v);
        nvntr->vntr.add_associated_indel(indel);

        insert(nvntr);
    }

}