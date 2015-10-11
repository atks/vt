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

#include "consolidate.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::string ref_fasta_file;
    std::vector<GenomeInterval> intervals;
    bool debug;
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    ////////////////
    //variant buffer
    ////////////////
    std::list<Variant *> variant_buffer; //front is most recent

    ////////////
    //filter ids
    ////////////
    char* overlap_snp;
    char* overlap_indel;
    char* overlap_vntr;
    int32_t overlap_snp_id;
    int32_t overlap_indel_id;
    int32_t overlap_vntr_id;

    /////////
    //stats//
    /////////
    int32_t no_total_variants;
    int32_t no_nonoverlap_variants;
    int32_t no_overlap_variants;
    int32_t no_new_multiallelic_snps;
    int32_t no_new_multiallelic_indels;
    int32_t no_new_multiallelic_vntr_indels;
    int32_t no_overlap_vntrs;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    faidx_t * fai;
    LogTool lt;

    Igor(int argc, char **argv)
    {
        version = "0.57";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Consolidates variants in the following ways\n"
                 "              1.Updates each variant's filter field if it overlaps with another variant.\n"
                 "                    SNP   - overlap at position \n"
                 "                    Indel - overlap with region bounded by exact flanks\n"
                 "                    VNTR  - overlap with region bounded by fuzzy flanks\n"
                 "              2.Adds candidate multiallelic SNPs that do not overlap with Indels or VNTRs\n"
                 "              3.Adds candidate multiallelic Indels if and only if\n"
                 "                    a. they do not overlap with SNPs and VNTRs\n"
                 "                    b. the exact flanks and fuzzy flanks are equal to one another\n"
                 "              4.Adds INFO fields indicating the number of variants that overlap with this variant";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
            debug = arg_debug.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    /**
     * Compute GLFSingle single variant likelihood ratio P(Non Variant)/P(Variant) for Indel
     */
    double compute_glfsingle_llr(uint32_t e, uint32_t n)
    {
        double ln_theta = -6.907755; // theta = 0.001;
        double ln_one_minus_theta = -0.0010005; // 1-theta = 0.999;
        double ln_one_third = -1.098612; 
        double ln_two_thirds = -0.4054651; 
        double ln_0_001 = -6.907755;                    
        double ln_0_999 = -0.0010005; 
        double ln_0_5 = -0.6931472;         
        
        double ln_pRR = (n-e) * ln_0_999 + e * ln_0_001;
        double ln_pRA = n * ln_0_5;
        double ln_pAA = e * ln_0_999 + (n-e) * ln_0_001;
        
        double ln_lr = ln_one_minus_theta + ln_pRR;
        ln_lr = logspace_add(ln_lr, ln_one_third+ln_theta+ln_pRA);
        ln_lr = logspace_add(ln_lr, ln_two_thirds+ln_theta+ln_pAA);
        
//        std::cerr << "\t\t" << ln_pRR << " " <<  ln_lr <<  " ";
        
        ln_lr = ln_pRR - ln_lr;
        
//         std::cerr <<  ln_lr << "\n";
        
        return ln_lr;
    }
    
    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_snp,Description=\"Overlaps with SNP.\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_indel,Description=\"Overlaps with Indel.\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_vntr,Description=\"Overlaps with VNTR.\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=shorter_vntr,Description=\"Another VNTR overlaps with this VNTR.\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=on_vntr_boundary,Description=\"This variant lies near a VNTR boundary.\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=fail_lr,Description=\"Fail likelihood ratio cutoff.\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=OVERLAPS,Number=3,Type=Integer,Description=\"Number of SNPs, Indels and VNTRs overlapping with this variant.\">");
        odw->write_hdr();      

        overlap_snp = const_cast<char*>("overlap_snp");
        overlap_indel = const_cast<char*>("overlap_indel");
        overlap_vntr = const_cast<char*>("overlap_vntr");

        overlap_snp_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_snp");
        overlap_indel_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_indel");
        overlap_vntr_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr");

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_total_variants = 0;
        no_nonoverlap_variants = 0;
        no_overlap_variants = 0;
        no_new_multiallelic_snps = 0;
        no_new_multiallelic_indels = 0;
        no_new_multiallelic_vntr_indels = 0;
        no_overlap_vntrs = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip();
        fai = fai_load(ref_fasta_file.c_str());
    }

    /**
     * Inserts a Variant record.
     */
    void insert_variant_record_into_buffer(Variant* variant)
    {
        char* tr = NULL;
        int32_t n = 0;

//        //this is for indels that are annotated with TR - you cannot
//        //detect by overlap of the reference sequence because in the normalized
//        //representation for insertions, only the anchor base is present.
//        //This is not an issue for deletions.
//        //reminder: consolidate is for a VCF file that is produced by annotate_indels.
//        if (bcf_get_info_string(odw->hdr, variant->v, "TR", &tr, &n)>0)
//        {
//            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
//            free(tr);
//        }

        //update filters
        std::list<Variant *>::iterator i = variant_buffer.begin();
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
                else if (variant->beg1 > cvariant->end1+1000) //does not overlap
                {
                    break;
                }
                else if (variant->end1 >= cvariant->beg1 && variant->beg1 <= cvariant->end1) //overlaps
                {
                    ///////
                    //SNP//
                    ///////
                    if (variant->type==VT_SNP)
                    {
                        if (cvariant->type==VT_SNP)
                        {
                            //make a new candidate multiallelic SNP 
                            //only if there are overlapping SNPs
                            if (bcf_get_n_filter(cvariant->v)==0)
                            {
                                Variant* multiallelic_variant = new Variant(cvariant, variant);
                                multiallelic_variant->no_overlapping_snps = 2;
                                variant_buffer.push_front(multiallelic_variant);
                            }

                            bcf_add_filter(odw->hdr, variant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_snp_id);
                            ++cvariant->no_overlapping_snps;
                        }
                        else if (cvariant->type==VT_INDEL)
                        {
                            bcf_add_filter(odw->hdr, variant->v, overlap_indel_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_snp_id);
                            ++cvariant->no_overlapping_snps;
                        }
                        else if (cvariant->type==VT_VNTR)
                        {
                            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
                            ++variant->no_overlapping_vntrs;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                        }
                        else if (cvariant->type==VT_UNDEFINED)
                        {
                            //bcf1_t* v = bcf_dup(variant->v);
                            bcf1_t* v = variant->v;
                            cvariant->vs.push_back(v);
                            cvariant->snp_vs.push_back(v);
                            ++cvariant->no_overlapping_snps;
                        }
                    }
                    /////////
                    //INDEL//
                    /////////
                    else if (variant->type==VT_INDEL)
                    {
                        if (cvariant->type==VT_SNP)
                        {
                            bcf_add_filter(odw->hdr, variant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_indel_id);
                            ++cvariant->no_overlapping_indels;
                        }
                        else if (cvariant->type==VT_INDEL)
                        {
                            //make a new a multiallelic only if 2 indels overlap with one another
                            if (bcf_get_n_filter(cvariant->v) == 0)
                            {
                                Variant* multiallelic_variant = new Variant(cvariant, variant);
                                multiallelic_variant->no_overlapping_indels = 2;
                                variant_buffer.push_front(multiallelic_variant);
                            }

                            bcf_add_filter(odw->hdr, variant->v, overlap_indel_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_indel_id);
                            ++cvariant->no_overlapping_indels;
                        }
                        else if (cvariant->type==VT_VNTR)
                        {
                            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
                            ++variant->no_overlapping_vntrs;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_indel_id);
                            ++cvariant->no_overlapping_indels;
                        }
                        else if (cvariant->type==VT_UNDEFINED)
                        {
                            //bcf1_t* v = bcf_dup(variant->v);
                            bcf1_t* v = variant->v;
                            cvariant->vs.push_back(v);
                            cvariant->indel_vs.push_back(v);
                            ++cvariant->no_overlapping_indels;
                        }
                    }
                    ////////
                    //VNTR//
                    ////////
                    else if (variant->type==VT_VNTR)
                    {
                        if (cvariant->type==VT_SNP)
                        {
                            bcf_add_filter(odw->hdr, variant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_vntr_id);
                            ++cvariant->no_overlapping_vntrs;
                        }
                        else if (cvariant->type==VT_INDEL)
                        {
                            bcf_add_filter(odw->hdr, variant->v, overlap_indel_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_vntr_id);
                            ++cvariant->no_overlapping_vntrs;
                        }
                        else if (cvariant->type==VT_VNTR)
                        {
                            //can mark the overlapping VNTRs
                            //write over
                            //this should perhaps only be handled in vntr overlaps with indels?
                            //this is a bit tricky.
                            //consider ONLY overlapping VNTRs
                            
                            
                            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_vntr_id);
                            ++cvariant->no_overlapping_vntrs;

                            //mark a better VNTR?? by score and by length?
                        }
                        else if (cvariant->type==VT_UNDEFINED)
                        {
                            //bcf1_t* v = bcf_dup(variant->v);
                            bcf1_t* v = variant->v;
                            cvariant->vs.push_back(v);
                            cvariant->vntr_vs.push_back(v);
                            ++cvariant->no_overlapping_vntrs;
                            
                            ++no_overlap_vntrs;
                        }
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

        variant_buffer.push_front(variant);
    }

    /**
     * Flush variant buffer.
     */
    void flush_variant_buffer(Variant* var)
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
                if (variant->type==VT_UNDEFINED)
                {
                    std::cerr << "PRINT CONSOLIDATED VARIANT\n";

                    if (consolidate_multiallelic(variant))
                    {
                        int32_t overlaps[3] = {variant->no_overlapping_snps, variant->no_overlapping_indels, variant->no_overlapping_vntrs};
                        bcf_update_info_int32(odw->hdr, variant->v, "OVERLAPS", &overlaps, 3);
                        odw->write(variant->v);
                        variant->v = NULL;
                        delete variant;
                        variant_buffer.pop_back();
                    }
                }
                else
                {
                    int32_t overlaps[3] = {variant->no_overlapping_snps, variant->no_overlapping_indels, variant->no_overlapping_vntrs};
                    bcf_update_info_int32(odw->hdr, variant->v, "OVERLAPS", &overlaps, 3);
                    odw->write(variant->v);
                    variant->v = NULL;
                    delete variant;
                    variant_buffer.pop_back();
                }
            }
            else if (variant->rid == rid)
            {
                if (variant->beg1 < beg1-1000)
                {
                    if (variant->type==VT_UNDEFINED)
                    {
                        if (consolidate_multiallelic(variant))
                        {
                            int32_t overlaps[3] = {variant->no_overlapping_snps, variant->no_overlapping_indels, variant->no_overlapping_vntrs};
                            bcf_update_info_int32(odw->hdr, variant->v, "OVERLAPS", &overlaps, 3);
                            odw->write(variant->v);
                            variant->v = NULL;
                        }

                        delete variant;
                        variant_buffer.pop_back();
                    }
                    else
                    {
                        int32_t overlaps[3] = {variant->no_overlapping_snps, variant->no_overlapping_indels, variant->no_overlapping_vntrs};
                        bcf_update_info_int32(odw->hdr, variant->v, "OVERLAPS", &overlaps, 3);
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
    bool consolidate_multiallelic(Variant* variant)
    {
        if (variant->no_overlapping_snps !=0 &&
            variant->no_overlapping_indels ==0 &&
            variant->no_overlapping_vntrs ==0)
        {
            bcf1_t* v = bcf_init1();
            bcf_clear(v);
            std::vector<bcf1_t*> vs = variant->vs;
            
            bcf_set_rid(v, bcf_get_rid(vs[0]));
            bcf_set_pos1(v, bcf_get_pos1(vs[0]));

            kstring_t s = {0,0,0};
            kputc(bcf_get_snp_ref(vs[0]), &s);
            kputc(',', &s);
            int32_t no_alleles = vs.size();


            char alts[no_alleles];
            for (uint32_t i=0; i<no_alleles; ++i)
            {
                bcf_unpack(vs[i], BCF_UN_STR);
                alts[i] = bcf_get_snp_alt(vs[i]);
            }
            //selection sort
            for (uint32_t i=0; i<no_alleles-1; ++i)
            {
                for (uint32_t j=i+1; j<no_alleles; ++j)
                {
                    if (alts[j]<alts[i])
                    {
                        char tmp = alts[j];
                        alts[j] = alts[i];
                        alts[i] = tmp;
                    }
                    alts[i] = bcf_get_snp_alt(vs[i]);
                }

                kputc(alts[i], &s);
                kputc(',', &s);
            }
            kputc(alts[no_alleles-1], &s);
            bcf_update_alleles_str(odw->hdr, v, s.s);

            variant->v = v;

            ++no_new_multiallelic_snps;
            return true;
        }
        else 
        {
            //////////////////////
            //One overlapping VNTR
            //////////////////////
            if (variant->no_overlapping_vntrs==1)
            {
                if (debug)
                {
                    std::cerr  << "###########################\n";
                    std::cerr  << "#1 VNTR and multiple Indels\n";
                    std::cerr  << "###########################\n";
                    std::cerr << "no overlapping SNPs   " << variant->no_overlapping_snps << "\n";
                    std::cerr << "no overlapping Indels " << variant->no_overlapping_indels << "\n";
                    std::cerr << "no overlapping VNTRs  " << variant->no_overlapping_vntrs << "\n";
                    std::cerr << "consolidating: " << (variant->vs.size()+1) << " alleles\n";
                }
                
                bcf1_t* vntr_v = variant->vntr_vs[0];
                if (debug)
                {    
                    bcf_print(odw->hdr, vntr_v);
                    std::cerr << "\tQUAL = " << bcf_get_qual(vntr_v) << "\n";
                }
                
                char* motif = NULL;
                int32_t n_motif = 0;
                double* fuzzy_concordance = NULL;
                int32_t n_fuzzy_concordance = 0;
                int32_t* flanks = NULL;
                int32_t n_flanks = 0;
                int32_t* fuzzy_flanks = NULL;
                int32_t n_fuzzy_flanks = 0;
                if (bcf_get_info_string(odr->hdr, vntr_v, "MOTIF", &motif, &n_motif)>0 &&
                    bcf_get_info_float(odr->hdr, vntr_v, "FZ_CONCORDANCE", &fuzzy_concordance, &n_fuzzy_concordance)>0 &&
                    bcf_get_info_int32(odr->hdr, vntr_v, "FLANKS", &flanks, &n_flanks)>0 &&
                    bcf_get_info_int32(odr->hdr, vntr_v, "FZ_FLANKS", &fuzzy_flanks, &n_fuzzy_flanks)>0)
                    
                {
                    std::vector<bcf1_t *>& indel_vs = variant->indel_vs;
                    std::vector<uint32_t> ref_beg1;
                    std::vector<uint32_t> ref_end1;              
                    uint32_t min_beg1 = INT_MAX;
                    uint32_t max_end1 = 0;
                    std::vector<bcf1_t*> candidate_indel_vs;
                        
                    for (uint32_t i=0; i<indel_vs.size(); ++i)
                    {
                        if (debug)
                        {    
                            bcf_print(odw->hdr, indel_vs[i]);
                            std::cerr << "\tQUAL = " << bcf_get_qual(indel_vs[i]) << "\n";
                        }
                        
                        if (bcf_get_qual(indel_vs[i])>30)
                        {
                            if (bcf_get_pos1(indel_vs[i])==flanks[0])
                            {
                                char* gmotif = NULL;
                                int32_t n_gmotif = 0;
                                if (bcf_get_info_string(odr->hdr, indel_vs[i], "GMOTIF", &gmotif, &n_gmotif)>0)
                                {
                                    if (strcmp(motif, gmotif)==0)
                                    {
                                        uint32_t beg1 = bcf_get_pos1(indel_vs[i]);
                                        uint32_t end1 = bcf_get_end1(indel_vs[i]);
                                        ref_beg1.push_back(beg1);
                                        ref_end1.push_back(end1);
                                        min_beg1 = beg1<min_beg1 ? beg1 : min_beg1;
                                        max_end1 = end1>max_end1 ? end1 : max_end1;
                                        candidate_indel_vs.push_back(indel_vs[i]);
                                    }    
                                
                                    free(gmotif);    
                                }
                            }
                            else
                            {
                                //not on correct position
                            }
                        }  
                    }
                    
                    for (uint32_t i=0; i<variant->snp_vs.size(); ++i)
                    {
                        if (debug)
                        {
                            bcf_print(odw->hdr, variant->snp_vs[i]);
                            std::cerr << "\tQUAL = " << bcf_get_qual(variant->snp_vs[i]) << "\n";
                        }
                    }
                    
                    
                    if (candidate_indel_vs.size()>1)
                    {
                        bcf1_t* v = bcf_init1();
                        bcf_clear(v);
                        bcf_set_rid(v, bcf_get_rid(candidate_indel_vs[0]));
                        bcf_set_pos1(v, min_beg1);
                        variant->v = v;
                        
                        kstring_t alleles = {0,0,0};
                        
                        int32_t seq_len = 0;
                        char* seq = faidx_fetch_seq(fai, variant->chrom.c_str(), min_beg1-1, max_end1-1, &seq_len);
                        kputs(seq, &alleles);    
                        
                        for (uint32_t i=0; i<candidate_indel_vs.size(); ++i)
                        {
                            kputc(',', &alleles);
                            std::string alt;
                            alt.append(seq, ref_beg1[i]-min_beg1);
                            alt.append(bcf_get_alt(candidate_indel_vs[i], 1));
                            alt.append(seq, ref_end1[i]-min_beg1+1, max_end1-ref_end1[i]);
                            kputs(alt.c_str(), &alleles);
                            
                        }
                        
                        bcf_update_alleles_str(odw->hdr, v, alleles.s);
                        
                        if (seq_len) free(seq);
                        if (alleles.m) free(alleles.s);
                        
                        if (debug)
                        {    
                            bcf_print(odw->hdr, v);   
                        }
                        
                        ++no_new_multiallelic_vntr_indels;
                        ++no_new_multiallelic_indels;
                        
                        free(motif);
                        free(fuzzy_concordance);
                        free(flanks);
                        free(fuzzy_flanks); 
                        
                        return true;     
                    }    
                    else
                    {
                        free(motif);
                        free(fuzzy_concordance);
                        free(flanks);
                        free(fuzzy_flanks);
                        
                        return false;
                    }                   
                } 
            }  
            //////////////////////////////
            //Two or more overlapping VNTR
            //////////////////////////////
            else if (variant->no_overlapping_vntrs > 1)
            {
                if (true)
                {
                    std::cerr  << "###################################\n";
                    std::cerr  << "#2 or more VNTR and multiple Indels\n";
                    std::cerr  << "###################################\n";
                    std::cerr << "no overlapping SNPs   " << variant->no_overlapping_snps << "\n";
                    std::cerr << "no overlapping Indels " << variant->no_overlapping_indels << "\n";
                    std::cerr << "no overlapping VNTRs  " << variant->no_overlapping_vntrs << "\n";
                    std::cerr << "consolidating: " << (variant->vs.size()+1) << " alleles\n";
                }
                
                for (uint32_t i=0; i<variant->vs.size(); ++i)
                {
                    if (true)
                    {    
                        bcf_print(odw->hdr, variant->vs[i]);
                        std::cerr << "\tQUAL = " << bcf_get_qual(variant->vs[i]) << "\n";
                    }
                    
                   
                }
                
                
                return false;
            }
            else //no VNTRs
            {
                if (debug)
                {
                    std::cerr  << "###########\n";
                    std::cerr  << "OTHER TYPES\n";
                    std::cerr  << "###########\n";
                }
                    
                bcf1_t* v = bcf_init1();
                bcf_clear(v);
    
                std::vector<bcf1_t*>& vs = variant->vs;
    
    //            std::cerr << "no overlapping SNPs " << variant->no_overlapping_snps << "\n";
    //            std::cerr << "consolidating: " << (vs.size()+1) << " alleles\n";
    
    
                bcf_set_rid(v, bcf_get_rid(vs[0]));
                bcf_set_pos1(v, bcf_get_pos1(vs[0]));
    
                kstring_t s = {0,0,0};
                kputc(bcf_get_snp_ref(vs[0]), &s);
                kputc(',', &s);
                int32_t no_alleles = vs.size();
    
    
                char alts[no_alleles];
                for (uint32_t i=0; i<no_alleles; ++i)
                {
                    if (debug)
                    {
                        bcf_unpack(vs[i], BCF_UN_STR);
                        bcf_print(odw->hdr, vs[i]);
                        
                        std::cerr << "\tQUAL = " << bcf_get_qual(vs[i]) << "\n";
                    }
                }
                
                bcf_destroy(v);
            }
           
            
            return false;
        }

        return false;
    }

    /**
     * Flush variant buffer.
     */
    void flush_variant_buffer()
    {
        while (!variant_buffer.empty())
        {
            Variant* variant = variant_buffer.back();

            if (variant->type==VT_UNDEFINED)
            {
                if (consolidate_multiallelic(variant))
                {
                    int32_t overlaps[3] = {variant->no_overlapping_snps, variant->no_overlapping_indels, variant->no_overlapping_vntrs};
                    bcf_update_info_int32(odw->hdr, variant->v, "OVERLAPS", &overlaps, 3);
                    odw->write(variant->v);
                    variant->v = NULL;
                }

                delete variant;
                variant_buffer.pop_back();


            }
            else
            {
                int32_t overlaps[3] = {variant->no_overlapping_snps, variant->no_overlapping_indels, variant->no_overlapping_vntrs};
                bcf_update_info_int32(odw->hdr, variant->v, "OVERLAPS", &overlaps, 3);
                odw->write(variant->v);
                variant->v = NULL;
                delete variant;
                variant_buffer.pop_back();
            }
        }
    }

    void consolidate()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();

        Variant* variant;
        while (odr->read(v))
        {
            variant = new Variant(odw->hdr, v);
            flush_variant_buffer(variant);

            vm->classify_variant(odw->hdr, v, *variant);

            //compute likelihood ratio and store in QUAL field
            if (variant->type&VT_INDEL)
            {
                double ln_lr = 0;
                double max_ln_lr = 0;
                int32_t* e = NULL;
                int32_t n_e = 0;
                int32_t* n = NULL;
                int32_t n_n = 0;
                
                if (bcf_get_info_int32(odr->hdr, v, "E", &e, &n_e)>0 &&
                    bcf_get_info_int32(odr->hdr, v, "N", &n, &n_n)>0)
                {
                    for (uint32_t i=0; i<n_e; ++i)
                    {
                        ln_lr  = compute_glfsingle_llr(e[i], n[i]);
                        ln_lr  = ln_lr>0 ? 0 : -10*(ln_lr-M_LOG10E);
                        max_ln_lr = ln_lr > max_ln_lr ? ln_lr : max_ln_lr;
                    }
                    
                    
                    bcf_set_qual(v, max_ln_lr);
                               
                    free(e);
                    free(n);
                }    
            }    

            insert_variant_record_into_buffer(variant);

            v = odw->get_bcf1_from_pool();

            ++no_total_variants;
        }

        flush_variant_buffer();

        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "consolidate v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        if (intervals.size()!=0)
        {
            std::clog << "         [i] intervals             " << intervals.size() <<  " intervals\n";
        }
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: Total number of observed variants        " << no_total_variants << "\n";
        std::clog << "       Total number of nonoverlap variants      " << no_nonoverlap_variants << "\n";
        std::clog << "       Total number of new multiallelic SNPs    " << no_new_multiallelic_snps << "\n";
        std::clog << "       Total number of new multiallelic Indels  " << no_new_multiallelic_indels << "\n";
        std::clog << "            VNTR                                     " << no_new_multiallelic_vntr_indels << "\n";
        std::clog << "       Total number of overlap VNTRs            " << no_overlap_vntrs << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void consolidate(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.consolidate();
    igor.print_stats();
};
