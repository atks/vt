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

#include "consolidate_vntrs.h"

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

    int32_t no_vntrs;

    int32_t no_overlap_vntrs;
    std::vector<int32_t> overlapping_vntr_hist;

    int32_t no_dropped_vntrs;

    //exactness refers to purity of sequence
    //exact   - concordance is 1
    //inexact - concordance is <1

    //isolate is a relative measure and refers to a VNTR that does not overlap with any other detected VNTR.
    //isolated  - does not overlap with another VNTR
    //clustered - overlaps with other VNTRs
    int32_t no_isolated_exact_vntrs;

    int32_t no_perfect_concordance_isolated_exact_vntrs;
    int32_t no_imperfect_concordance_isolated_exact_vntrs;

    int32_t no_perfect_concordance_isolated_inexact_vntrs;
    int32_t no_imperfect_concordance_isolated_inexact_vntrs;


    int32_t no_isolated_inexact_vntrs;

    int32_t no_isolated_complete_overlap_vntrs;
    int32_t no_isolated_incomplete_overlap_vntrs;
    int32_t no_isolated_partial_overlap_vntrs;
    int32_t no_isolated_no_overlap_vntrs;

    int32_t no_clustered_exact_vntrs;
    int32_t no_clustered_inexact_vntrs;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.57";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Consolidates VNTRs by\n"
                 "              1. removing overlapping VNTRs and leaving behind the most complete VNTR\n"
                 "              2. resolves adjacent VNTRs\n"
                 "              3. marking VNTRs that are \n"
                 "              4. Adds INFO fields indicating the number of variants that overlap with this variant";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            debug = arg_debug.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 3000);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=shorter_vntr,Description=\"Another VNTR overlaps with this VNTR.\">");
        odw->write_hdr();

        overlap_vntr = const_cast<char*>("overlap_vntr");
        overlap_vntr_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr");

        ////////////////////////
        //stats initialization//
        ////////////////////////

        no_total_variants = 0;
        no_vntrs = 0;
        no_overlap_vntrs = 0;
        no_dropped_vntrs = 0;

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
        vm = new VariantManip();
    }

    /**
     * Update distribution of overlapping VNTRs
     */
    void update_overlapping_vntr_hist(int32_t no_overlapping_vntrs)
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
     */
    void insert_variant_record_into_buffer(Variant* variant)
    {
        std::list<Variant *>::iterator i = variant_buffer.begin();

        if (variant->type==VT_VNTR)
        {
            ++no_vntrs;
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
                //after most recent variant, we need to have the 1000 window, because the variants are roughly
                //ordered by start.  It is possible to have the start positions changed when merging VNTRs
                //resulting in unordered variants.
                else if (variant->beg1 > cvariant->end1 + 1000)
                {
                    break;
                }
                else if (variant->end1 >= cvariant->beg1 && variant->beg1 <= cvariant->end1) //overlaps
                {
                    if (variant->type==VT_VNTR && cvariant->type==VT_VNTR)
                    {
                        bcf1_t* v = variant->v;
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
                if (variant->beg1 < beg1-1000)
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
     * Compute purity by sequence content.
     */
    float compute_purity_by_sequence_content(char* repeat_tract, char* motif)
    {
        uint32_t motif_count[20];
        motif_count[0] = 0;
        motif_count[2] = 0;
        motif_count[6] = 0;
        motif_count[19] = 0;

        //count bases
        char* motif_ptr = motif;
        while (*motif_ptr)
        {
            ++motif_count[*motif_ptr-'A'];
            ++motif_ptr;
        }

        uint32_t dmb = 0;
        if (motif_count[0]) ++dmb;
        if (motif_count[2]) ++dmb;
        if (motif_count[6]) ++dmb;
        if (motif_count[19]) ++dmb;

        uint32_t repeat_tract_count[20];
        repeat_tract_count[0] = 0;
        repeat_tract_count[2] = 0;
        repeat_tract_count[6] = 0;
        repeat_tract_count[19] = 0;

        uint32_t len = 0;
        char* repeat_tract_ptr = repeat_tract;
        while (*repeat_tract_ptr)
        {
            ++repeat_tract_count[*repeat_tract_ptr-'A'];
            ++len;
            ++repeat_tract_ptr;
        }

        uint32_t db = 0;

        db += motif_count[0] ? repeat_tract_count[0] : 0;
        db += motif_count[2] ? repeat_tract_count[2] : 0;
        db += motif_count[6] ? repeat_tract_count[6] : 0;
        db += motif_count[19] ? repeat_tract_count[19] : 0;

        return (float) db / (float) len;
    }

    /**
     * Consolidate multiallelic variant based on associated biallelic records
     * stored in vs.  Updates v which is to be the consolidated multiallelic
     * variant.
     *
     * returns true if the multiallelic variant is good to go.
     */
    bool consolidate_multiple_overlapping_vntrs(Variant* variant)
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

                bcf_print(odw->hdr, variant->v);
            }

            bcf1_t* vntr_v = variant->v;

            char* motif = NULL;
            int32_t n_motif = 0;
            float* fuzzy_concordance = NULL;
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
//                std::cerr << "1" << ") " << motif << "\t" << fuzzy_concordance[0] << "\t" << fuzzy_flanks[0] << "," << fuzzy_flanks[1] << "\n";
//                std::cerr << "\t" << bcf_get_ref(variant->v) << "\n";

                float impurity = compute_purity_by_sequence_content(bcf_get_ref(variant->v), motif);

                if (flanks[0]==fuzzy_flanks[0]  &&
                    flanks[1]==fuzzy_flanks[1])
                {
                    ++no_isolated_exact_vntrs;

                    if (fuzzy_concordance[0]==1)
                    {
                        ++no_perfect_concordance_isolated_exact_vntrs;
                    }
                    else
                    {
                        std::cerr  << "################\n";
                        std::cerr  << "#1 isolated VNTR\n";
                        std::cerr  << "################\n";
                        std::cerr << "no overlapping SNPs   " << variant->no_overlapping_snps << "\n";
                        std::cerr << "no overlapping Indels " << variant->no_overlapping_indels << "\n";
                        std::cerr << "no overlapping VNTRs  " << variant->no_overlapping_vntrs << "\n";
                        std::cerr << "consolidating: " << variant->vs.size() << " alleles\n";

                        bcf_print(odw->hdr, variant->v);

//                        std::cerr << "1" << ") " << motif << "\t" << fuzzy_concordance[0] << "\t" << flanks[0] << "," << flanks[1] << "\t" << fuzzy_flanks[0] << "," << fuzzy_flanks[1] << "\n";
//                        std::cerr << "\t" << bcf_get_ref(variant->v) << "\n";

                        ///////////////////////////////////////////////////////////////
                        //large deletions OR repeat tract contains inexact repeat units
                        ///////////////////////////////////////////////////////////////

//TTTG and TTA sandwiches a perfect 4 copies of TTGTTGTTGTTG
//20    48231646    .   TTTGTTGTTGTTGTTGTTA T   .   .   NSAMPLES=1;E=10;N=16;ESUM=10;NSUM=16;FLANKS=48231646,48231678;FZ_FLANKS=48231646,48231678;FLANKSEQ=TTTGATTGGT[TTGTTGTTGTTGTTGTTATTGTTGTTGTTGT]CGTCATTGTT;GMOTIF=GTT;TR=20:48231647:TTGTTGTTGTTGTTGTTATTGTTGTTGTTGT:<VNTR>:GTT
//20    48231647    .   TTGTTGTTGTTGTTGTTATTGTTGTTGTTGT <VNTR>  .   .   MOTIF=GTT;RU=TTG;FUZZY;FZ_CONCORDANCE=0.969697;FZ_RL=31;FZ_LL=0;FLANKS=48231646,48231678;FZ_FLANKS=48231646,48231678;FZ_RU_COUNTS=10,11;FLANKSEQ=TTTGATTGGT[TTGTTGTTGTTGTTGTTATTGTTGTTGTTGT]CGTCATTGTT

//2 copies of TTTTAG sandwiches a TTAAC.  i.e. TTTTAG[TTAAC]TTTTAG
//20    10546879    .   GTTTTAGTTAAC    G   .   .   NSAMPLES=1;E=21;N=48;ESUM=21;NSUM=48;FLANKS=10546879,10546897;FZ_FLANKS=10546879,10546897;FLANKSEQ=ATTGCCATTG[TTTTAGTTAACTTTTAG]CACTGGGTAT;GMOTIF=AGTTTT;TR=20:10546880:TTTTAGTTAACTTTTAG:<VNTR>:AGTTTT
//20    10546880    .   TTTTAGTTAACTTTTAG   <VNTR>  .   .   MOTIF=AGTTTT;RU=TTTTAG;FUZZY;FZ_CONCORDANCE=0.833333;FZ_RL=17;FZ_LL=0;FLANKS=10546879,10546897;FZ_FLANKS=10546879,10546897;FZ_RU_COUNTS=2,3;FLANKSEQ=ATTGCCATTG[TTTTAGTTAACTTTTAG]CACTGGGTAT

                        ++no_imperfect_concordance_isolated_exact_vntrs;
                    }

//                    odw->write(variant->v);
//                    variant->v = NULL;
//                    delete variant;
//                    variant_buffer.pop_back();
                }
                else
                {
                    //complete overlap
                    //most should have imperfect concordance
                    //those that have perfect concordance implies that the alternate allele resulted in a imperfect VNTR
                    if (flanks[0]>=fuzzy_flanks[0]  &&
                        flanks[1]<=fuzzy_flanks[1])
                    {
//                        std::cerr << "1" << ") " << motif << "\t" << fuzzy_concordance[0] << "\t" << flanks[0] << "," << flanks[1] << "\t" << fuzzy_flanks[0] << "," << fuzzy_flanks[1] << "\n";
//                        std::cerr << "\t" << impurity << "\t" << bcf_get_ref(variant->v) << "\n";
                       ++no_isolated_complete_overlap_vntrs;
                    }
                    //partial overlaps
                    //these are induced possibly by errors at the boundary of VNTRs
                    //
                    //
                    else if (flanks[0]<=fuzzy_flanks[1]  &&
                             flanks[1]>=fuzzy_flanks[0])
                    {
//                      std::cerr << "1" << ") " << motif << "\t" << fuzzy_concordance[0] << "\t" << flanks[0] << "," << flanks[1] << "\t" << fuzzy_flanks[0] << "," << fuzzy_flanks[1] << "\n";
//                      std::cerr << "\t" << impurity << "\t" << bcf_get_ref(variant->v) << "\n";
                        ++no_isolated_partial_overlap_vntrs;
                    }
                    else
                    {
//                      std::cerr << "1" << ") " << motif << "\t" << fuzzy_concordance[0] << "\t" << flanks[0] << "," << flanks[1] << "\t" << fuzzy_flanks[0] << "," << fuzzy_flanks[1] << "\n";
//                      std::cerr << "\t" << impurity << "\t" << bcf_get_ref(variant->v) << "\n";
                        ++no_isolated_no_overlap_vntrs;
                    }

                    if (fuzzy_concordance[0]==1)
                    {
                        ++no_perfect_concordance_isolated_inexact_vntrs;
                    }
                    else
                    {
                        ++no_imperfect_concordance_isolated_inexact_vntrs;
                    }

                    ++no_perfect_concordance_isolated_inexact_vntrs;

                    ++no_isolated_inexact_vntrs;
                }

                free(motif);
                free(fuzzy_concordance);
                free(flanks);
                free(fuzzy_flanks);
            }
        }
        else if (variant->no_overlapping_vntrs >= 1)
        {
            if ((true && variant->vntr_vs.size()>6) || debug)
            {
                std::cerr  << "###################################\n";
                std::cerr  << "#2 or more VNTR and multiple Indels\n";
                std::cerr  << "###################################\n";
                std::cerr << "no overlapping SNPs   " << variant->no_overlapping_snps << "\n";
                std::cerr << "no overlapping Indels " << variant->no_overlapping_indels << "\n";
                std::cerr << "no overlapping VNTRs  " << variant->no_overlapping_vntrs << "\n";
                std::cerr << "consolidating: " << (variant->vs.size()+1) << " alleles\n";

                for (uint32_t i=0; i<variant->vntr_vs.size(); ++i)
                {
                    bcf1_t* vntr_v = variant->vntr_vs[i];

                    char* motif = NULL;
                    int32_t n_motif = 0;
                    float* fuzzy_concordance = NULL;
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
                        float impurity = compute_purity_by_sequence_content(bcf_get_ref(vntr_v), motif);

                        std::cerr << (i+1) << ") " << motif << "\t" << fuzzy_concordance[0] << "\t" << flanks[0] << "," << flanks[1] << "\t" << fuzzy_flanks[0] << "," << fuzzy_flanks[1] << "\n";
                        std::cerr << "\t" << bcf_get_ref(vntr_v) << "\n";
                        bcf_print(odw->hdr, vntr_v);
                        
                        free(motif);
                        free(fuzzy_concordance);
                        free(flanks);
                        free(fuzzy_flanks);


                    }

//                    bcf_print(odw->hdr, variant->vntr_vs[i]);
                }
            }


            //if all motifs are consistent, take the one with the largest region
            
//            bool is_motif_consistent




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

    void consolidate_vntrs()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();

        Variant* variant;
        while (odr->read(v))
        {
            variant = new Variant(odw->hdr, v);

            flush_variant_buffer(variant);
            insert_variant_record_into_buffer(variant);
            v = odw->get_bcf1_from_pool();

            ++no_total_variants;
        }

        flush_variant_buffer();

        //you have to close the writer first as the header used by the writer is linked to the reader
        odw->close();
        odr->close();
    };

    void print_options()
    {
        std::clog << "consolidate_vntrs v" << version << "\n\n";

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
        std::clog << "       VNTR Classification          \n";
        std::clog << "       Number of VNTRs      " << no_vntrs << "\n";
        std::clog << "\n";
        std::clog << "       Number of isolated exact VNTRs      " << no_isolated_exact_vntrs << "\n";
        std::clog << "                     perfect concordance   " << no_perfect_concordance_isolated_exact_vntrs << "\n";
        std::clog << "                     imperfect concordance " << no_imperfect_concordance_isolated_exact_vntrs << "\n";
        std::clog << "       Number of isolated inexact VNTRs    " << no_isolated_inexact_vntrs << "\n";
        std::clog << "                     complete overlaps     " << no_isolated_complete_overlap_vntrs << "\n";
        std::clog << "                     partial overlaps      " << no_isolated_partial_overlap_vntrs << "\n";
        std::clog << "                     no overlaps           " << no_isolated_no_overlap_vntrs << "\n";
        std::clog << "       Number of clustered exact VNTRs     " << no_clustered_exact_vntrs << "\n";

        for (uint32_t i=0; i<overlapping_vntr_hist.size(); ++i)
        {
            if (overlapping_vntr_hist[i])
            {
                std::clog << "       " << i << "    "  << overlapping_vntr_hist[i] << "\n";
            }
        }


        std::clog << "       Number of clustered inexact VNTRs   " << no_clustered_inexact_vntrs << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void consolidate_vntrs(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.consolidate_vntrs();
    igor.print_stats();
};
