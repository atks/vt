/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

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

#include "annotate_indels.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string ref_fasta_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    std::string method;           //methods of detection
    std::string annotation_mode;  //modes of annotation
    int32_t vntr_classification;  //vntr classification schemas
    bool override_tag;

    //INFO tags
    std::string TR;               //annotating indels with overlapping tandem repeat

    std::string MOTIF;            //canonical motif of tandem repeat
    std::string RU;               //repeat unit on reference sequence that is in phase from the start position of the variant

    std::string RL;               //repeat tract length
    std::string LL;               //repeat tract length of the longest allele
    std::string CONCORDANCE;      //concordance of the repeat unit
    std::string RU_COUNTS;        //repeat unit counts - exact and inexact
    std::string FLANKS;           //flank positions

    std::string FZ_RL;            //repeat tract length
    std::string FZ_LL;            //repeat tract length of the longest allele
    std::string FZ_CONCORDANCE;   //concordance of the repeat unit
    std::string FZ_RU_COUNTS;     //repeat unit counts - exact and inexact
    std::string FZ_FLANKS;        //flank positions

    //helper variables for populating additional VNTR records
    uint32_t no_samples;
    int32_t* gts;

    //convenience kstring so that we do not have to free memory ALL the time.
    kstring_t s;

    bool debug;

    /////////////
    //vntr buffer
    /////////////
    std::list<VNTR> vntr_buffer; //front is most recent

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    int32_t no_variants_annotated;

    ////////////////
    //common tools//
    ////////////////
    VariantManip* vm;
    VNTRAnnotator* va;
    faidx_t* fai;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "annotates indels with VNTR information - repeat tract length, repeat motif, flank information";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_annotation_mode("a", "a", "annotation type [v]\n"
                 "              v : a. output VNTR variant (defined by classification).\n"
                 "                     RU    repeat unit on reference sequence (CA)\n"
                 "                     MOTIF canonical representation (AC)\n"
                 "                     RL    repeat tract length in bases (11)\n"
                 "                  b. mark indels with overlapping VNTR.\n"
                 "                     TR    position and alleles of VNTR (20:23413:CACACACACAC:<VNTR>)\n"
                 "                     END   end position if allelic region (23423)\n"
                 "              a : annotate each indel with RU, RL, MOTIF, REF.",
                 false, "v", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_vntr_classification("c", "c", "classification schemas of tandem repeat [6]\n"
                 "              1 : lai2003     \n"
                 "              2 : kelkar2008  \n"
                 "              3 : fondon2012  \n"
                 "              4 : ananda2013  \n"
                 "              5 : willems2014 \n"
                 "              6 : tan_kang2015",
                 false, 6, "integer", cmd);
            TCLAP::ValueArg<std::string> arg_method("m", "m", "mode [f]\n"
                 "              e : by exact alignment"
                 "              f : by fuzzy alignment",
                 false, "f", "str", cmd);

            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_override_tag("x", "x", "override tags [false]", cmd, false);
            TCLAP::SwitchArg arg_add_vntr_record("v", "v", "add vntr record [false]", cmd, false);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            method = arg_method.getValue();
            annotation_mode = arg_annotation_mode.getValue();
            vntr_classification = arg_vntr_classification.getValue();
            override_tag = arg_override_tag.getValue();
            fexp = arg_fexp.getValue();
            debug = arg_debug.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
            abort();
        }
    };

    ~Igor() 
    {
        if (s.m) free(s.s);
    };

    void initialize()
    {
        ///////////
        //options//
        ///////////
        if (method!="e" && method!="f")
        {
            fprintf(stderr, "[%s:%d %s] Not a valid mode of VNTR detection: %s\n", __FILE__,__LINE__,__FUNCTION__, method.c_str());
            exit(1);
        }

        if (annotation_mode!="v" && annotation_mode!="a")
        {
            fprintf(stderr, "[%s:%d %s] Not a valid mode of annotation: %s\n", __FILE__,__LINE__,__FUNCTION__, annotation_mode.c_str());
            exit(1);
        }

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 10000);
        odw->link_hdr(odr->hdr);

        //////////////////////////////
        //INFO header adding for VCF//
        //////////////////////////////
        MOTIF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "MOTIF", "1", "String", "Canonical motif in an VNTR or homopolymer", true);
        RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU", "1", "String", "Repeat unit in a VNTR or homopolymer", true);

        RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RL", "1", "Float", "Reference repeat unit length", true);
        LL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "LL", "1", "Float", "Longest repeat unit length", true);
        CONCORDANCE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "CONCORDANCE", "1", "Float", "Concordance of repeat unit.", true);
        RU_COUNTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU_COUNTS", "2", "Integer", "Number of exact repeat units and total number of repeat units.", true);
        FLANKS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FLANKS", "2", "Integer", "Left and right flank positions of the Indel, left/right alignment invariant, not necessarily equal to POS.", true);

        if (method=="f")
        {
            FZ_RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_RL", "1", "Float", "Fuzzy reference repeat unit length", true);
            FZ_LL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_LL", "1", "Float", "Fuzzy longest repeat unit length", true);
            FZ_CONCORDANCE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_CONCORDANCE", "1", "Float", "Fuzzy concordance of repeat unit.", true);
            FZ_RU_COUNTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_RU_COUNTS", "2", "Integer", "Fuzzy number of exact repeat units and total number of repeat units.", true);
            FZ_FLANKS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_FLANKS", "2", "Integer", "Fuzzy left and right flank positions of the Indel, left/right alignment invariant, not necessarily equal to POS.", true);
        }

        TR = bcf_hdr_append_info_with_backup_naming(odw->hdr, "TR", "1", "String", "Tandem repeat associated with this indel.", true);

        bcf_hdr_append(odw->hdr, "##INFO=<ID=LARGE_REPEAT_REGION,Number=0,Type=Flag,Description=\"Very large repeat region, vt only detects up to 1000bp long regions.\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=FLANKSEQ,Number=1,Type=String,Description=\"Flanking sequence 10bp on either side of detected repeat region.\">");
                
        //helper variable initialization for adding genotype fields for additional vntr records
        if (annotation_mode=="v")
        {
            no_samples = bcf_hdr_nsamples(odw->hdr);
            gts = (int32_t*) malloc(no_samples*sizeof(int32_t));
            for (uint32_t i=0; i<no_samples; ++i)
            {
                gts[i] = 0;
            }
        }
        else
        {
            no_samples = 0;
            gts = NULL;
        }

        s = {0,0,0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants_annotated = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip(ref_fasta_file);
        va = new VNTRAnnotator(ref_fasta_file, debug);
        fai = fai_load(ref_fasta_file.c_str());

    }

    void print_options()
    {
        std::clog << "annotate_indels v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file          " << output_vcf_file << "\n";
        std::clog << "         [m] method of VNTR detection " << method << "\n";
        std::clog << "         [a] mode of annotation       " << annotation_mode << "\n";
        print_boo_op("         [d] debug                    ", debug);
        print_ref_op("         [r] ref FASTA file           ", ref_fasta_file);
        print_boo_op("         [x] override tag             ", override_tag);
        print_int_op("         [i] intervals                ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants annotated   " << no_variants_annotated << "\n";
        std::clog << "\n";
    }

    /**
     * Inserts a VNTR record.
     * Returns true if successful.
     */
    bool insert_vntr_record_into_buffer(VNTR& vntr)
    {
        if (method=="e")
        {
            std::list<VNTR>::iterator i = vntr_buffer.begin();
            while(i!=vntr_buffer.end())
            {
                VNTR& cvntr = *i;
    
                if (vntr.rid > cvntr.rid)
                {
                    vntr_buffer.insert(i, vntr);
                    return true;
                }
                else if (vntr.rid == cvntr.rid)
                {
                    if (vntr.rbeg1 > cvntr.rbeg1)
                    {
                        vntr_buffer.insert(i, vntr);
                        return true;
                    }
                    else if (vntr.rbeg1 == cvntr.rbeg1)
                    {
                        if (vntr.rend1 > cvntr.rend1)
                        {
                            vntr_buffer.insert(i, vntr);
                            return true;
                        }
                        else if (cvntr.rend1 == vntr.rend1)
                        {
                            if (cvntr.motif > vntr.motif)
                            {
                                vntr_buffer.insert(i, vntr);
                                return true;
                            }
                            else if (cvntr.motif == vntr.motif)
                            {
                                //do not insert
                                return false;
                            }
                            else // cvntr.motif > vntr.motif
                            {
                                ++i;
                            }
                        }
                        else // cvntr.rend1 > vntr.rend1
                        {
                            ++i;
                        }
                    }
                    else //vntr.rbeg1 < cvntr.rbeg1
                    {
                        ++i;
                    }
                }
                else //vntr.rid < cvntr.rid is impossible if input file is ordered.
                {
                    fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
                    exit(1);
                }
            }
        }
        else if (method=="f")
        {
            std::list<VNTR>::iterator i = vntr_buffer.begin();
            while(i!=vntr_buffer.end())
            {
                VNTR& cvntr = *i;
    
                if (vntr.rid > cvntr.rid)
                {
                    vntr_buffer.insert(i, vntr);
                    return true;
                }
                else if (vntr.rid == cvntr.rid)
                {
                    if (vntr.fuzzy_rbeg1 > cvntr.fuzzy_rbeg1)
                    {
                        vntr_buffer.insert(i, vntr);
                        return true;
                    }
                    else if (vntr.fuzzy_rbeg1 == cvntr.fuzzy_rbeg1)
                    {
                        if (vntr.fuzzy_rend1 > cvntr.fuzzy_rend1)
                        {
                            vntr_buffer.insert(i, vntr);
                            return true;
                        }
                        else if (cvntr.fuzzy_rend1 == vntr.fuzzy_rend1)
                        {
                            if (cvntr.motif > vntr.motif)
                            {
                                vntr_buffer.insert(i, vntr);
                                return true;
                            }
                            else if (cvntr.motif == vntr.motif)
                            {
                                //do not insert
                                return false;
                            }
                            else // cvntr.motif > vntr.motif
                            {
                                ++i;
                            }
                        }
                        else // cvntr.fuzzy_rend1 > vntr.fuzzy_rend1
                        {
                            ++i;
                        }
                    }
                    else //vntr.fuzzy_rbeg1 < cvntr.fuzzy_rbeg1
                    {
                        ++i;
                    }
                }
                else //vntr.rid < cvntr.rid is impossible if input file is ordered.
                {
                    fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
                    exit(1);
                }
            }
        }
        
        vntr_buffer.push_back(vntr);
        return true;
    }

    /**
     * Flush variant buffer.
     */
    void flush_vntr_buffer(bcf1_t* v)
    {
        if (vntr_buffer.empty())
        {
            return;
        }

        int32_t rid = bcf_get_rid(v);
        int32_t pos1 = bcf_get_pos1(v);

        //search for vntr to start deleting from.
        std::list<VNTR>::iterator i = vntr_buffer.begin();
        while(i!=vntr_buffer.end())
        {
            VNTR& vntr = *i;

//            std::cerr << vntr.rid << " " << rid << "\n";

            if (vntr.rid < rid)
            {
                break;
            }
            else if (vntr.rid == rid)
            {
                if (method=="e")
                {
                    if (vntr.rend1 < pos1-2000)
                    {
                        break;
                    }
                }
                else if (method=="f")
                {
                    if (vntr.fuzzy_rend1 < pos1-2000)
                    {
                        break;
                    }
                }
            }
            else //rid < vntr.rid is impossible
            {
                fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
                exit(1);
            }

            ++i;
        }

        while (i!=vntr_buffer.end())
        {
            i = vntr_buffer.erase(i);
        }
    }

    /**
     * Creates a VNTR record.
     */
    void create_vntr_record(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
    {
        VNTR& vntr = variant.vntr;

        //shared fields
        bcf_set_rid(v, variant.rid);
        bcf_update_info_string(h, v, MOTIF.c_str(), vntr.motif.c_str());
        bcf_update_info_string(h, v, RU.c_str(), vntr.ru.c_str());


        if (method=="e")
        {
            bcf_set_pos1(v, vntr.rbeg1);

            s.l = 0;
            kputs(vntr.repeat_tract.c_str(), &s);
            kputc(',', &s);
            kputs("<VNTR>", &s);
            bcf_update_alleles_str(h, v, s.s);
            bcf_update_info_string(h, v, MOTIF.c_str(), vntr.motif.c_str());
            bcf_update_info_string(h, v, RU.c_str(), vntr.ru.c_str());

            bcf_update_info_float(h, v, RL.c_str(), &vntr.rl, 1);
            bcf_update_info_float(h, v, LL.c_str(), &vntr.ll, 1);
            int32_t ru_count[2] = {vntr.no_exact_ru, vntr.total_no_ru};
            bcf_update_info_float(h, v, CONCORDANCE.c_str(), &vntr.motif_concordance, 1);
            bcf_update_info_int32(h, v, RU_COUNTS.c_str(), &ru_count, 2);
            
            //update flanking sequences
            std::string flanks;
            int32_t seq_len;
            char* seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.rbeg1-10-1, variant.vntr.rbeg1-1-1, &seq_len);
            flanks.assign(seq);
            free(seq);
            flanks.append(1, '[');
            seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.rbeg1-1, variant.vntr.rend1-1, &seq_len);
            flanks.append(seq);
            free(seq);
            flanks.append(1, ']');
            seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.rend1+1-1, variant.vntr.rend1+10-1, &seq_len);
            flanks.append(seq);
            free(seq);
            bcf_update_info_string(h, v, "FLANKSEQ", flanks.c_str()); 
            
            if (vntr.is_large_repeat_tract)
            {
                bcf_update_info_flag(h, v, "LARGE_REPEAT_REGION", NULL, 1);
            }
        }
        else if (method=="f")
        {
            bcf_set_pos1(v, vntr.fuzzy_rbeg1);

            s.l = 0;
            kputs(vntr.fuzzy_repeat_tract.c_str(), &s);
            kputc(',', &s);
            kputs("<VNTR>", &s);
            bcf_update_alleles_str(h, v, s.s);
            bcf_update_info_string(h, v, MOTIF.c_str(), vntr.motif.c_str());
            bcf_update_info_string(h, v, RU.c_str(), vntr.ru.c_str());

            bcf_update_info_float(h, v, FZ_CONCORDANCE.c_str(), &vntr.fuzzy_motif_concordance, 1);

            bcf_update_info_float(h, v, FZ_RL.c_str(), &vntr.fuzzy_rl, 1);
            bcf_update_info_float(h, v, FZ_LL.c_str(), &vntr.fuzzy_ll, 1);
            
            int32_t flank_pos1[2] = {variant.vntr.rbeg1-1, variant.vntr.rend1+1};
            bcf_update_info_int32(h, v, FLANKS.c_str(), &flank_pos1, 2);
            
            int32_t fuzzy_flank_pos1[2] = {variant.vntr.fuzzy_rbeg1-1, variant.vntr.fuzzy_rend1+1};
            bcf_update_info_int32(h, v, FZ_FLANKS.c_str(), &fuzzy_flank_pos1, 2);
            int32_t ru_count[2] = {vntr.fuzzy_no_exact_ru, vntr.fuzzy_total_no_ru};
            bcf_update_info_int32(h, v, FZ_RU_COUNTS.c_str(), &ru_count, 2);
            
            //update flanking sequences
            std::string flanks;
            int32_t seq_len;
            char* seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.fuzzy_rbeg1-10-1, variant.vntr.fuzzy_rbeg1-1-1, &seq_len);
            flanks.assign(seq);
            free(seq);
            flanks.append(1, '[');
            seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.fuzzy_rbeg1-1, variant.vntr.fuzzy_rend1-1, &seq_len);
            flanks.append(seq);
            free(seq);
            flanks.append(1, ']');
            seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.fuzzy_rend1+1-1, variant.vntr.fuzzy_rend1+10-1, &seq_len);
            flanks.append(seq);
            free(seq);
            bcf_update_info_string(h, v, "FLANKSEQ", flanks.c_str()); 
                                    
            if (vntr.is_large_repeat_tract)
            {
                bcf_update_info_flag(h, v, "LARGE_REPEAT_REGION", NULL, 1);
            }
        }

        //individual fields - just set GT
        bcf_update_genotypes(h, v, gts, no_samples);
    }

    bool genotype_str(char* repeat_tract, char* ru, float& score, float& discordance, bool& exact)
    {
        //check if start is not the same as the RU.
        int32_t rlen = strlen(repeat_tract);
        int32_t mlen = strlen(ru);

        if (strncmp(repeat_tract, ru, mlen)!=0)
        {
            return false;
        }

        return true;
    }

    void annotate_indels()
    {
        odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odw->hdr;
        Variant variant;
        kstring_t old_alleles = {0,0,0};

        int32_t no_exact = 0;
        int32_t no_inexact = 0;

        while (odr->read(v))
        {
            if (filter_exists)
            {
                vm->classify_variant(h, v, variant);
                if (!filter.apply(h, v, &variant, false))
                {
                    continue;
                }
            }

            bcf_unpack(v, BCF_UN_STR);
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            if (vtype&VT_INDEL)
            {
//                bcf_print_liten(h,v);
                
                flush_vntr_buffer(v);

                //  bcf_print(odr->hdr, v);
                va->annotate(odr->hdr, v, variant, method);

                if (annotation_mode=="v")
                {
                    if (method=="e")
                    {
                        int32_t flank_pos1[2] = {variant.vntr.rbeg1-1, variant.vntr.rend1+1};
                        bcf_update_info_int32(h, v, FLANKS.c_str(), &flank_pos1, 2);

                        //update flanking sequences
                        std::string flanks;
                        int32_t seq_len;
                        char* seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.rbeg1-10-1, variant.vntr.rbeg1-1-1, &seq_len);
                        flanks.assign(seq);
                        free(seq);
                        flanks.append(1, '[');
                        seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.rbeg1-1, variant.vntr.rend1-1, &seq_len);
                        flanks.append(seq);
                        free(seq);
                        flanks.append(1, ']');
                        seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.rend1+1-1, variant.vntr.rend1+10-1, &seq_len);
                        flanks.append(seq);
                        free(seq);
                        bcf_update_info_string(h, v, "FLANKSEQ", flanks.c_str());  
                    }
                    else if (method=="f")
                    {
                        int32_t flank_pos1[2] = {variant.vntr.rbeg1-1, variant.vntr.rend1+1};
                        bcf_update_info_int32(h, v, FLANKS.c_str(), &flank_pos1, 2);

                        int32_t fuzzy_flank_pos1[2] = {variant.vntr.fuzzy_rbeg1-1, variant.vntr.fuzzy_rend1+1};
                        bcf_update_info_int32(h, v, FZ_FLANKS.c_str(), &fuzzy_flank_pos1, 2);
                        
                        //update flanking sequences
                        std::string flanks;
                        int32_t seq_len;
                        char* seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.fuzzy_rbeg1-10-1, variant.vntr.fuzzy_rbeg1-1-1, &seq_len);
                        flanks.assign(seq);
                        free(seq);
                        flanks.append(1, '[');
                        seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.fuzzy_rbeg1-1, variant.vntr.fuzzy_rend1-1, &seq_len);
                        flanks.append(seq);
                        free(seq);
                        flanks.append(1, ']');
                        seq = faidx_fetch_seq(fai, variant.chrom.c_str(), variant.vntr.fuzzy_rend1+1-1, variant.vntr.fuzzy_rend1+10-1, &seq_len);
                        flanks.append(seq);
                        free(seq);
                        bcf_update_info_string(h, v, "FLANKSEQ", flanks.c_str());                        
                    }

                    if (va->is_vntr(variant, vntr_classification, method))
                    {
                        if (method=="e")
                        {
                            variant.get_vntr_string(&s);
                        }
                        else if (method=="f")
                        {
                            variant.get_fuzzy_vntr_string(&s);
                        }
                        
                        bcf_update_info_string(h, v, TR.c_str(), s.s);
                        odw->write(v);
                        v = odw->get_bcf1_from_pool();
                            
                        if (insert_vntr_record_into_buffer(variant.vntr))
                        {
                            create_vntr_record(odr->hdr, v, variant);
                            odw->write(v);
                            v = odw->get_bcf1_from_pool();
                        }
                    }
                    else
                    {
                        odw->write(v);
                        v = odw->get_bcf1_from_pool();
                    }
                }

//                std::cerr << "vntr_buffer size " << vntr_buffer.size() << "\n";

                ++no_variants_annotated;
            }
            else if (vtype&VT_VNTR)
            {
                char** allele = bcf_get_allele(v);
                char* ru = NULL;
                int32_t n = 0;
                if (bcf_get_info_string(odw->hdr, v, "RU", &ru, &n)>0)
                {
                    float genotype = 0;
                    float discordance = 0;
                    bool exact = false;

                    bool res = genotype_str(allele[0], ru, genotype, discordance, exact);

                    if (!res)
                    {
                        ++no_inexact;

                        bcf_print(odw->hdr, v);
                        std::cerr << "genotype    : " << genotype << "\n";
                        std::cerr << "discordance : " << discordance << "\n";
                        std::cerr << "exact       : " << exact << "\n";
                    }

                    ++no_exact;

                    free(ru);
                }
            }
            else // SNP
            {
                //update flanking sequences
                std::string flanks;
                int32_t seq_len;
                std::string chrom = bcf_get_chrom(h, v);
                int32_t pos1 = bcf_get_pos1(v);
                    
                char* seq = faidx_fetch_seq(fai, chrom.c_str(), pos1-10-1, pos1-1-1, &seq_len);
                flanks.assign(seq);
                free(seq);
                flanks.append(1, '[');
                char** alleles = bcf_get_allele(v);
                flanks.append(alleles[0]);
                flanks.append(1, '/');
                flanks.append(alleles[1]);
                flanks.append(1, ']');
                seq = faidx_fetch_seq(fai, chrom.c_str(), pos1+1-1, pos1+10-1, &seq_len);
                flanks.append(seq);
                free(seq);
                bcf_update_info_string(h, v, "FLANKSEQ", flanks.c_str()); 
                
                odw->write(v);
                v = odw->get_bcf1_from_pool();
            }
        }

//        std::cerr << "no inexact : " << no_inexact << "\n";
//        std::cerr << "no exact : " << no_exact << "\n";

        odw->close();
        odr->close();
    };

    private:
};
}

void annotate_indels(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.annotate_indels();
    igor.print_stats();
};
