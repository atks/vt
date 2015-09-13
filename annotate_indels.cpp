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

    std::string method;  //methods of detection
    std::string annotation_mode;  //modes of annotation

    bool override_tag;
    uint32_t alignment_penalty;
    bool add_vntr_record;

    std::string MOTIF;
    std::string RU;
    std::string RL;
    std::string REF;
    std::string REFPOS;
    std::string SCORE;
    std::string TR;

    int32_t vntr_classification;

    //helper variables for populating additional VNTR records
    uint32_t no_samples;
    int32_t* gts;

    kstring_t s;

    bool debug;

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
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file [e]", true, "e", "str", cmd);
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
            TCLAP::ValueArg<std::string> arg_method("m", "m", "mode [e]\n"
                 "              exact : determine by exact alignment.\n"
                 "              fuzzy : determine by fuzzy alignment.",
                 false, "e", "str", cmd);

            TCLAP::ValueArg<uint32_t> arg_alignment_penalty("p", "p", "alignment penalty [0]", false, 0, "int", cmd);
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
            //add_vntr_record = arg_add_vntr_record.getValue();
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

    ~Igor() {};

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

        add_vntr_record = annotation_mode=="v" ? true: false;
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

        MOTIF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "MOTIF", "1", "String", "Canonical motif in an VNTR or homopolymer", true);
        RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU", "1", "String", "Repeat unit in a VNTR or homopolymer", true);
        RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RL", "1", "Float", "Repeat unit length", true);
//        REF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "REF", "1", "String", "Repeat tract on the reference sequence", true);
//        REFPOS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "REFPOS", "1", "Integer", "Start position of repeat tract", true);
//        SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "SCORE", "1", "Float", "Score of repeat unit", true);
        TR = bcf_hdr_append_info_with_backup_naming(odw->hdr, "TR", "1", "String", "Tandem repeat representation", true);

        TR = bcf_hdr_append_info_with_backup_naming(odw->hdr, "END", "1", "Integer", "End position of the variant described in this record", true);

//        bcf_hdr_append(odw->hdr, "##INFO=<ID=OLD_VARIANT,Number=1,Type=String,Description=\"Original chr:pos:ref:alt encoding\">\n");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_LFLANK,Number=1,Type=String,Description=\"Right Flank Sequence\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RFLANK,Number=1,Type=String,Description=\"Left Flank Sequence\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_LFLANKPOS,Number=2,Type=Integer,Description=\"Positions of left flank\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RFLANKPOS,Number=2,Type=Integer,Description=\"Positions of right flank\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_MOTIF_DISCORDANCE,Number=1,Type=Integer,Description=\"Descriptive Discordance for each reference repeat unit.\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_MOTIF_COMPLETENESS,Number=1,Type=Integer,Description=\"Descriptive Discordance for each reference repeat unit.\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_STR_CONCORDANCE,Number=1,Type=Float,Description=\"Overall discordance of RUs.\">");

        //helper variable initialization for adding additional vntr records
        if (add_vntr_record)
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
        print_boo_op("         [v] add vntr record          ", add_vntr_record);
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
     */
    void insert_vntr_record(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
    {
        VNTR& vntr = variant.vntr;

        if (va->is_vntr(variant, vntr_classification))
        {
            //shared fields
            bcf_set_rid(v, variant.rid);
            bcf_set_pos1(v, vntr.rbeg1);
            s.l = 0;
            kputs(vntr.repeat_tract.c_str(), &s);
            kputc(',', &s);
            kputs("<VNTR>", &s);
            bcf_update_alleles_str(h, v, s.s);
            bcf_update_info_string(h, v, MOTIF.c_str(), vntr.motif.c_str());
            bcf_update_info_string(h, v, RU.c_str(), vntr.ru.c_str());
            bcf_update_info_float(h, v, RL.c_str(), &vntr.rl, 1);

            //individual fields - just set GT
            bcf_update_genotypes(h, v, gts, no_samples);
        }
    }

    /**
     * Updates an Indel record with an overlapping VNTR and end of allelelic region.
     */
    void update_indel_record(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
    {
        VNTR& vntr = variant.vntr;

        if (va->is_vntr(variant, vntr_classification))
        {
            variant.get_vntr_string(&s);
            bcf_update_info_string(h, v, "TR", s.s);
        }

        bcf_update_info_int32(h, v, "END", &vntr.rend1, 1);
    }

    void annotate_indels()
    {
        odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odw->hdr;
        Variant variant;
        kstring_t old_alleles = {0,0,0};

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
            if (vtype&VT_INDEL || vtype&VT_VNTR)
            {
//                bcf_print(odr->hdr, v);
                va->annotate(odr->hdr, v, variant, method);


                //setup a pass filter for the VNTR classification

                //add a filter to indicate VNTR fitness.
                update_indel_record(h, v, variant);
                odw->write(v);
                v = odw->get_bcf1_from_pool();

                if (add_vntr_record)
                {
                    //check if the record exists before ...
                    //should be done with a little buffer?
                    //done with hash?
                    //10% overlap as observed

                    insert_vntr_record(odr->hdr, v, variant);
                    odw->write(v);
                    v = odw->get_bcf1_from_pool();
                }

                ++no_variants_annotated;
            }
            else
            {
                odw->write(v);
                v = odw->get_bcf1_from_pool();
            }
        }

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
