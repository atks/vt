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
    //methods of detection 
    //1. e: exact alignment, motif tree - VMOTIF, VSCORE
    //2. f: fuzzy alignment, motif tree - VMOTIF, VSCORE
    //3. x: fuzzy alignment, motif tree, robust detection of flanks - VMOTIF, VSCORE, LFLANK, RFLANK
    std::string method;
    //modes of annotation
    //1. a: update ALLELES fields; update INFO fields; RL, MOTIF, RU, FLANKS, OLD_VARIANT
    //2. i: update INFO fields; RL, MOTIF, RU, FLANKS, OLD_VARIANT
    //3. n: insert new record and update OLD_VARIANT
    //      a. NEW RECORD: update ALLELES fields; update INFO fields; RL, MOTIF, RU, FLANKS, OLD_VARIANT
    //      b. update INFO fields; VNTR
    std::string annotation_mode;
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
            TCLAP::ValueArg<std::string> arg_method("m", "m", "mode [x]\n"
                 "              e : determine by exact alignment.\n"
                 "              f : determine by fuzzy alignment.\n"
                 "              p : determine by penalized fuzzy alignment.\n"
                 "              h : using HMMs"
                 "              x : integrated models",
                 false, "e", "str", cmd);
            TCLAP::ValueArg<std::string> arg_annotation_mode("a", "a", "mode [x]\n"
                 "              e : output VNTR variant .\n"
                 "              f : determine by fuzzy alignment.\n"
                 "              p : determine by penalized fuzzy alignment.\n"
                 "              h : using HMMs"
                 "              x : integrated models",
                 false, "e", "str", cmd);
//            TCLAP::ValueArg<std::string> arg_classification("a", "a", "mode [x]\n"g
//                 "              e : output VNTR variant .\n"
//                 "              f : determine by fuzzy alignment.\n"
//                 "              p : determine by penalized fuzzy alignment.\n"
//                 "              h : using HMMs"
//                 "              x : integrated models",
//                 false, "e", "str", cmd);     
            TCLAP::ValueArg<uint32_t> arg_alignment_penalty("p", "p", "alignment penalty [0]\n", false, 0, "int", cmd);
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
            override_tag = arg_override_tag.getValue();
            add_vntr_record = arg_add_vntr_record.getValue();
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
        if (method!="e" && method!="f" && method!="p" && method!="x")
        {
            fprintf(stderr, "[%s:%d %s] Not a valid mode of annotation: %s\n", __FILE__,__LINE__,__FUNCTION__, method.c_str());
            exit(1);
        }

        if (annotation_mode!="e" && annotation_mode!="f" && annotation_mode!="p" && annotation_mode!="x")
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
        REF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "REF", "1", "String", "Repeat tract on the reference sequence", true);
        REFPOS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "REFPOS", "1", "Integer", "Start position of repeat tract", true);
        SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "SCORE", "1", "Float", "Score of repeat unit", true);
        TR = bcf_hdr_append_info_with_backup_naming(odw->hdr, "TR", "1", "String", "Tandem repeat representation", true);

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
        std::clog << "options:     input VCF file(s)     " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [m] method                " << method << "\n";
        std::clog << "         [a] annotation_mode       " << annotation_mode << "\n";
        print_boo_op("         [d] debug                 ", debug);
        print_ref_op("         [r] ref FASTA file        ", ref_fasta_file);
        print_boo_op("         [x] override tag          ", override_tag);        
        print_boo_op("         [v] add vntr record       ", add_vntr_record);        
        
        print_int_op("         [i] intervals             ", intervals);
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

    /**
     * Updates an Indel record with:
     * a. MOTIF
     * b. RU
     * c. RL
     */
    void update_indel_record(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
    {
        VNTR& vntr = variant.vntr;
         
        bcf_set_rid(v, variant.rid);
        bcf_set_pos1(v, vntr.rbeg1);
        s.l = 0;
        kputs(vntr.repeat_tract.c_str(), &s);
        kputc(',', &s);
        kputs("<VNTR>", &s);
        
        if ((vntr.rl-vntr.motif.size())>=6)
        {    
            bcf_update_alleles_str(h, v, s.s);
        }
        bcf_update_info_string(h, v, MOTIF.c_str(), vntr.motif.c_str());
        bcf_update_info_string(h, v, RU.c_str(), vntr.ru.c_str());
        bcf_update_info_float(h, v, RL.c_str(), &vntr.rl, 1);
    }

    /**
     * Updates an Indel record with an overlapping VNTR:
     */
    void update_indel_record_with_vntr_annotation(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
    {
        VNTR& vntr = variant.vntr;
         
        bcf_update_info_string(h, v, MOTIF.c_str(), vntr.motif.c_str());
        bcf_update_info_string(h, v, RU.c_str(), vntr.ru.c_str());
        bcf_update_info_float(h, v, RL.c_str(), &vntr.rl, 1);
        
        variant.get_vntr_string(&s);
        bcf_update_info_string(h, v, "TR", s.s);    
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
