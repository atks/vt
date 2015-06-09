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
    //modes of annotation
    //1. e: exact alignment, motif tree - VMOTIF, VSCORE
    //2. f: fuzzy alignment, motif tree - VMOTIF, VSCORE
    //3. x: fuzzy alignment, motif tree, robust detection of flanks - VMOTIF, VSCORE, LFLANK, RFLANK
    std::string mode;
    bool override_tag;
    uint32_t alignment_penalty;

    std::string MOTIF;
    std::string SCORE;        
    std::string RU;
    std::string RL;
            
    bool debug;

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
            std::string desc = "annotates indels with STR type information - repeat tract length, repeat motif, flank information";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_mode("m", "m", "mode [x]\n"
                 "              e : determine by exact alignment.\n"
                 "              f : determine by fuzzy alignment.\n"
                 "              x : using HMMs",
                 false, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_alignment_penalty("p", "p", "alignment penalty [0]\n", false, 0, "int", cmd);     
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
            TCLAP::SwitchArg arg_override_tag("x", "x", "override tags [false]", cmd, false);
            
            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            mode = arg_mode.getValue();
            override_tag = arg_override_tag.getValue();
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
        if (mode!="e" && mode!="f" && mode!="x")
        {
            fprintf(stderr, "[%s:%d %s] Not a valid mode of annotation: %s\n", __FILE__,__LINE__,__FUNCTION__, mode.c_str());
            exit(1);
        }   
       
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);
        
        MOTIF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "MOTIF", "1", "String", "Canonical Motif in an VNTR or Homopolymer", true);
        SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "SCORE", "1", "Float", "Score of repeat unit", true);
        RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU", "1", "String", "Repeat unit in a VNTR or Homopolymer", true);
        RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RL", "1", "Float", "Repeat Unit Length", true);

        bcf_hdr_append(odw->hdr, "##INFO=<ID=OLD_VARIANT,Number=1,Type=String,Description=\"Original chr:pos:ref:alt encoding\">\n");

//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_LFLANK,Number=1,Type=String,Description=\"Right Flank Sequence\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RFLANK,Number=1,Type=String,Description=\"Left Flank Sequence\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_LFLANKPOS,Number=2,Type=Integer,Description=\"Positions of left flank\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_RFLANKPOS,Number=2,Type=Integer,Description=\"Positions of right flank\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_MOTIF_DISCORDANCE,Number=1,Type=Integer,Description=\"Descriptive Discordance for each reference repeat unit.\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_MOTIF_COMPLETENESS,Number=1,Type=Integer,Description=\"Descriptive Discordance for each reference repeat unit.\">");
//        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_STR_CONCORDANCE,Number=1,Type=Float,Description=\"Overall discordance of RUs.\">");

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
        print_ref_op("         [r] ref FASTA file        ", ref_fasta_file);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants annotated   " << no_variants_annotated << "\n";
        std::clog << "\n";
    }

    void update_vntr_info(bcf_hdr_t* h, bcf1_t *v, Variant& variant)
    {
        if (variant.vntr.motif!="") 
        {
            bcf_update_info_string(odw->hdr, v, MOTIF.c_str(), variant.vntr.motif.c_str());
        }
        
        std::cerr << "RU: " << variant.vntr.ru << "\n";
        
        if (variant.vntr.ru!="") 
        {
            std::cerr << "RU: " << variant.vntr.ru << "\n";
                
            bcf_update_info_string(odw->hdr, v, RU.c_str(), variant.vntr.ru.c_str());
        }
        if (variant.vntr.motif_score>=0) 
        {    
            bcf_update_info_float(odw->hdr, v, SCORE.c_str(), &variant.vntr.motif_score, 1);
        }

        
//        bcf_update_info_string(odw->hdr, v, "VRU", variant.eru.c_str());
//        int32_t rl = variant.eregion.end1-variant.eregion.beg1-1;
//        bcf_update_info_int32(odw->hdr, v, "VRL", &rl, 1);
//        int32_t irl = variant.iregion.end1-variant.iregion.beg1-1;
//        bcf_update_info_int32(odw->hdr, v, "IRL", &irl, 1);
//
//        if (irl!=rl)
//        {
//            int32_t irg[2] = {variant.iregion.beg1, variant.iregion.end1};
//            bcf_update_info_int32(odw->hdr, v, "IRG", &irg, 2);
//            int32_t len = 0;
//            const char* chrom = bcf_get_chrom(odr->hdr, v);
//            char* seq = faidx_fetch_seq(fai, chrom, irg[0]-1, irg[1]-1, &len);
//            bcf_update_info_string(odw->hdr, v, "ISQ", seq);
//            if (len) free(seq);
//        }
//
//                //annotate old alleles
//        old_alleles.l = 0;
//        bcf_variant2string(odw->hdr, v, &old_alleles);
//        bcf_update_info_string(odw->hdr, v, "OLD_VARIANT", old_alleles.s);
//
//        //update alleles
//        bcf_set_pos1(v, variant.vntr.pos1);
//        std::string new_alleles = variant.vntr.ref;
//        new_alleles += ",<VNTR>";
//        bcf_update_alleles_str(odw->hdr, v, new_alleles.c_str());

    }

    void annotate_indels()
    {
        odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        Variant variant;
        kstring_t old_alleles = {0,0,0};
        
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
            if (vtype&VT_INDEL || vtype&VT_VNTR)
            {
//                bcf_print(odr->hdr, v);
                //annotate indel like variant
                va->annotate(odr->hdr, v, variant, mode);
                update_vntr_info(odr->hdr, v, variant);
   
                ++no_variants_annotated;
            }

            odw->write(v);
            v = odw->get_bcf1_from_pool();
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
