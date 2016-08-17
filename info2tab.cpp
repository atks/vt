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

#include "info2tab.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_text_file;
    std::vector<GenomeInterval> intervals;
    std::string ref_fasta_file;
    std::vector<std::string> info_tags;
    bool debug;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    uint32_t no_variants;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "extracts info fields to a tab delimited file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_text_file("o", "o", "output tab delimited file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_info_tags("t", "t", "list of info tags to be removed []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::SwitchArg arg_smart("s", "s", "smart decomposition [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_text_file = arg_output_text_file.getValue();
            fexp = arg_fexp.getValue();
            parse_string_list(info_tags, arg_info_tags.getValue());
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

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip();
    }

    void info2tab()
    {
        bcf_hdr_t* h = odr->hdr;
        bcf1_t* v = bcf_init();
        Variant variant;

        FILE *out = fopen(output_text_file.c_str(), "w");

//            fprintf(out, "\n");
//            fprintf(out, "\\begin{frame}{Data set summary}\n");
//            fprintf(out, "\\resizebox{\\linewidth}{!}{\n");
//            fprintf(out, "\\rowcolors{2}{blue!25}{blue!10}\n");
//            fprintf(out, "\\begin{tabular}{rrrrr}\n");
//            fprintf(out, "\\rowcolor{blue!50}\n");
//            fprintf(out, "%s & no. indels & ins/del & ins & del\\\\ \n", dataset_labels[i].c_str());
//            fprintf(out, "A-B & %d & %.1f & %d & %d\\\\ \n",  stats[i].a, (float)stats[i].a_ins/(stats[i].a_del), stats[i].a_ins, stats[i].a_del);
//            fprintf(out, "A\\&B & %d & %.1f & %d & %d\\\\ \n",  stats[i].ab, (float)stats[i].ab_ins/(stats[i].ab_del), stats[i].ab_ins, stats[i].ab_del);
//            fprintf(out, "B-A & %d & %.1f & %d & %d\\\\ \n",  stats[i].b, (float)stats[i].b_ins/(stats[i].b_del), stats[i].b_ins, stats[i].b_del);
//            fprintf(out, " &  &  & &  \\\\ \n");
//            fprintf(out, " Precision & %.2f\\%% &  &  & \\\\ \n", 100*(float)stats[i].ab/(stats[i].a+stats[i].ab));
//            fprintf(out, " Sensitivity & %.2f\\%% &  &  &  \\\\ \n", 100*(float)stats[i].ab/(stats[i].b+stats[i].ab));
//            fprintf(out, "\\end{tabular}}\n");
//            fprintf(out, "\\end{frame}\n");

        //get the types of each info field
        std::vector<std::string> info_tag_str;
        std::vector<int32_t> info_tag_id;
        std::vector<int32_t> info_tag_vlen;
        std::vector<int32_t> info_tag_type;
        std::vector<int32_t> info_tag_num;
            
        for (uint32_t i=0; i<info_tags.size(); ++i)
        {
            int32_t id = bcf_hdr_id2int(h, BCF_DT_ID, info_tags[i].c_str());
            int32_t vlen = bcf_hdr_id2length(h, BCF_HL_INFO, info_tag_id[i]);
            int32_t type = bcf_hdr_id2type(h, BCF_HL_INFO, info_tag_id[i]);
            int32_t num = bcf_hdr_id2number(h, BCF_HL_INFO, info_tag_id[i]);               
       
            if (vlen==BCF_VL_FIXED && num!=1)
            {
                info_tag_id.push_back(id);
                info_tag_vlen.push_back(vlen);
                info_tag_type.push_back(type);
                info_tag_num.push_back(num);        
            }    
            else
            {
                notice("info2tab does not support id:%s, vlen=%s, type=%s, num=%d", bcf_hdr_vl2str(id).c_str(), bcf_hdr_vl2str(id).c_str(), bcf_hdr_vl2str(id).c_str(), num);
            }
        }

        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_INFO);

            if (filter_exists)
            {
                vm->classify_variant(h, v, variant);
                if (!filter.apply(h, v, &variant, false))
                {
                    continue;
                }
            }

            int32_t ret = 0;
            for (uint32_t i=0; i<info_tag_str.size(); ++i)
            {
                if (i)
                {
                    fprintf(out, "\t");
                }
                
                int32_t id = info_tag_id[i];
                int32_t vlen = info_tag_vlen[i];
                int32_t type = info_tag_type[i];
                int32_t num = info_tag_num[i];

                if (vlen==BCF_VL_G)
                {
                    if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                    {
                        //to be implemented
                    }
                    else if (type==BCF_BT_FLOAT)
                    {
                        //to be implemented
                    }
                    else if (type==BCF_BT_CHAR)
                    {
                        //to be implemented
                    }
                }
                else if (vlen == BCF_VL_A)
                {
                    if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                    {
                        
                    }
                    else if (type==BCF_BT_FLOAT)
                    {
                       
                    }
                    else if (type==BCF_BT_CHAR)
                    {
                        //to be implemented
                    }
                }
                else if (vlen == BCF_VL_R)
                {
                    if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                    {
                     
                    }
                    else if (type==BCF_BT_FLOAT)
                    {
                
                    }
                    else if (type==BCF_BT_CHAR)
                    {
                        //to be implemented
                    }
                }
                else if (vlen == BCF_VL_FIXED)
                {
                    
                }
                else if (vlen == BCF_VL_VAR)
                {
                    
                }
                
                if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                {
                 
                }
                else if (type==BCF_BT_FLOAT)
                {
            
                }
                else if (type==BCF_BT_CHAR)
                {
                    //to be implemented
                }
                
                                    
            }

            ++no_variants;

        }

        fclose(out);
        odr->close();
    };

    void print_options()
    {
        std::clog << "info2tab v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output text file      " << output_text_file << "\n";
        print_strvec("         [t] info tags                        ", info_tags);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. variants                 : " << no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void info2tab(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.info2tab();
    igor.print_stats();
};
