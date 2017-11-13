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
    std::vector<std::string> filter_tags;
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
            std::string desc = "extracts INFO and FILTER fields to a tab delimited file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_text_file("o", "o", "output tab delimited file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_info_tags("t", "t", "list of info tags to be extracted []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_filter_tags("u", "u", "list of filter tags to be extracted []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_text_file = arg_output_text_file.getValue();
            fexp = arg_fexp.getValue();
            std::cerr << arg_info_tags.getValue() << "\n";
            parse_string_list(filter_tags, arg_filter_tags.getValue());
            parse_string_list(info_tags, arg_info_tags.getValue());
            if (filter_tags.size()==0 && info_tags.size()==0)
            {
                notice("no filter or info tags specified\n");
                exit(1);
            }
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

        FILE *out = NULL;
        if (output_text_file=="-")
        {
            out = stdout;
        }
        else
        {
            out = fopen(output_text_file.c_str(), "w");
        }

        fprintf(out, "CHROM\tPOS\tREF\tALT\tN_ALLELE");

        /////////////////////////
        //headers for filter tags
        /////////////////////////
        std::vector<std::string> filter_tag_str;
        for (uint32_t i=0; i<filter_tags.size(); ++i)
        {
            //filter tag is not examined to ensure it is
            //described in the header for VCF compatibility
            filter_tag_str.push_back(filter_tags[i]);
            fprintf(out, "\t%s", filter_tags[i].c_str());
        }

        if (info_tags.size()!=0) fprintf(out, "\t");

        //get the types of each info field
        std::vector<std::string> info_tag_str;
        std::vector<int32_t> info_tag_id;
        std::vector<int32_t> info_tag_vlen;
        std::vector<int32_t> info_tag_type;
        std::vector<int32_t> info_tag_num;

        ///////////////////////
        //headers for info tags
        ///////////////////////
        for (uint32_t i=0; i<info_tags.size(); ++i)
        {
            int32_t id = bcf_hdr_id2int(h, BCF_DT_ID, info_tags[i].c_str());

            if (id==-1)
            {
                notice("%s info tag does not exist", info_tags[i].c_str());
                continue;
            }

            int32_t vlen = bcf_hdr_id2length(h, BCF_HL_INFO, id);
            int32_t type = bcf_hdr_id2type(h, BCF_HL_INFO, id);
            int32_t num = bcf_hdr_id2number(h, BCF_HL_INFO, id);

            if (vlen==BCF_VL_FIXED)
            {
                info_tag_str.push_back(info_tags[i]);
                info_tag_id.push_back(id);
                info_tag_vlen.push_back(vlen);
                info_tag_type.push_back(type);
                info_tag_num.push_back(num);

                if (type==BCF_HT_FLAG)
                {
                    if (info_tag_str.size()>0) fprintf(out, "\t");
                    fprintf(out, "%s", info_tags[i].c_str());
                }
                else if (type==BCF_HT_INT)
                {
                    for (uint32_t j=0; j<num; ++j)
                    {
                        if (info_tag_str.size()!=1 || j!=0) fprintf(out, "\t");

                        if (num>1)
                        {
                            fprintf(out, "%s_%d", info_tags[i].c_str(), j+1);
                        }
                        else
                        {
                            fprintf(out, "%s", info_tags[i].c_str());
                        }
                    }
                }
                else if (type==BCF_HT_REAL)
                {
                    for (uint32_t j=0; j<num; ++j)
                    {
                        if (info_tag_str.size()!=1 || j!=0) fprintf(out, "\t");

                        if (num>1)
                        {
                            fprintf(out, "%s_%d", info_tags[i].c_str(), j+1);
                        }
                        else
                        {
                            fprintf(out, "%s", info_tags[i].c_str());
                        }
                    }
                }
            }
            else if (vlen==BCF_VL_VAR)
            {
                info_tag_str.push_back(info_tags[i]);
                info_tag_id.push_back(id);
                info_tag_vlen.push_back(vlen);
                info_tag_type.push_back(type);
                info_tag_num.push_back(num);

                if (info_tag_str.size()!=1) fprintf(out, "\t");

                if (type==BCF_HT_INT)
                {
                    fprintf(out, "%s", info_tags[i].c_str());
                }
                else if (type==BCF_HT_STR)
                {
                    fprintf(out, "%s", info_tags[i].c_str());
                }
            }
            else
            {
                notice("info2tab does not support id:%s, vlen=%s, type=%s, num=%d",
                               info_tags[i].c_str(), bcf_hdr_vl2str(vlen).c_str(), bcf_hdr_ht2str(type).c_str(), num);
            }
        }

        fprintf(out, "\n");

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

            //print variant information
            fprintf(out, "%s\t%d", bcf_get_chrom(h, v), bcf_get_pos1(v));
            char** alleles = bcf_get_allele(v);
            fprintf(out, "\t%s", alleles[0]);
            int32_t no_alleles = bcf_get_n_allele(v);
            for (uint32_t i=1; i<no_alleles; ++i)
            {
                fprintf(out, "%c", i==1?'\t':',');
                fprintf(out, "%s", alleles[i]);
            }
            fprintf(out, "\t%d\t", no_alleles);

//            bcf_print(h, v);

            for (uint32_t i=0; i<filter_tag_str.size(); ++i)
            {
                if (i)
                {
                    fprintf(out, "\t");
                }
/**
 *  Returns 1 if present, 0 if absent, or -1 if filter does not exist. "PASS" and "." can be used interchangeably.
 *   int bcf_has_filter(const bcf_hdr_t *hdr, bcf1_t *line, char *filter);
 */
                if (bcf_has_filter(h, v, const_cast<char*>(filter_tag_str[i].c_str()))==1)
                {
                    fprintf(out, "%d",  1);
                }
                else
                {
                    fprintf(out, "%d",  0);
                }
            }

            for (uint32_t i=0; i<info_tag_str.size(); ++i)
            {
                if (i || filter_tag_str.size()!=0)
                {
                    fprintf(out, "\t");
                }

                int32_t id = info_tag_id[i];
                int32_t vlen = info_tag_vlen[i];
                int32_t type = info_tag_type[i];
                int32_t num = info_tag_num[i];

//                notice("processing id:%s, vlen=%s, type=%s, num=%d",
//                               info_tag_str[i].c_str(), bcf_hdr_vl2str(vlen).c_str(), bcf_hdr_ht2str(type).c_str(), num);

                if (vlen==BCF_VL_FIXED)
                {
                    if (type==BCF_HT_FLAG)
                    {
                        bool present = bcf_get_info_flg(h, v, info_tag_str[i].c_str());
                        fprintf(out, "%d", (present ? 1 : 0));
                    }
                    else if (type==BCF_HT_INT)
                    {
                        for (uint32_t j=0; j<num; ++j)
                        {
                            if (j) fprintf(out, "\t");

                            int32_t val = bcf_get_info_int(h, v, info_tag_str[i].c_str());
                            fprintf(out, "%d", val);
                        }
                    }
                    else if (type==BCF_HT_REAL)
                    {
                        for (uint32_t j=0; j<num; ++j)
                        {
                            if (j) fprintf(out, "\t");

                            float val = bcf_get_info_flt(h, v, info_tag_str[i].c_str());
                            fprintf(out, "%.2f", val);
                        }
                    }
                    else if (type==BCF_HT_STR)
                    {
                        for (uint32_t j=0; j<num; ++j)
                        {
                            if (j) fprintf(out, "\t");

                            std::string val = bcf_get_info_str(h, v, info_tag_str[i].c_str());
                            fprintf(out, "%s", val.c_str());
                        }
                    }
                }
                else if (vlen == BCF_VL_VAR)
                {
                    //print as comma delimited in one column
                    if (type==BCF_HT_INT)
                    {
                        std::vector<int32_t> vals = bcf_get_info_int_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "\t");
                            fprintf(out, "%d", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_REAL)
                    {
                        std::vector<float> vals = bcf_get_info_flt_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, ",");
                            fprintf(out, "%.2f", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_STR)
                    {
                        std::vector<std::string> vals = bcf_get_info_str_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            if (j) fprintf(out, ",");
                            fprintf(out, "%s", vals[j].c_str());
                        }
                    }
                }
                else if (vlen==BCF_VL_G)
                {
                    //seems like the assumption here depends very much on ploidy.
                    //does not really seem reasonable to use in INFO field.
                    int32_t no_alleles = bcf_get_n_allele(v);
                    int32_t no_genotypes = no_alleles;

                    if (type==BCF_HT_INT)
                    {
                        std::vector<int32_t> vals = bcf_get_info_int_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%d", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_REAL)
                    {
                        std::vector<float> vals = bcf_get_info_flt_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%f", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_STR)
                    {
                        std::vector<std::string> vals = bcf_get_info_str_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%s", vals[j].c_str());
                        }
                    }
                }
                else if (vlen == BCF_VL_A)
                {
                    int32_t no_alt_alleles = bcf_get_n_allele(v) - 1;

                    if (type==BCF_HT_INT)
                    {
                        std::vector<int32_t> vals = bcf_get_info_int_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%d", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_REAL)
                    {
                        std::vector<float> vals = bcf_get_info_flt_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%f", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_STR)
                    {
                        std::vector<std::string> vals = bcf_get_info_str_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%s", vals[j].c_str());
                        }
                    }
                }
                else if (vlen == BCF_VL_R)
                {
                    int32_t no_alleles = bcf_get_n_allele(v);

                    if (type==BCF_HT_INT)
                    {
                        std::vector<int32_t> vals = bcf_get_info_int_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%d", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_REAL)
                    {
                        std::vector<float> vals = bcf_get_info_flt_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%f", vals[j]);
                        }
                    }
                    else if (type==BCF_HT_STR)
                    {
                        std::vector<std::string> vals = bcf_get_info_str_vec(h, v, info_tag_str[i].c_str());

                        for (uint32_t j=0; j<vals.size(); ++j)
                        {
                            fprintf(out, "%c", j==0?'\t':',');
                            fprintf(out, "%s", vals[j].c_str());
                        }
                    }
                }
            }

            fprintf(out, "\n");
            ++no_variants;
        }

        if (output_text_file!="-") fclose(out);
        odr->close();
    };

    void print_options()
    {
        std::clog << "info2tab v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output text file      " << output_text_file << "\n";
        print_strvec("         [u] filter tags           ", filter_tags);
        print_strvec("         [t] info tags             ", info_tags);
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
