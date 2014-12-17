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

#include "remove_overlap.h"

KHASH_MAP_INIT_STR(xdict, bcf1_t*);

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
    bool merge_by_pos;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    std::vector<bcf1_t*> pool;

    /////////
    //stats//
    /////////
    int32_t no_total_variants;
    int32_t no_nonoverlap_variants; 
    int32_t no_overlap_variants; 

    /////////
    //tools//
    /////////
    VariantManip *var_manip;

    Igor(int argc, char **argv)
    {
        version = "0.57";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Removes overlapping variants in a VCF file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
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
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        odw->link_hdr(odr->hdr);
        
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=PASS,Description=\"Passed variant\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap,Description=\"Overlapping variant\">");
        
        odw->write_hdr();
                
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_total_variants = 0;
        no_nonoverlap_variants = 0; 
        no_overlap_variants = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void remove_overlap()
    {
        kstring_t variant = {0, 0, 0};

        int32_t crid = -1;
        int32_t cpos1 = -1;
        int32_t cepos1 = -1;
        bcf1_t* cv = NULL;
        
        bcf1_t *v = odw->get_bcf1_from_pool();
        
        int32_t tpass_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "PASS");
        int32_t overlap_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap");
     
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            int32_t rid = bcf_get_rid(v);
            int32_t pos1 = bcf_get_pos1(v);
            int32_t epos1 = bcf_get_end_pos1(v);

            //does this overlap            
            if (crid==rid && cepos1>=pos1 && cpos1<=epos1)
            {
                //update overlap range
                cpos1 = std::min(cpos1, pos1);
                cepos1 = std::max(cepos1, epos1);
                
                if (cv)
                {
                    ++no_overlap_variants;
                    bcf_add_filter(odw->hdr, cv, overlap_id);
                    odw->write(cv);
                    cv = NULL;
                }    
                
                ++no_overlap_variants;       
            }
            else
            {
                crid = rid;
                cpos1 = pos1;
                cepos1 = epos1;
                
                if (cv)
                {    
                    bcf_add_filter(odw->hdr, cv, tpass_id);
                    odw->write(cv);
                    ++no_nonoverlap_variants;
                }
                cv = v;
                v = odw->get_bcf1_from_pool();
            }

            ++no_total_variants;
        }
        
        //the last non overlapping variant
        if (cv)
        {    
            bcf_add_filter(odw->hdr, cv, tpass_id);
            odw->write(cv);
            ++no_nonoverlap_variants;
        }
        
        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "remove_overlap v" << version << "\n\n";

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
        std::clog << "stats: Total number of observed variants    " << no_total_variants << "\n";
        std::clog << "       Total number of nonoverlap variants  " << no_nonoverlap_variants << "\n";
        std::clog << "       Total number of overlap variants     " << no_overlap_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void remove_overlap(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.remove_overlap();
    igor.print_stats();
};
