/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#include "uniq.h"

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
    
    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    BCFOrderedWriter *odw;

    std::vector<bcf1_t*> pool;

    /////////
    //stats//
    /////////
    uint32_t no_total_variants;
    uint32_t no_unique_variants;

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
            std::string desc = "Drops duplicate variants that appear later in the the VCF file.";

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
        std::vector<std::string> input_vcf_files;
        input_vcf_files.push_back(input_vcf_file);  
        sr = new BCFSyncedReader(input_vcf_files, intervals, SYNC_BY_VAR);
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        odw->link_hdr(sr->hdrs[0]);
        odw->write_hdr();

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_total_variants = 0;
        no_unique_variants = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void uniq()
    {
        std::vector<bcfptr*> crecs;
        while (sr->read_next_position(crecs))
        {
            if (crecs.size()>1)
            {
                //check if the first variant has OLD_VARIANT
                //proceed to merge the others.
                char* dst = 0;
                int32_t ndst = 0;    
                std::map<std::string, uint32_t> ov_map;
                std::string old_vars;
                for (uint32_t i=0; i<crecs.size(); ++i)
                {
                    if(bcf_get_info_string(sr->hdrs[0], crecs[i]->v, "OLD_VARIANT", &dst, &ndst)>0)
                    {
                        std::string old_var(dst);
                        if (ov_map.find(dst)==ov_map.end())
                        {
                            ov_map[old_var] = 1;
                            if (old_vars!="") old_vars.append(1, ',');
                            old_vars.append(old_var);
                        }
                    }
                }
                if (dst!=0) free(dst);
                
                //update OLD_VARIANT in the first record.
                if (old_vars != "")
                {
                    bcf_update_info_string(sr->hdrs[0], crecs[0]->v, "OLD_VARIANT", old_vars.c_str());
                }   
            }
            
            odw->write(crecs[0]->v);
            
            ++no_unique_variants;
            no_total_variants += crecs.size();
        }

        sr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "uniq v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
    }

    void print_stats()
    {
        std::clog << "\nstats: Total number of observed variants   " << no_total_variants << "\n";
        std::clog <<   "       Total number of unique variants     " << no_unique_variants << "\n\n";
    };

    ~Igor() {};

    private:
};

}

void uniq(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.uniq();
    igor.print_stats();
};
