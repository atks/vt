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

#include "index.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    kstring_t output_vcf_index_file;
    bool print;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Indexes a VCF.GZ or BCF file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            print = arg_print.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
    }

    void index()
    {
        htsFile *file = hts_open(input_vcf_file.c_str(), "r");
        if (file==NULL)
        {
            exit(1);
        }    
        
        htsFormat ftype = file->format;
        if (ftype.compression!=bgzf&&ftype.format!=vcf&&ftype.format!=bcf)
        {
            fprintf(stderr, "[%s:%d %s] Not a BGZF VCF/BCF file: %s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
            exit(1);
        }

        int32_t min_shift;
        output_vcf_index_file = {0,0,0};
        int32_t ret;
        if (ftype.format==bcf)
        {
            kputs(input_vcf_file.c_str(), &output_vcf_index_file);
            kputs(".csi", &output_vcf_index_file);
            min_shift = 14;

            ret = bcf_index_build(input_vcf_file.c_str(), min_shift);
        }
        else if (ftype.format==vcf)
        {
            kputs(input_vcf_file.c_str(), &output_vcf_index_file);
            kputs(".tbi", &output_vcf_index_file);
            min_shift = 0;
            tbx_conf_t conf = tbx_conf_vcf;

            ret = tbx_index_build(input_vcf_file.c_str(), min_shift, &conf);
        }

        if (ret)
        {
            exit(1);
        }
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "index v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "             output index file     " << output_vcf_index_file.s << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
    };

    ~Igor() {};

    private:
};
}

bool index(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.index();
    igor.print_stats();
    return igor.print;
};