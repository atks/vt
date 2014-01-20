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

#include "concat.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_vcf_files;
    std::string input_vcf_file_list;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    bool print;
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    
    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Merge VCF. Includes only GT in format field by default.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", true, "", "str", cmd);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            
            cmd.parse(argc, argv);

            input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print = arg_print.getValue();

            ///////////////////////
            //parse input VCF files
            ///////////////////////
            htsFile *file = hts_open(input_vcf_file_list.c_str(), "r");
            if (file==NULL)
            {
                std::cerr << "cannot open " << input_vcf_file_list.c_str() << "\n";
                exit(1);
            }
            kstring_t *s = &file->line;
            while (hts_getline(file, KS_SEP_LINE, s) >= 0)
            {
                if (s->s[0]!='#')
                {
                    input_vcf_files.push_back(std::string(s->s));
                }
            }
            hts_close(file);
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
        odr = new BCFOrderedReader(input_vcf_files[0], intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        odw->link_hdr(odr->hdr);

        ///////////////
        //general use//
        ///////////////

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        /////////
        //tools//
        /////////
    }

    void concat()
    {
        odw->write_hdr();   
        
        for (int32_t i=1; i<input_vcf_files.size(); ++i)
        {
            if (i)
            {
                odr = new BCFOrderedReader(input_vcf_files[i], intervals);
            }
            
            bcf1_t *v = odw->get_bcf1_from_pool();
            while(odr->read(v))
            {
               odw->write(v);
               v =  odw->get_bcf1_from_pool();
            }
            
            odr->close();
        }
        
        odw->close();
    };

    void print_options()
    {
        if (!print) return;
        
        std::clog << "concat v" << version << "\n\n";
        std::clog << "options: [L] input VCF file list   " << input_vcf_file_list << " (" << input_vcf_files.size() << " files)\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;
        
        std::clog << "\n";
        std::cerr << "stats: Total Number of variants   " << no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {
    };

    private:
};

}

bool concat(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.concat();
    igor.print_stats();
    return igor.print;
}

