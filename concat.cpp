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
    bool print_sites_only;
    int32_t no_subset_samples;
    
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

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
    
    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Concatenate VCF files.  Assumes individuals are in the same order and files share the same header.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::SwitchArg arg_print_sites_only("s", "s", "print site information only without genotypes [false]", cmd, false);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);
                        
            cmd.parse(argc, argv);

            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            fexp = arg_fexp.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            no_subset_samples = arg_print_sites_only.getValue() ? 0 : -1;
            print = arg_print.getValue();
            
            if (input_vcf_files.size()==0)
            {
                fprintf(stderr, "[E:%s:%d %s] no input vcf files.\n", __FILE__, __LINE__, __FUNCTION__);
                exit(1);
            }
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
        if (no_subset_samples==-1)
        {
            odw->set_hdr(odr->hdr);
        }
        //perform subsetting
        else if (no_subset_samples==0)
        {
            odw->link_hdr(bcf_hdr_subset(odr->hdr, 0, 0, 0));
        }

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp=="" ? false : true;
      
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");    
    }

    void concat()
    {
        odw->write_hdr();   
        bcf1_t *v = odw->get_bcf1_from_pool();
        Variant variant;
        
        for (int32_t i=0; i<input_vcf_files.size(); ++i)
        {   
            std::cerr << "processing " << input_vcf_files[i] << "\n"; 
            
            if (i)
            {
                odr = new BCFOrderedReader(input_vcf_files[i], intervals);
            }
            
            bcf_hdr_t *h = odr->hdr;
            
            while(odr->read(v))
            {
                if (filter_exists)
                {
                    vm->classify_variant(h, v, variant);
                    if (!filter.apply(h, v, &variant))
                    {
                        continue;
                    }
                }
                
                if (no_subset_samples==0)
                {
                    bcf_subset(odw->hdr, v, 0, 0);
                    //maybe add some additional adhoc fixing for BCF files that do not have a complete header.
                
                }
         
                odw->write(v);
                ++no_variants;
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
        print_ifiles("options:     input VCF file        ", input_vcf_files);
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_str_op("         [f] filter                ", fexp);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;
        
        std::clog << "\n";
        std::cerr << "stats: no. of variants   " << no_variants << "\n";
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

