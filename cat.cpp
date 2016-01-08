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

#include "cat.h"

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
    int32_t compression_level;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    uint32_t sort_window_size;
    bool naive;
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
            std::string desc = "Concatenate VCF files. Individuals must be in the same order.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<int32_t> arg_compression_level("c", "c", "compression level 0-9, 0 and -1 denotes uncompressed with the former being wrapped in bgzf.[6]", false, 6, "int", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_sort_window_size("w", "w", "local sorting window size [0]", false, 0, "int", cmd);
            TCLAP::SwitchArg arg_naive("n", "n", "naive, assumes that headers are the same. [false]", cmd, false);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary [false]", cmd, false);
            TCLAP::SwitchArg arg_print_sites_only("s", "s", "print site information only without genotypes [false]", cmd, false);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);

            cmd.parse(argc, argv);

            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            compression_level = arg_compression_level.getValue();
            fexp = arg_fexp.getValue();
            sort_window_size = arg_sort_window_size.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            no_subset_samples = arg_print_sites_only.getValue() ? 0 : -1;
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
        //////////////////////
        //i/o initialization//
        //////////////////////
        if (input_vcf_files.size()==0)
        {
            fprintf(stderr, "[E:%s:%d %s] no input vcf files.\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
        odr = new BCFOrderedReader(input_vcf_files[0], intervals);
        odw = new BCFOrderedWriter(output_vcf_file, sort_window_size, compression_level);
        if (no_subset_samples==-1)
        {
            odw->set_hdr(odr->hdr);
        }
        else if (no_subset_samples==0)
        {
            odw->set_hdr(bcf_hdr_subset(odr->hdr, 0, 0, 0));
        }

        if (!naive)
        {
            if (input_vcf_files.size()>1)
            {
                uint32_t no_samples = bcf_hdr_get_n_sample(odw->hdr);
    
                //inspect headers in second file and above, insert missing headers
                for (uint32_t i=1; i<input_vcf_files.size(); ++i)
                {
                    vcfFile* file = bcf_open(input_vcf_files[i].c_str(), "r");
                    if (!file)
                    {
                        fprintf(stderr, "[%s:%d %s] Cannot open %s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files[i].c_str());
                        exit(1);
                    }
                    bcf_hdr_t* h = bcf_hdr_read(file);
    
                    //check samples consistency
                    if (no_subset_samples==-1)
                    {
                        if (no_samples!=bcf_hdr_get_n_sample(h))
                        {
                            fprintf(stderr, "[%s:%d %s] Different number of samples%s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files[i].c_str());
                            exit(1);
                        }
    
                        for (uint32_t j=0; j<no_samples; ++j)
                        {
                            if (strcmp(bcf_hdr_get_sample_name(odw->hdr, j), bcf_hdr_get_sample_name(h, j)))
                            {
                                fprintf(stderr, "[%s:%d %s] Different samples %s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files[i].c_str());
                                exit(1);
                            }
                        }
                    }
    
                    bcf_hdr_merge(odw->hdr, h);
    
                    bcf_close(file);
                    bcf_hdr_destroy(h);
                }
            }
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

    void cat()
    {
        if (!naive)
        {
            odw->write_hdr();
            bcf1_t *v = bcf_init();
            Variant variant;
    
            for (size_t i=0; i<input_vcf_files.size(); ++i)
            {
                if (i)
                {
                    delete odr;
                    odr = new BCFOrderedReader(input_vcf_files[i], intervals);
                }
    
                while(odr->read(v))
                {
                    if (filter_exists)
                    {
                        vm->classify_variant(odr->hdr, v, variant);
                        if (!filter.apply(odr->hdr, v, &variant))
                        {
                            continue;
                        }
                    }
    
                    if (no_subset_samples==0)
                    {
                        bcf_subset(odw->hdr, v, 0, 0);
                    }
    
                    bcf_translate(odw->hdr, odr->hdr, v);
                    odw->write(v);
                    if (sort_window_size)
                    {
                        v = odw->get_bcf1_from_pool();
                    }
                    ++no_variants;
                }
    
                odr->close();
            }
    
            odw->close();
        }
        else
        {
            odw->write_hdr();
            bcf1_t *v = bcf_init();
            Variant variant;
    
            for (size_t i=0; i<input_vcf_files.size(); ++i)
            {
                if (i)
                {
                    delete odr;
                    odr = new BCFOrderedReader(input_vcf_files[i], intervals);
                }
    
                while(odr->read(v))
                {
                    if (filter_exists)
                    {
                        vm->classify_variant(odr->hdr, v, variant);
                        if (!filter.apply(odr->hdr, v, &variant))
                        {
                            continue;
                        }
                    }
    
                    if (no_subset_samples==0)
                    {
                        bcf_subset(odw->hdr, v, 0, 0);
                    }
    
                    odw->write(v);
                    if (sort_window_size)
                    {
                        v = odw->get_bcf1_from_pool();
                    }
                    ++no_variants;
                }
    
                odr->close();
            }
    
            odw->close();
        }
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "cat v" << version << "\n\n";
        print_ifiles("options:     input VCF file        ", input_vcf_files);
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_num_op("         [w] sorting window size   ", sort_window_size);
        print_str_op("         [f] filter                ", fexp);
        std::clog << "         [s] print sites only      " << (no_subset_samples==-1?"no":"yes") << "\n";
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

bool cat(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.cat();
    igor.print_stats();
    return igor.print;
}