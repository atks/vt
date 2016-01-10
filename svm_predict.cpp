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

#include "svm_predict.h"

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
    std::vector<std::string> input_vcf_files;
    std::string ref_fasta_file;
    std::string ref_data_sets_list;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    bcf1_t *v;
    kstring_t line;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;
   
    /////////
    //stats//
    /////////
    int32_t no_snps;
    
    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;
   
    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "profile SNPs";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_model_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "VTYPE==SNP", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            input_vcf_file = arg_input_vcf_file.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip(ref_fasta_file);
        
        ////////////////////////
        //stats initialization//
        ////////////////////////

    }

    void svm_predict()
    {
        Variant variant;
        bcf1_t *v = bcf_init1();
        while(odr->read(v))
        {
            //apply model
            
            //add SVM score
         
            //write out
        }
    };

    void print_options()
    {
        std::clog << "svm_predict v" << version << "\n\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File                 " << input_vcf_file << "\n";
        print_str_op("         [f] filter                         ", fexp);
        print_int_op("         [i] intervals                      ", intervals);
        std::clog << "\n\n";
   }

    void print_stats()
    {
      

    };

    ~Igor()
    {
    };

    private:
};

}

void svm_predict(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.svm_predict();
    igor.print_stats();
}
