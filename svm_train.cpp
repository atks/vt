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

#include "svm_train.h"

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
    std::string positive_training_set_vcf_file;
    std::string negative_training_set_vcf_file;
    std::string model;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;

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
    int32_t no_variants;

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
            std::string desc = "svm train";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_positive_training_set_vcf_file("p", "p", "file containing list of reference datasets []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_negative_training_set_vcf_file("p", "p", "file containing list of reference datasets []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "VTYPE==SNP", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            positive_training_set_vcf_file = arg_positive_training_set_vcf_file.getValue();
            negative_training_set_vcf_file = arg_negative_training_set_vcf_file.getValue();
            input_vcf_file = arg_input_vcf_file.getValue();

            ///////////////////////
            //parse input VCF files
            ///////////////////////
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
        std::vector<std::string> input_vcf_files;
        sr = new BCFSyncedReader(input_vcf_files, intervals, false);

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_snps = 0;

    }

    void svm_train()
    {
        //extract common variants with training sets
           //store covariates
        std::vector<bcfptr*> current_recs;
        std::vector<Interval*> overlaps;
        Variant variant;

        while(sr->read_next_position(current_recs))
        {
            //check first variant
            bcf1_t *v = current_recs[0]->v;
            bcf_hdr_t *h = current_recs[0]->h;
            int32_t vtype = vm->classify_variant(h, v, variant);

        }

        //inverse normalize covariates

        //train

        //output model trained

    };

    void print_options()
    {
        std::clog << "svm_train v" << version << "\n\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File                 " << input_vcf_file << "\n";
        print_str_op("         [f] filter                         ", fexp);
        print_int_op("         [i] intervals                      ", intervals);
        std::clog << "\n\n";
   }

    void print_stats()
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "  %s\n", "data set");
        fprintf(stderr, "\n");

    };

    ~Igor()
    {
    };

    private:
};

}

void svm_train(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.svm_train();
    igor.print_stats();
}
