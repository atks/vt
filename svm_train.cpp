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
            TCLAP::ValueArg<std::string> arg_positive_training_set_vcf_file("p", "p", "positive training data set []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_negative_training_set_vcf_file("n", "n", "negative training data set []", false, "", "file", cmd);
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
        no_variants = 0;
    }

    /**
     * Inverse normalizes a set of values.
     */
    void inverse_normalize(std::vector<float>& raw, std::vector<float>& normalized)
    {
        //sort vector
        std::vector<float*> raw_ptr(raw.size());
        for (int32_t i=0; i<raw.size(); ++i)
        {
            raw_ptr[i] = &raw[i];
        }

        struct
        {
            bool operator()(float *a, float *b)
            {
                return *a < *b;
            }
        } ptr_compare;

        //print unsorted values
        std::cout << "Unsorted Values\n";
        for (int32_t i=1; i<raw_ptr.size(); ++i)
        {
            std::cerr << *raw_ptr[i] << "\n";
        }

        std::sort(raw_ptr.begin(), raw_ptr.end(), ptr_compare);

        //print sorted values
        std::cout << "Sorted Values\n";
        for (int32_t i=1; i<raw_ptr.size(); ++i)
        {
            std::cerr << *raw_ptr[i] << "\n";
        }

        //count unique values
        int32_t no_uniq_values = raw_ptr.size()==0 ? 0 : 1;
        for (int32_t i=1; i<raw.size(); ++i)
        {
            if (*raw_ptr[i]!=*raw_ptr[i-1])
            {
                ++no_uniq_values;
            }
        }

        std::cout << "No of unique values : " << no_uniq_values << "\n";

        //inverse normalize values
        //update values
        float* offset = &raw[0];
        float step = 0.9999999/(float)no_uniq_values;
        float c_p = step;
        normalized.resize(raw.size());
        int32_t pos = (raw_ptr[0]-offset);
        std::cerr << "updating position " << pos << "\n";
        normalized[pos] = qnorm(c_p, 0, 1, 1, 0);
        for (int32_t i=1; i<raw_ptr.size(); ++i)
        {
            if (*raw_ptr[i]!=*raw_ptr[i-1])
            {
                c_p += step;
            }

            pos = (raw_ptr[i]-offset);
            std::cerr << "updating position " << pos << "\n";

            normalized[pos] = qnorm(c_p, 0, 1, 1, 0);
        }


    }

    void svm_train()
    {


        //extract common variants with training sets
           //store covariates
//        std::vector<bcfptr*> current_recs;
//        std::vector<Interval*> overlaps;
//        Variant variant;
//
//        while(sr->read_next_position(current_recs))
//        {
//            //check first variant
//            bcf1_t *v = current_recs[0]->v;
//            bcf_hdr_t *h = current_recs[0]->h;
//            int32_t vtype = vm->classify_variant(h, v, variant);
//
//        }

        //inverse normalize covariates in memory
        std::vector<float> raw;
        std::vector<float> normalized;

        for (int32_t i=0; i<10; ++i)
        {
            raw.push_back(runif(-100, 100));
        }

        inverse_normalize(raw, normalized);

        for (int32_t i=0; i<10; ++i)
        {
            printf("%.4f\t%.4f\n", raw[i], normalized[i]);
        }


        //train



        //output model trained


        //output inverse normalized features to allow for inspection

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
