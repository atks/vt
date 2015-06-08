/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#include "annotate_1000g.h"

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
    std::string l000g_vcf_file;
    std::vector<std::string> input_vcf_files;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
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
    int32_t no_variants;
    int32_t no_annotated_variants;

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
            std::string desc = "annotates variants that are present in 1000 Genomes variant set";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_1000g_vcf_file("d", "d", "1000G data set VCF file []", true, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "file", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            l000g_vcf_file = arg_1000g_vcf_file.getValue();
            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
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
        input_vcf_files.push_back(input_vcf_file);
        input_vcf_files.push_back(l000g_vcf_file);
        sr = new BCFSyncedReader(input_vcf_files, intervals, false);
        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(sr->hdrs[0]);
        bcf_hdr_append(sr->hdrs[0], "##INFO=<ID=1000G,Number=0,Type=Flag,Description=\"1000 Genomes variant\">");
        odw->write_hdr();

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp!="";

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip();

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_annotated_variants = 0;
    }

    void annotate_1000g()
    {
        std::vector<bcfptr*> crecs;
        int32_t presence[2] = {0,0};
        Variant variant;

        while(sr->read_next_position(crecs))
        {
            bcf1_t *v = NULL;
            bcf1_t *dv = NULL;

            //check existence
            for (size_t i=0; i<crecs.size(); ++i)
            {
                int32_t index = crecs[i]->file_index;

                if (filter_exists)
                {
                    if (!filter.apply(crecs[i]->h,crecs[i]->v,&variant))
                    {
                        continue;
                    }
                }

                if (index==0)
                {
                    v = crecs[i]->v;
                }

                if (index==1)
                {
                    dv = crecs[i]->v;
                }

                ++presence[index];
            }

            if (presence[0] && presence[1])
            {
                bcf_update_info_flag(odw->hdr, v, "1000G", "", 1);

                ++no_annotated_variants;
            }

            if (presence[0])
            {
                odw->write(v);

                ++no_variants;
            }

            presence[0] = 0;
            presence[1] = 0;
        }

        odw->close();
        sr->close();
    };

    void print_options()
    {
        std::clog << "annotate_1000g v" << version << "\n\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF file     " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file    " << output_vcf_file << "\n";
        print_str_op("         [f] filter             ", fexp);
        std::clog << "         [d] 1000G VCF file     " << l000g_vcf_file << "\n";
        print_int_op("         [i] intervals          ", intervals);
        std::clog << "\n";
   }

    void print_stats()
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "    no. variants             : %10d\n", no_variants);
        fprintf(stderr, "    no. annotated variants   : %10d\n", no_annotated_variants);
        fprintf(stderr, "\n");
    };

    ~Igor()
    {
    };

    private:
};

}

void annotate_1000g(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.annotate_1000g();
    igor.print_stats();
}
