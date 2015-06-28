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

#include "subset.h"

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
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    std::string arg_sample_list;
    char** samples;
    int32_t *imap;
    int32_t nsamples;

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
    int32_t no_samples;
    int32_t no_subset_samples;
    int32_t no_variants;
    int32_t no_subset_variants;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char ** argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Subsets a VCF file to a set of variants that are polymorphic on a selected set of individuals.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "intervals", "Intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "interval-list", "file containing list of intervals", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_sample_list("s", "s", "file containing list of samples []", true, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF/VCF.GZ/BCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            samples = read_sample_list(nsamples, arg_sample_list.getValue());
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
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
        odw = new BCFOrderedWriter(output_vcf_file);

        imap = (int32_t*) malloc(sizeof(int32_t)*nsamples);
        odw->link_hdr(bcf_hdr_subset(odr->hdr, nsamples, samples, imap));

        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n");

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_samples = bcf_hdr_nsamples(odr->hdr);
        no_subset_samples = bcf_hdr_nsamples(odw->hdr);
        no_variants = 0;
        no_subset_variants = 0;

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");
    }

    void subset()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odr->hdr;
        Variant variant;

        int32_t *gts = NULL;
        int32_t n = 0;

        odw->write_hdr();

        while(odr->read(v))
        {
            variant.clear();
            bool printed = false;

            ++no_variants;

            if (filter_exists && !filter.apply(h,v,&variant))
            {
                vm->classify_variant(h, v, variant);
                continue;
            }

            if (no_subset_samples)
            {
                bcf_subset(odw->hdr, v, nsamples, imap);

                //update AC
                bcf_unpack(v, BCF_UN_ALL);
                int32_t ploidy = bcf_get_genotypes(odw->hdr, v, &gts, &n)/no_subset_samples;
                int32_t n_allele = bcf_get_n_allele(v);

                int32_t g[ploidy];
                for (int32_t i=0; i<ploidy; ++i) g[i]=0;
                int32_t AC[n_allele];
                for (int32_t i=0; i<n_allele; ++i) AC[i]=0;
                int32_t AN=0;

                for (int32_t i=0; i<no_subset_samples; ++i)
                {
                    for (int32_t j=0; j<ploidy; ++j)
                    {
                        g[j] = bcf_gt_allele(gts[i*ploidy+j]);

                        if (g[j]>=0)
                        {
                            ++AC[g[j]];
                            ++AN;
                        }
                    }
                }

                if (AC[0]<AN)
                {
                    int32_t* AC_PTR = &AC[1];
                    bcf_update_info_int32(odw->hdr,v,"VT_AC",AC_PTR,n_allele-1);
                    bcf_update_info_int32(odw->hdr,v,"VT_AN",&AN,1);
                    odw->write(v);
                    ++no_subset_variants;
                }

                //check if the alleles used are a subset of the report alleles in shared data
                //reduce observed alleles

            }
        }

        odw->close();
    };

    void print_options()
    {
        std::clog << "subset v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File    " << input_vcf_file << "\n";
        std::clog << "         [s] sample file list  " << nsamples << " samples\n";
        print_str_op("         [f] filter            ", fexp);
        print_int_op("         [i] Intervals         ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: samples    : " << no_subset_samples << "/" << no_samples << "\n";
        std::clog << "       variants   : " << no_subset_variants << "/" << no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {

    };

    private:
};

}

void subset(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.subset();
    igor.print_stats();
}
