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

#include "compute_features.h"

#include <Rmath/Rmath.h>

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
    bool print_sites_only;
    
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
    int32_t no_variants;
    
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
            TCLAP::ValueArg<std::string> arg_interval_list("I", "interval-list", "File containing list of intervals", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF/VCF.GZ/BCF file [-]", false, "-", "str", cmd);
            TCLAP::SwitchArg arg_print_sites_only("s", "s", "print site information only without genotypes [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            print_sites_only = arg_print_sites_only.getValue();
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
        if (print_sites_only)
        {
            odw->link_hdr(bcf_hdr_subset(odr->hdr, 0, 0, 0));
        }
        else
        {
            odw->link_hdr(odr->hdr);
        }
        
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
        no_variants = 0;
        
        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");
    }

    void compute_features()
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

            

//            //update AC
//            bcf_unpack(v, BCF_UN_ALL);
//            int32_t ploidy = bcf_get_genotypes(odw->hdr, v, &gts, &n)/no_samples;        
//            int32_t n_allele = bcf_get_n_allele(v); 
//            
//            int32_t g[ploidy];
//            for (int32_t i=0; i<ploidy; ++i) g[i]=0;
//            int32_t AC[n_allele];
//            for (int32_t i=0; i<n_allele; ++i) AC[i]=0;
//            int32_t AN=0;
//            
//            for (int32_t i=0; i<no_samples; ++i)
//            {
//                for (int32_t j=0; j<ploidy; ++j)
//                {   
//                    g[j] = bcf_gt_allele(gts[i*ploidy+j]);
//                    
//                    if (g[j]>=0)
//                    {
//                        ++AC[g[j]];
//                        ++AN;
//                    }
//                }
//            }
//                
//            if (AC[0]<AN)
//            {   
//                int32_t* AC_PTR = &AC[1];
//                bcf_update_info_int32(odw->hdr,v,"VT_AC",AC_PTR,n_allele-1); 
//                bcf_update_info_int32(odw->hdr,v,"VT_AN",&AN,1);  
//            }
            
             std::cerr << "PVAL for 12 : " << pchisq(12, 1,0,0) << "\n";
            
            if (print_sites_only)
            {
                bcf_subset(odw->hdr, v, 0, 0);
            }

            odw->write(v);        
            ++no_variants;
        }

        odw->close();
    };

    void print_options()
    {
        std::clog << "compute_features v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File    " << input_vcf_file << "\n";
        print_str_op("         [f] filter            ", fexp);
        print_int_op("         [i] Intervals         ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: variants   : " << no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {

    };

    private:
};

}

void compute_features(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.compute_features();
    igor.print_stats();
}

