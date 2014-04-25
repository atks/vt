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
    Estimator *est;

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
        
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_AF,Number=A,Type=Float,Description=\"MLE Allele Frequency assuming HWE\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=VT_GF,Number=G,Type=Float,Description=\"MLE Genotype Frequency assuming HWE\">\n");
      
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
        est = new Estimator();
    }

    void compute_features()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odr->hdr;
        Variant variant;

        int32_t *gts = NULL;
        int32_t *pls = NULL;
        int32_t n_gts = 0;
        int32_t n_pls = 0;
        
        odw->write_hdr();    
            
        while(odr->read(v))
        {
            variant.clear();
            bool printed = false;

            ++no_variants;
            vm->classify_variant(h, v, variant);
            if (filter_exists && !filter.apply(h,v,&variant))
            {
                
                continue;
            }

            //update AC
            bcf_unpack(v, BCF_UN_ALL);
            int32_t ploidy = bcf_get_genotypes(odr->hdr, v, &gts, &n_gts);   
            ploidy /= no_samples;
            
            bcf_get_format_int32(odr->hdr, v, "PL", &pls, &n_pls);        
            int32_t n_allele = bcf_get_n_allele(v); 
            
            std::cerr << "ploidy " << ploidy << " no_samples " << no_samples << "\n";  
            
            int32_t g[ploidy];
            for (int32_t i=0; i<ploidy; ++i) g[i]=0;
            int32_t AC[n_allele];
            for (int32_t i=0; i<n_allele; ++i) AC[i]=0;
            int32_t AN=0;
            
            float MLE_HWE_AF[n_allele];
            int32_t no_genotypes = ((n_allele+1)*n_allele)>>1;
            float MLE_HWE_GF[no_genotypes];
            int32_t n = 0;
            est->compute_hwe_af(gts, pls, no_samples, ploidy,n_allele, MLE_HWE_AF, MLE_HWE_GF,  n, 1e-20);
            
            std::cerr << "no_allele    " << n_allele << " " << n << "\n";    
            std::cerr << "no_genotypes " << no_genotypes << " " << n << "\n";    
            if (n)
            {   
                float* MLE_HWE_AF_PTR = &MLE_HWE_AF[1];
                bcf_update_info_float(odw->hdr, v, "VT_AF", MLE_HWE_AF_PTR, n_allele-1); 
                bcf_update_info_float(odw->hdr, v, "VT_GF", &MLE_HWE_GF, no_genotypes);  
            }
            
    
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

