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

#include "set_ref.h"

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
    std::string ref_fasta_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    bool debug;
    std::string ref_region_tag;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    int32_t no_variants_with_REF_set;
    int32_t no_variants_within_bounds;
    int32_t no_variants_partial_overlap;
    int32_t no_variants_no_overlap;
    int32_t no_variants;

    ////////////////
    //common tools//
    ////////////////
    VariantManip* vm;
    ReferenceSequence* rs;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "annotates indels with VNTR information and adds a VNTR record.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_region_tag("t", "t", "reference region tag []", true, "", "str", cmd);
            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression []", false, "", "str", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            fexp = arg_fexp.getValue();
            debug = arg_debug.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
            ref_region_tag = arg_ref_region_tag.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
            abort();
        }
    };

    ~Igor()
    {
    };

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 10000);
        odw->link_hdr(odr->hdr);

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants_with_REF_set = 0;
        no_variants_within_bounds = 0;
        no_variants_partial_overlap = 0;
        no_variants_no_overlap = 0;
        no_variants = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip(ref_fasta_file);
        rs = new ReferenceSequence(ref_fasta_file);
    }

    void print_options()
    {
        std::clog << "set_ref v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file          " << output_vcf_file << "\n";
        std::clog << "         [t] reference region tag     " << ref_region_tag << "\n";
        print_boo_op("         [d] debug                    ", debug);
        print_str_op("         [f] filter                   ", fexp);
        print_ref_op("         [r] ref FASTA file           ", ref_fasta_file);
        print_int_op("         [i] intervals                ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: no. of variants with REF set   " << no_variants_with_REF_set << "\n";
        std::cerr << "           within bounds              " << no_variants_within_bounds << "\n";
        std::cerr << "           partial overlap            " << no_variants_partial_overlap << "\n";
        std::cerr << "           no overlap                 " << no_variants_no_overlap << "\n";
        std::clog << "\n";
        std::cerr << "       total no. of variants          " << no_variants << "\n";
        std::clog << "\n";
    }

    void set_ref()
    {
        odw->write_hdr();

        bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odw->hdr;
        Variant variant;

        while (odr->read(v))
        {
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);

            if (filter_exists)
            {
                if (!filter.apply(h, v, &variant, false))
                {
                    continue;
                }
            }

            if (vtype&VT_INDEL)
            {
                int32_t n = 0;
                int32_t* region1 = NULL;
                if (bcf_get_info_int32(h, v, ref_region_tag.c_str(), &region1, &n)>0)
                {
                    if (n==2)
                    {
                        int32_t pos1 = bcf_get_pos1(v);
                        int32_t end1 = bcf_get_end1(v);
                        
                        if (debug)
                        {    
                            std::cerr << "+++++++++++++++++++++++++++++++\n";
                            std::cerr << "region : [" << (region1[0]-1) << "," << (region1[1]+1) << "]\n";
                            std::cerr << "REF    : [" << pos1 << "," << end1 << "]\n";
                        }
                        --region1[0];
                        ++region1[1];
                        
                        if (pos1>=region1[0] && end1<=region1[1])
                        {
                            char* lflank = 0;
                            char* rflank = 0;
             
                            if (region1[0]<pos1) lflank = rs->fetch_seq(bcf_get_chrom(h,v), region1[0], pos1-1);
                            if (end1<region1[1]) rflank = rs->fetch_seq(bcf_get_chrom(h,v), end1+1, region1[1]);
                            kstring_t new_alleles = {0,0,0};
                            int32_t n_allele = bcf_get_n_allele(v);
                            char** alleles = bcf_get_allele(v);
                            
                            if (debug)
                            {    
                                std::cerr << "===============================\n";
                                std::cerr << "OLD alleles\n";
                                std::cerr << "pos1 : " << pos1 << "\n";
                                for (uint32_t i=0; i<n_allele; ++i)
                                {
                                    std::cerr << "       " << alleles[i] << "\n";
                                }
                            }
                            
                            if (debug)
                            {    
                                std::cerr << "add flanks\n";
                                if (lflank) std::cerr << "lflank: " << lflank << "\n";
                                if (rflank) std::cerr << "rflank: " << rflank << "\n";
                            }
                                       
                            for (uint32_t i=0; i<n_allele; ++i)
                            {
                                if (i) kputc(',', &new_alleles);
                                if (lflank) kputs(lflank, &new_alleles);
                                kputs(alleles[i], &new_alleles);
                                if (rflank) kputs(rflank, &new_alleles);
                            }
                            
                            if (lflank) free(lflank);
                            if (rflank) free(rflank);
                                
                            bcf_set_pos1(v, region1[0]);
                            bcf_update_alleles_str(h, v, new_alleles.s);
                                
                            if (debug)
                            {    
                                std::cerr << "NEW alleles\n";
                                std::cerr << "pos1 : " << region1[0] << "\n";
                                for (uint32_t i=0; i<n_allele; ++i)
                                {
                                    std::cerr << "       " << alleles[i] << "\n";
                                }
                            }
                            
                            ++no_variants_within_bounds;
                        }
                        else
                        {
                            if (end1>=region1[0] && pos1<=region1[1])
                            {
                                ++no_variants_partial_overlap;
                            }
                            else
                            {    
//                                std::cerr << "===============================\n";
//                                std::cerr << "weird\n";
//                                    
//                                bcf_print(h,v);  
                                
                                ++no_variants_no_overlap;
                            }
                        }
                    }
                }

                ++no_variants_with_REF_set;
            }

            ++no_variants;
            
            odw->write(v);
            v = odw->get_bcf1_from_pool();
        }

        odw->close();
        odr->close();
    };

    private:
};
}

void set_ref(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.set_ref();
    igor.print_stats();
};
