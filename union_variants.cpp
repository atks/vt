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

#include "union_variants.h"

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
    std::vector<std::string> labels;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    
    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    ///////////////
    //general use//
    ///////////////
    kstring_t variant;
    
    /////////
    //stats//
    /////////
    int32_t no_unique_variants;
    int32_t no_variants;

    /////////
    //tools//
    /////////
    VariantManip * vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Union variants from a set of VCF files.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_labels("l", "l", "Comma delimited labels for the files", true, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);
            
            cmd.parse(argc, argv);

            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            split(labels, ",", arg_labels.getValue());
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            
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
        sr = new BCFSyncedReader(input_vcf_files, intervals, false);

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bcf_hdr_set_version(odw->hdr, "VCFv4.1");
        bcf_hdr_transfer_contigs(sr->hdrs[0], odw->hdr);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NCENTERS,Number=1,Type=Integer,Description=\"Number of centers with variant evidence.\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=CENTERS,Number=.,Type=String,Description=\"List of centers where variant is found.\">");
        odw->write_hdr();

        ///////////////
        //general use//
        ///////////////
        variant = {0,0,0};
        
        ////////////////
        //check labels//
        ////////////////
        if (labels.size()!=input_vcf_files.size())
        {
            fprintf(stderr, "[%s:%d %s] Number of labels has match the number of input files: %zu != %zu\n", __FILE__, __LINE__, __FUNCTION__, labels.size(), input_vcf_files.size());
            exit(1);
        }
        
        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_unique_variants = 0;
        no_variants = 0;
        
        /////////
        //tools//
        /////////
        vm = new VariantManip();
    }

    void union_variants()
    {
        int32_t nfiles = sr->get_nfiles();
        kstring_t centers = {0,0,0};
        int32_t ncenter = 0;
        int32_t *presence = new int32_t[nfiles];
        int32_t pass_filter = 0;
        
        std::vector<bcfptr*> current_recs;
        while(sr->read_next_position(current_recs))
        {
            ncenter = 0;
            centers.l = 0;
            for (size_t i=0; i<nfiles; ++i)
            {
                presence[i] = false;            
            }
                 
            for (size_t i=0; i<current_recs.size(); ++i)
            {
                int32_t file_index = current_recs[i]->file_index;    
                
                if (!presence[file_index])
                {
                    if (ncenter) kputc(',', &centers); 
                    kputs(labels[file_index].c_str(), &centers);     
                    presence[file_index] = true;
                    
                    ++no_variants;  
                    ++ncenter;            
                }
            }
                        
            bcf1_t *v = current_recs[0]->v;
            bcf_hdr_t *h = current_recs[0]->h;  
            
            //update variant information
            bcf1_t* nv = odw->get_bcf1_from_pool();
            bcf_set_chrom(odw->hdr, nv, bcf_get_chrom(h, v));
            bcf_set_pos1(nv, bcf_get_pos1(v));
            bcf_update_alleles(odw->hdr, nv, const_cast<const char**>(bcf_get_allele(v)), bcf_get_n_allele(v));
            bcf_update_info_string(odw->hdr, nv, "CENTERS", centers.s);
            bcf_update_info_int32(odw->hdr, nv, "NCENTERS", &ncenter, 1);
            bcf_update_filter(odw->hdr, nv, &pass_filter, 1);
            
            odw->write(nv);
            
            ++no_unique_variants;
        }

        sr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "union_variants v" << version << "\n\n";
        std::clog << "options: [L] input VCF file list   " << input_vcf_file_list << " (" << input_vcf_files.size() << " files)\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "stats: Total Number of variants                 " << no_variants << "\n";
        std::cerr << "       Total Number of unique variants          " << no_unique_variants << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {
    };

    private:
};

}

void union_variants(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.union_variants();
    igor.print_stats();
}

