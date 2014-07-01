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

#include "config.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string resource_bundle_dir;

    ///////
    //i/o//
    ///////

    /////////
    //stats//
    /////////

    /////////
    //tools//
    /////////

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "configures the reference list in the vt bundle directory";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_resource_bundle_dir("r", "r", "resource bundle directory []", true, "", "str", cmd);
            //TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<c>", "input VCF file", true, "","file", cmd);
            cmd.parse(argc, argv);

            resource_bundle_dir = arg_resource_bundle_dir.getValue();
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


        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    void config()
    {
        //check directory?
        
        std::string output_reference_file = resource_bundle_dir + "/indel.reference.txt";

        htsFile *file = hts_open(output_reference_file.c_str(), "w");

        std::string hdr = "# This file contains information on how to process reference data sets.\n"
                          "# dataset - name of data set, this label will be printed.\n"
                          "# type    - True Positives (TP) and False Positives (FP).\n"
                          "#           overlap percentages labeled as (Precision, Sensitivity) and (False Discovery Rate, Type I Error) respectively.\n"
                          "#         - annotation.\n"
                          "#           file is used for GENCODE annotation of frame shift and non frame shift Indels.\n"
                          "# filter  - filter applied to variants for this particular data set.\n"
                          "# path    - path of indexed BCF file.\n";
                          
        size_t ret = hwrite(file->fp.hfile, hdr.c_str(), hdr.size());

        std::string dataset = "#dataset     type            filter                       path\n"
                              "1000g        TP              N_ALLELE==2&&VTYPE==INDEL    " + resource_bundle_dir + "1000G.snps_indels.sites.bcf\n"
                              "mills        TP              N_ALLELE==2&&VTYPE==INDEL    " + resource_bundle_dir + "mills.208620indels.sites.bcf\n"
                              "dbsnp        TP              N_ALLELE==2&&VTYPE==INDEL    " + resource_bundle_dir + "dbsnp.13147541variants.sites.bcf\n"
                              "GENCODE_V19  cds_annotation  .                            " + resource_bundle_dir + "gencode.cds.bed.gz\n"
                              "DUST         cplx_annotation .                            " + resource_bundle_dir + "mdust.bed.gz\n";


//#mills.chip.common      TP           N_ALLELE==2&&VTYPE==INDEL&&INFO.AF>0.005   mills.chip.158samples.8904indels.sites.bcf
//#affy.exome.chip        TP           N_ALLELE==2&&VTYPE==SNP                    affy.exome.chip.1249samples.316520variants.sites.bcf
//#affy.exome.chip.poly   TP           N_ALLELE==2&&VTYPE==INDEL&&INFO.AC!=0      affy.exome.chip.1249samples.316520variants.sites.bcf
//#affy.exome.chip.mono   FP           N_ALLELE==2&&VTYPE==INDEL&&INFO.AC==0      affy.exome.chip.1249samples.316520variants.sites.bcf

        ret = hwrite(file->fp.hfile, dataset.c_str(), dataset.size());

        std::cerr << "\twrote to " << output_reference_file << "\n";



        hts_close(file);
    };

    void print_options()
    {
        std::clog << "config v" << version << "\n";
        std::clog << "\n";
        std::clog << "options: [r] resource bundle directory  " << resource_bundle_dir << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {

    };

    ~Igor() {};

    private:
};

}

void config(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.config();
    igor.print_stats();
};
