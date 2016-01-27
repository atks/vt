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

#include "profile_vntrs.h"

#define REF     0
#define TRUE    1

namespace
{

class VNTROverlapStats
{
    public:

    uint32_t a, ap, ab, bp, b;
    
    VNTROverlapStats()
    {
        a = 0;
        ap= 0;
        ab = 0;
        bp = 0;
        b = 0;
    };
};

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_vcf_files;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////////////////////
    //reference data sets//
    ///////////////////////
    std::string ref_data_sets_list;
    std::vector<std::string> dataset_labels;
    std::vector<int32_t> dataset_types;
    std::vector<std::string> dataset_fexps;
    std::vector<std::string> dataset_info_site_tags;
    std::vector<std::string> dataset_info_gt_tags;
    std::string cds_bed_file;
    std::string cplx_bed_file;
    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    std::vector<OrderedBCFOverlapMatcher *> oboms;
    OrderedRegionOverlapMatcher *orom_lcplx;
    OrderedRegionOverlapMatcher *orom_gencode_cds;

    //////////
    //filter//
    //////////
    std::string fexp;
    std::vector<Filter> filters;
    std::vector<bool> filter_exists;
    int32_t no_filters;

    /////////
    //stats//
    /////////
    std::vector<VNTROverlapStats> stats;
    uint32_t no_lcplx;
    uint32_t no_cds;
    uint32_t no_vntrs;
    uint32_t no_redundant;
                
    ////////////////
    //common tools//
    ////////////////
    VariantManip *vm;
    VNTRTree* vntr_tree;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "profile VNTRs.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_data_sets_list("g", "g", "file containing list of reference datasets []", false, "", "file", cmd);
            TCLAP::SwitchArg arg_write_partition("w", "w", "write partitioned variants to file", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            fexp = arg_fexp.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_data_sets_list = arg_ref_data_sets_list.getValue();
            input_vcf_files.push_back(arg_input_vcf_file.getValue());
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
        //reference data set//
        //////////////////////
//      # This file contains information on how to process reference data sets.
//      # dataset - name of data set, this label will be printed.
//      # type    - True Positives (TP) and False Positives (FP).
//      #           overlap percentages labeled as (Precision, Sensitivity) and (False Discovery Rate, Type I Error) respectively.
//      #         - annotation.
//      #           file is used for GENCODE annotation of frame shift and non frame shift Indels.
//      # filter  - filter applied to variants for this particular data set.
//      # path    - path of indexed BCF file.
//      #dataset     type            filter                       path
//      trf_lobstr   TP              VTYPE==VNTR                  /net/fantasia/home/atks/ref/vt/grch37/trf.lobstr.sites.bcf
//      GENCODE_V19  cds_annotation  .                            /net/fantasia/home/atks/ref/vt/grch37/gencode.cds.bed.gz
//      DUST         cplx_annotation .                            /net/fantasia/home/atks/ref/vt/grch37/mdust.bed.gz

        dataset_labels.push_back("data");
        dataset_types.push_back(REF);
        dataset_fexps.push_back(fexp);
        dataset_info_site_tags.push_back("");
        dataset_info_gt_tags.push_back("");

        htsFile *hts = hts_open(ref_data_sets_list.c_str(), "r");
        if (!hts)
        {
            fprintf(stderr, "[E:%s:%d %s] Reference file cannot be opened %s\n", __FILE__, __LINE__, __FUNCTION__, ref_data_sets_list.c_str());
            exit(1);
        }

        kstring_t s = {0,0,0};
        std::vector<std::string> vec;
        while (hts_getline(hts, '\n', &s)>=0)
        {
            if (s.s[0] == '#')
                continue;

            std::string line(s.s);
            split(vec, " ", line);

            if (vec[1] == "TP")
            {
                dataset_labels.push_back(vec[0]);
                dataset_types.push_back(TRUE);
                dataset_fexps.push_back(vec[2]);
                input_vcf_files.push_back(vec[3]);
            }
            else if (vec[1] == "cds_annotation")
            {
                cds_bed_file = vec[3];
            }
            else if (vec[1] == "cplx_annotation")
            {
                cplx_bed_file = vec[3];
            }
            else
            {
                fprintf(stderr, "[E:%s:%d %s] Reference data set type %s not recognized\n", __FILE__, __LINE__, __FUNCTION__, vec[1].c_str());
                exit(1);
            }
        }
        hts_close(hts);
        if (s.m) free(s.s);

        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_files[0], intervals);
        for (uint32_t i=0; i<input_vcf_files.size(); ++i)
        {
            oboms.push_back(new OrderedBCFOverlapMatcher(input_vcf_files[i], intervals));
        }

        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");
        vntr_tree = new VNTRTree();
        orom_gencode_cds = new OrderedRegionOverlapMatcher(cds_bed_file);
        orom_lcplx = new OrderedRegionOverlapMatcher(cplx_bed_file);

        /////////////////////////
        //filter initialization//
        /////////////////////////
        for (size_t i=0; i<dataset_fexps.size(); ++i)
        {
            filters.push_back(Filter(dataset_fexps[i]));
            filter_exists.push_back(dataset_fexps[i]!="");
        }
        no_filters = filters.size();

        ////////////////////////
        //stats initialization//
        ////////////////////////
        VNTROverlapStats vstats;
        stats.resize(input_vcf_files.size(), vstats);
        no_lcplx = 0;
        no_cds = 0;
        no_vntrs = 0;
        no_redundant = 0;
    }

    void profile_vntrs()
    {
        //for combining the alleles
        bcf1_t* v = bcf_init1();
        Variant variant;
        std::vector<int32_t> presence(2, 0);
        std::vector<bcf1_t*> overlap_vars;
        //for measure of redudancy        
        int32_t last_start1=0, last_end1=0, last_rid=0;
        
        while(odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            int32_t vtype = vm->classify_variant(odr->hdr, v, variant);

            if (vtype!=VT_VNTR)
            {
                continue;
            }

            ++no_vntrs;
            
            std::string chrom = bcf_get_chrom(odr->hdr,v);
            int32_t rid = bcf_get_rid(v);
            int32_t start1 = bcf_get_pos1(v);
            int32_t end1 = bcf_get_end_pos1(v);

            //last condition always true for an ordered file
            if (rid==last_rid && start1<=last_end1 && end1>=last_start1)
            {   
//                std::cerr << "REDUNDANT\n";
//                std::cerr << "\t" << rid << ":" << last_start1 << "-" << last_end1 << "\n";    
//                std::cerr << "\t" << rid << ":" << start1 << "-" << end1 << "\n";    
                
                ++no_redundant;
            
                last_rid = rid;
                last_start1 = std::min(last_start1, start1);
                last_end1 = std::max(last_end1, end1);    
            }
            else
            {
                last_rid = rid;
                last_start1 = start1;
                last_end1 = end1;
            }

            if (orom_lcplx->overlaps_with(chrom, start1, end1))
            {
                ++no_lcplx;
            }

            if (orom_gencode_cds->overlaps_with(chrom, start1, end1))
            {
                ++no_cds;
            }
            
            for (uint32_t i = 0; i<oboms.size(); ++i)
            {
                bool exact = false;
                if (oboms[i]->overlaps_with(chrom, start1, end1, overlap_vars))
                {
                    //check all overlapping variants
                    for (uint32_t j=0; j<overlap_vars.size(); ++j)
                    {
                        //check for exactness
                        
                        if (oboms[i]->is_exact_match(overlap_vars[j]))
                        {
                            exact = true;
                        }
                    }

                    if (exact)
                    {
                        ++stats[i].ab;
                    }
                    else
                    {
                        ++stats[i].ap;
                    }
                }
                else
                {
                    ++stats[i].a;
                }
            }

            vntr_tree->count(variant);
        }

        //for processing variants in the reference data sets that occur after the last record of the data set
        for (uint32_t i = 0; i<oboms.size(); ++i)
        {
            oboms[i]->flush();
            
            stats[i].bp += oboms[i]->get_no_fuzzy_overlaps();
            stats[i].b += oboms[i]->get_no_nonoverlaps();
        }

        odr->close();
    };

    void print_options()
    {
        std::clog << "profile_vntrs v" << version << "\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF file                 " << input_vcf_files[0] << "\n";
        print_str_op("         [f] filter                         ", fexp);
        std::clog << "         [g] reference data sets list file  " << ref_data_sets_list << "\n";
        print_int_op("         [i] intervals                      ", intervals);
        std::clog << "\n";
   }

    void print_stats()
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "    no VNTRs           %d\n", no_vntrs);
        fprintf(stderr, "    no low complexity  %d\n", no_lcplx);
        fprintf(stderr, "    no coding          %d\n", no_cds);
        fprintf(stderr, "    no redundant       %d\n", no_redundant);
        fprintf(stderr, "\n");
     
        for (int32_t i=0; i<dataset_labels.size(); ++i)
        {
            fprintf(stderr, "  %s\n", dataset_labels[i].c_str());
            fprintf(stderr, "    A-B  %10d \n", stats[i].a);
            fprintf(stderr, "    A-B~ %10d \n", stats[i].ap);
            fprintf(stderr, "    A&B  %10d \n", stats[i].ab);
            fprintf(stderr, "    B-A~ %10d \n", stats[i].bp);
            fprintf(stderr, "    B-A  %10d \n", stats[i].b);
            fprintf(stderr, "\n");
        }

        vntr_tree->print(BASIS);
                
        std::clog << "\n";
    };

    ~Igor()
    {
        for (uint32_t i = 0; i<oboms.size(); ++i)
        {
            delete oboms[i];
        }

        delete odr;

    };

    private:

};

}

void profile_vntrs(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_vntrs();
    igor.print_stats();
}
