/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

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

#include "profile_indels.h"

namespace
{
class BEDRecord: public Interval
{
    public:
    std::string chrom;

    BEDRecord(std::string& chrom, uint32_t start, uint32_t end)
    {
        this->chrom = chrom;
        this->start = start;
        this->end = end;
    };

    void print()
    {
        std::cerr << "chrom   : " << chrom << "\n";
        std::cerr << "[" << start << "," << end << "]\n";
    };

    private:
};

class GTFRecord: public Interval
{
    public:
    std::string gene;
    std::string feature;
    std::string chrom;
    char strand;
    int32_t frame;
    int32_t exonNo;
    bool fivePrimeConservedEssentialSpliceSite;
    bool threePrimeConservedEssentialSpliceSite;
    bool containsStartCodon;
    bool containsStopCodon;
    uint32_t level;
    std::string attrib;

    GTFRecord(std::string& chrom, uint32_t start, uint32_t end, char strand,
              std::string& gene, std::string& feature, int32_t frame, int32_t exonNo,
              bool fivePrimeConservedEssentialSpliceSite, bool threePrimeConservedEssentialSpliceSite,
              bool containsStartCodon, bool containsStopCodon,
              uint32_t level, std::string& attrib)
    {
        this->chrom = chrom;
        this->start = start;
        this->end = end;
        this->strand = strand;
        this->gene = gene;
        this->feature = feature;
        this->frame = frame;
        this->exonNo = exonNo;
        this->fivePrimeConservedEssentialSpliceSite = fivePrimeConservedEssentialSpliceSite;
        this->threePrimeConservedEssentialSpliceSite = threePrimeConservedEssentialSpliceSite;
        this->containsStartCodon = containsStartCodon;
        this->containsStopCodon = containsStopCodon;
        this->level = level;
        this->attrib = attrib;
    };

    void print()
    {
        std::cerr << "chrom   : " << chrom << "\n";
        std::cerr << "[" << start << "," << end << "]\n";
        std::cerr << "strand                    : " << strand << "\n";
        std::cerr << "address                   : " << this << "\n";
        std::cerr << "gene                      : " << gene << "\n";
        std::cerr << "feature                   : " << feature << "\n";
        std::cerr << "frame                     : " << frame << "\n";
        std::cerr << "exon number               : " << exonNo << "\n";
        std::cerr << "5' conserved splice site  : " << fivePrimeConservedEssentialSpliceSite << "\n";
        std::cerr << "3' conserved splice site  : " << threePrimeConservedEssentialSpliceSite << "\n";
        std::cerr << "contains start codon      : " << containsStartCodon << "\n";
        std::cerr << "contains stop codon       : " << containsStopCodon << "\n";
        std::cerr << "level                     : " << level << "\n";
        std::cerr << "attrib                    : " << attrib << "\n";
    };

    private:
};


class OverlapStats
{
    public:

    uint32_t a,ab,b,a_ins,ab_ins,b_ins,a_del,ab_del,b_del;

    OverlapStats()
    {
        a = 0;
        ab = 0;
        b = 0;

        a_ins = 0;
        a_del = 0;
        ab_ins = 0;
        ab_del = 0;
        b_ins = 0;
        b_del = 0;
    };
};

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string ref_fasta_file;
    std::string ref_data_sets_list;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    std::vector<OverlapStats> stats;
    std::string gencode_gtf_file;
    std::map<std::string, IntervalTree*> GENCODE;
    std::vector<Interval*> exons;
    
    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    bcf1_t *v;
    Filter *filter;
    kstring_t line;
    Filter *rare_filter;

    /////////
    //stats//
    /////////
    uint32_t no_candidate_snps;
    uint32_t no_candidate_indels;
    uint32_t nfs;
    uint32_t fs;
    uint32_t rare_nfs;
    uint32_t rare_fs;
    uint32_t common_nfs;
    uint32_t common_fs;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "profile Indels";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_ref_data_sets_list("g", "g", "file containing list of reference datasets []", false, "", "file", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            ref_fasta_file = arg_ref_fasta_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_data_sets_list = arg_ref_data_sets_list.getValue();
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
        //////////////////////
        //reference data set//
        //////////////////////        
        
        std::vector<std::string> input_vcf_files;
        input_vcf_files.push_back(input_vcf_file);
        
        std::vector<std::string> dataset_labels;
        dataset_labels.push_back("data");
        
        htsFile *hts = hts_open(ref_data_sets_list.c_str(), "r");
        kstring_t s = {0,0,0};
        std::vector<std::string> vec;
        while (hts_getline(hts, '\n', &s)>=0)
        {
            if (s.s[0] == '#')
                continue;    
            
            std::string line(s.s);
            split(vec, " ", line);
            
            //analysis purpose
            bool overlap = vec[1]=="overlap";
            
            //data set label
            dataset_labels.push_back(vec[0]);
            
            //filter
            
            //path
            if (overlap)
            {
                input_vcf_files.push_back(std::string(vec[3]));
            }
        }
        
        
//        //////////////////////
//        //i/o initialization//
//        //////////////////////
//        line = {0,0,0};
//
//        gencode_gtf_file = "/net/fantasia/home/atks/ref/encode/gencode.v15.annotation.gtf.gz";
//
//        //input_vcf_files.push_back(input_vcf_file);
//        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/1000g_phase1.indels.cplxsubs.sites.vcf.gz");
//        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/mills.doublehit.gatk.qc.indels.cplxsubs.sites.vcf.gz");
//        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/mills.chip.1000g_phase3.poly.indels.cplxsubs.sites.vcf.gz");
//        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/affymetrix.exome.chip.1000g_phase3.poly.indels.cplxsubs.sites.vcf.gz");
//        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/cg.mendel.errors.indels.cplxsubs.sites.vcf.gz");
//        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/affymetrix.exome.chip.mono.indels.cplxsubs.sites.vcf.gz");
//        input_vcf_files.push_back("/net/fantasia/home/atks/ref/all_chrom/cg.1000g_phase3.poly.indels.cplxsubs.sites.vcf.gz");
//
//        //input vcfs
//        sr = new BCFSyncedReader(input_vcf_files, intervals);
//
//        bcf_hdr_append(sr->hdrs[0], "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
//        bcf_hdr_append(sr->hdrs[0], "##INFO=<ID=FR,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
//
//        std::stringstream ss;
//        for (uint32_t i=0; i<intervals.size(); ++i)
//        {
//            ss.str("");
//            ss << "chr" << intervals[i].to_string();
//            populate_gencode_tree(ss.str().c_str());
//        }
        fs = 0;
        nfs = 0;
        rare_fs = 0;
        rare_nfs = 0;
        common_fs = 0;
        common_nfs = 0;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
        no_candidate_indels = 0;
    }

    void profile_indels()
    {
//        //map to contain the variant type
//        std::vector<bcf_hdr_t *> ivcf_hdrs;
//
//        //for combining the alleles
//        std::vector<bcfptr> current_recs;
//        std::map<std::string, std::map<std::string, int32_t> > variants;
//        std::stringstream ss;
//        std::stringstream centers;
//
//        dataset_labels = {"data", "1000g_phase1", "mills_doublehit", "mills_chip", "affymetrix", "cg_mendel", "affy_mono", "cg", "cg_rare" , "cg_common"};
//        stats.resize(dataset_labels.size());
//        std::vector<Interval*> intervals;
//
//        while(sr->read_next_position(current_recs))
//        {
//            variants.clear();
//            ss.str("");
//            centers.str("");
//
//            //for each file that contains the next record
//            char *ref = const_cast<char*>("N");
//            char *alt = const_cast<char*>("N");
//            bool in_ref = false;
//            for (uint32_t i=0; i<current_recs.size(); ++i)
//            {
//                int32_t d = current_recs[i].file_index;
//                bcf1_t *v = current_recs[i].v;
//                //bcf_set_variant_types(v);
//
//                if (d==0)
//                {
//                    bcf_hdr_t *h = sr->hdrs[0];
//                    if (!bcf_is_passed(h, v))
//                    {
////                        std::cerr << "FAIL PASS\n";
////                        vcf_format1(h, v, &line);
////                        std::cerr << line.s << "\n";
//                        continue;
//                    }
//
//                    if (filter!=NULL && !filter->apply(h, v))
//                    {
////                        std::cerr << "FAIL AF filter\n";
////                        vcf_format1(h, v, &line);
////                        std::cerr << line.s << "\n";
//                        continue;
//                    }
//                }
//
//                if (bcf_get_var_type(v)==VCF_SNP || bcf_get_var_type(v)==VCF_MNP || bcf_get_n_allele(v)!=2)
//                {
//                   continue;
//                }
//
//                if (d==0)
//                {
//                    ref = bcf_get_alt(v, 0);
//                    alt = bcf_get_alt(v, 1);
//
//                    int32_t pos1 = bcf_get_pos1(v);
//                    int32_t ref_len = strlen(ref);
//
//                    bcf_hdr_t *h = sr->hdrs[0];
//                    
//                    GENCODE[std::string(bcf_get_chrom(h, v))]->search(pos1+1, pos1+ref_len-1, exons);
//
////                  if ($start<=$F[$END]-1 && $end>=$F[$START])
////                  {
////                      my @alts = split(",", $alt);
////                      my $refLength = length($ref);
////                      my %diff = ();
////                      $diff{($refLength-length($_))%3} = 1 foreach @alts;
////
////                      if (exists($diff{1})||exists($diff{2}))
////                      {
////                          $GENCODE{"CDS_FS"} = 1;
////                      }
////                      else
////                      {
////                          $GENCODE{"CDS_NFS"} = 1;
////                      }
////                  }
//
//                    if (exons.size()!=0)
//                    {
//                        int32_t alt_len = strlen(alt);
//                        bool not_frame_shift = ((abs(ref_len-alt_len)%3)==0);
//                        if (not_frame_shift)
//                        {
//                            nfs += exons.size();
//                        }
//                        else
//                        {
//                            fs += exons.size();
//                        }
//
//                        if (rare_filter->apply(h, v))
//                        {
//                            if (not_frame_shift)
//                            {
//                                rare_nfs += exons.size();
//                            }
//                            else
//                            {
//                                rare_fs += exons.size();
//                            }
//                        }
//                        else
//                        {
//                            if (not_frame_shift)
//                            {
//                                common_nfs += exons.size();
//                            }
//                            else
//                            {
//                                common_fs += exons.size();
//                            }
//                        }
//                    }
//
//                    in_ref = true;
//                    int32_t ins = 0;
//                    int32_t del = 0;
//
//                    if (strlen(ref) > strlen(alt))
//                    {
//                        ++del;
//                    }
//                    else
//                    {
//                        ++ins;
//                    }
//
//                    for (uint32_t j=1; j<dataset_labels.size(); ++j)
//                    {
//                        ++stats[j].a;
//                        stats[j].a_ins += ins;
//                        stats[j].a_del += del;
//                    }
//                }
//                else
//                {
//                    update_stats(d, v, in_ref, ref, alt);
//
//                    if (d==7)
//                    {
//                        //rare
//                        bcf_hdr_t *h = sr->hdrs[7];
//                        if (rare_filter->apply(h, v))
//                        {
//                           update_stats(8, v, in_ref, ref, alt);
//                        }
//                        else
//                        {
//                           update_stats(9, v, in_ref, ref, alt);
//                        }
//                    }
//                }
//            }
//        }
    };

    void update_stats(int32_t d, bcf1_t *v, bool in_ref, char* ref, char* alt)
    {
        char* r = bcf_get_alt(v, 0);
        char* a = bcf_get_alt(v, 1);

        int32_t ins = 0;
        int32_t del = 0;

        if (strlen(r) > strlen(a))
        {
            ++del;
        }
        else
        {
            ++ins;
        }

        if (in_ref)
        {
            if (!strcmp(ref,r) && !strcmp(alt,a))
            {
                --stats[d].a;
                stats[d].a_ins -= ins;
                stats[d].a_del -= del;
                ++stats[d].ab;
                stats[d].ab_ins += ins;
                stats[d].ab_del += del;
            }
            else
            {
                ++stats[d].b;
                stats[d].b_ins += ins;
                stats[d].b_del += del;
            }
        }
        else
        {
            ++stats[d].b;
            stats[d].b_ins += ins;
            stats[d].b_del += del;
        }
    }

    void print_options()
    {
        std::clog << "profile_indels v" << version << "\n\n";
        std::clog << "\n";
        std::clog << "Options:     input VCF File        " << input_vcf_file << "\n";
        std::clog << "         [r] reference FASTA file  " << ref_fasta_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
   }

    void print_stats()
    {
//        for (uint32_t j=1; j<dataset_labels.size(); ++j)
//        {
//            double insdel_a_ab =  (stats[j].a_del+stats[j].ab_del)==0 ? -1 : ((double)(stats[j].a_ins+stats[j].ab_ins)/(double)(stats[j].a_del+stats[j].ab_del));
//            double insdel_a =  stats[j].a_del==0 ? -1 : ((double)stats[j].a_ins/(double)stats[j].a_del);
//            double insdel_ab =  stats[j].ab_del==0 ? -1 : ((double)stats[j].ab_ins/(double)stats[j].ab_del);
//            double insdel_b =  stats[j].b_del==0 ? -1 : ((double)stats[j].b_ins/(double)stats[j].b_del);
//            double insdel_b_ab =  (stats[j].b_del+stats[j].ab_del)==0 ? -1 : ((double)(stats[j].b_ins+stats[j].ab_ins)/(double)(stats[j].b_del+stats[j].ab_del));
//
//            uint32_t totala = stats[j].a+stats[j].ab;
//            double ab_of_a = totala==0 ? -1 : ((double)stats[j].ab/(double)(totala));
//
//            uint32_t totalb = stats[j].b+stats[j].ab;
//            double ab_of_b = totalb==0 ? -1 : ((double)stats[j].ab/(double)(totalb));
//
//            if (j==1)
//            {
//                std::cout << std::setprecision(3);
//                std::cout << "\n";
//                std::cout << dataset_labels[0] << "\n";
//                std::cout << "variants: " << totala << "\n";
//                std::cout << "ins/del ratio: " << insdel_a_ab << "\n";
//                std::cout << "FS Proportion: " << (float)fs/(float)(fs+nfs) << " (" << fs << "," << nfs << ")\n\n";
//                std::cout << "FS Proportion (Rare): " << (float)rare_fs/(float)(rare_fs+rare_nfs) << " (" << rare_fs << "," << rare_nfs << ")\n\n";
//                std::cout << "FS Proportion (Common): " << (float)common_fs/(float)(common_fs+common_nfs) << " (" << common_fs << "," << common_nfs << ")\n\n";
//
//            }
//
//            std::cout << std::setprecision(3);
//            std::cout << dataset_labels[j] << " (" << totalb << ") " << "[" << insdel_b_ab << "]\n";
//            std::cout << "TP\t" << ab_of_b << "(" << stats[j].ab << "/" << totalb << ") " << "[" << insdel_ab << "," << insdel_b << "] \n";
//            std::cout << "FP\t" << ab_of_a << "(" << stats[j].ab << "/" << totala << ") " << "[" << insdel_ab << "," << insdel_a << "] \n\n";
//        }
    };

    ~Igor()
    {
    };

    private:
    bool exists(std::map<std::string, IntervalTree*>& map, const std::string& key)
    {
        return map.end()!=map.find(key);
    }
};

}

void profile_indels(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.profile_indels();
    igor.print_stats();
}
