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

#include "compute_concordance.h"

namespace
{

#define RR 0;
#define RA 1;
#define AA 2;
#define NA 3;

class VariantConcordanceStats
{
    public:

    int32_t concordance[4][4];

    VariantConcordanceStats()
    {
        for (uint32_t i=0; i<4; ++i)
        {
            for (uint32_t j=0; j<4; ++j)
            {
                concordance[i][j] = 0;
            }
        }
    };

    void update_stats(bcf_hdr_t *h1, bcf1_t *v1, bcf_hdr_t *h2, bcf1_t *v2)
    {

    };

};

class ConcordanceStats
{
    public:

    int32_t concordance[4][4];

    ConcordanceStats()
    {
        for (uint32_t i=0; i<4; ++i)
        {
            for (uint32_t j=0; j<4; ++j)
            {
                concordance[i][j] = 0;
            }
        }
    };

    void reset()
    {
        for (uint32_t i=0; i<4; ++i)
        {
            for (uint32_t j=0; j<4; ++j)
            {
                concordance[i][j] = 0;
            }
        }
    }
};

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_vcf_files;
    std::string filters;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    std::string variant_concordance_txt_file;
    std::string sample_concordance_txt_file;

    ///////
    //i/o//
    ///////
    htsFile *variant_concordance_txt;
    htsFile *sample_concordance_txt;

    BCFSyncedReader *sr;
    bcf1_t *v;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;

    /////////
    //stats//
    /////////
    std::vector<ConcordanceStats> stats;
    uint32_t no_candidate_snps;
    uint32_t no_candidate_indels;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    
    void intersect_samples(bcf_hdr_t* h1, bcf_hdr_t *h2, std::vector<std::string>& s, std::vector<int32_t>& a, std::vector<int32_t>& b)
    {
        vdict_t *d2 = (vdict_t*)h2->dict[BCF_DT_SAMPLE];

        int32_t nsamples = bcf_hdr_get_n_sample(sr->hdrs[0]);
        char** samples = bcf_hdr_get_samples(sr->hdrs[0]);

        std::map<std::string, int32_t> m;

        for (int32_t i=0; i<nsamples; ++i)
        {
            khint_t k = kh_get(vdict, d2, samples[i]);
            if ( k != kh_end(d2) )
            {
                s.push_back(samples[i]);
                a.push_back(i);
                b.push_back(kh_val(d2, k).id);
            }
        }
    }

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Compute Concordance.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "Intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "File containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_variant_concordance_txt_file("m", "m", "Variant concordance text file [m.txt]", false, "m.txt", "str", cmd);
            TCLAP::ValueArg<std::string> arg_sample_concordance_txt_file("s", "s", "Sample concordance text file [s.txt]", false, "s.txt", "str", cmd);
            TCLAP::ValueArg<std::string> arg_filters("f", "f", "filter expression", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf><in2.vcf>", "Input VCF File", true, "str", cmd);

            cmd.parse(argc, argv);

            if (arg_input_vcf_files.getValue().size()!=2)
            {
                fprintf(stderr, "[%s:%d %s] Requires 2 VCF files as input\n", __FILE__, __LINE__, __FUNCTION__);
                exit(1);
            }

            input_vcf_files = arg_input_vcf_files.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            filters = arg_filters.getValue();
            variant_concordance_txt_file = arg_variant_concordance_txt_file.getValue();
            sample_concordance_txt_file = arg_sample_concordance_txt_file.getValue();
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
        sr = new BCFSyncedReader(input_vcf_files, intervals);
        variant_concordance_txt = hts_open(variant_concordance_txt_file.c_str(), "w");
        sample_concordance_txt = hts_open(sample_concordance_txt_file.c_str(), "w");

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str());
        filter_exists = fexp=="" ? false : true;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
        no_candidate_indels = 0;
        
        ///////////////////////
        //tool initialization//
        ///////////////////////
        vm = new VariantManip("");
    }

    void compute_concordance()
    {
        ///////////////////////
        //find common samples//
        ///////////////////////
        char** samples1 = bcf_hdr_get_samples(sr->hdrs[0]);
        char** samples2 = bcf_hdr_get_samples(sr->hdrs[1]);

        std::vector<std::string> s;
        std::vector<int32_t> a;
        std::vector<int32_t> b;

        intersect_samples(sr->hdrs[0], sr->hdrs[1], s, a, b);

        stats.resize(s.size());

        kstring_t *line = &variant_concordance_txt->line;

        ConcordanceStats variant_stat;
        kputs("variant\tRR_RR\tRR_RA\tRR_AA\tRR_NA\tRA_RR\tRA_RA\tRA_AA\tRA_NA\tAA_RR\tAA_RA\tAA_AA\tAA_NA\tNA_RR\tNA_RA\tNA_AA\tNA_NA\n", line);
        //std::cerr << line.s;

        //for combining the alleles
        std::vector<bcfptr*> current_recs;
        std::map<std::string, bcfptr*> variants;
        
        Variant variant;

        while(sr->read_next_position(current_recs))
        {
            if (current_recs.size()==2 && current_recs[0]->file_index!=current_recs[1]->file_index)
            {
                if (bcf_get_n_allele(v)!=2)
                {
                    continue;
                }
                
                if (filter_exists)
                {     
                    vm->classify_variant(current_recs[0]->h, current_recs[0]->v, variant);
                    if (!filter.apply(current_recs[0]->h, current_recs[0]->v, &variant))
                    {    
                        continue;
                    }
                }

                bcf1_t *v1 = current_recs[0]->v;
                bcf_hdr_t *h1 = sr->hdrs[0];
                bcf_fmt_t* f1 = bcf_get_fmt(h1, v1, "GT");

                bcf1_t *v2 = current_recs[1]->v;
                bcf_hdr_t *h2 = sr->hdrs[1];
                bcf_fmt_t* f2 = bcf_get_fmt(h2, v2, "GT");

                for (size_t i=0; i<s.size(); ++i)
                {
                    int8_t *x1 = (int8_t*)(f1->p + a[i] * f1->size);
                    int8_t *x2 = (int8_t*)(f2->p + b[i] * f2->size);

                    int32_t g1 = 0;
                    int32_t g2 = 0;
                    int32_t l = 0;

                    for (l = 0; l < f1->n && x1[l]!=INT8_MIN; ++l)
                    {
                        if (x1[l] == 0)
                        {
                            g1 = 3;
                            break;
                        }
                        else if (x1[l]>>1)
                        {
                            g1 += (x1[l]>>1) - 1;
                        }
                    }

                    for (l = 0; l < f2->n && x1[l]!=INT8_MIN; ++l)
                    {
                        if (x2[l] == 0)
                        {
                            g2 = 3;
                            break;
                        }
                        else if (x2[l]>>1)
                        {
                            g2 += (x2[l]>>1) - 1;
                        }
                    }

                    ++stats[i].concordance[g1][g2];
                    ++variant_stat.concordance[g1][g2];
                }

//                line->l = 0;
//                kputs(sr->get_current_sequence().c_str(), line);
//                kputc(':', line);
//                kputw(sr->get_current_pos1(), line);
//                kputc(':', line);
//                kputs(ss.str().c_str(), line);
//                kputc('\t', line);
//
//                for (uint32_t j=0; j<4; ++j)
//                {
//                    for (uint32_t k=0; k<4; ++k)
//                    {
//                        kputw(variant_stat.concordance[j][k], line);
//                        kputs("\t", line);
//
//                    }
//                }
//
//                kputs("\n", line);
            }
        }

        /////////////////////////////////////////////
        //output concordance matrix for individuals//
        /////////////////////////////////////////////
        line = &sample_concordance_txt->line;
        kputs("sample\tRR_RR\tRR_RA\tRR_AA\tRR_NA\tRA_RR\tRA_RA\tRA_AA\tRA_NA\tAA_RR\tAA_RA\tAA_AA\tAA_NA\tNA_RR\tNA_RA\tNA_AA\tNA_NA\n", line);
        std::cerr << line->s;
        for (uint32_t i=0; i<s.size(); ++i)
        {
            line->l = 0;
            kputs(s[i].c_str(), line);
            kputc('\t', line);
            for (uint32_t j=0; j<4; ++j)
            {
                for (uint32_t k=0; k<4; ++k)
                {
                    kputw(stats[i].concordance[j][k], line);
                    kputs("\t", line);
                }
            }

            kputs("\n", line);
            std::cerr << line->s;
        }
      
        hts_close(variant_concordance_txt);
        hts_close(sample_concordance_txt);
    };

    void print_options()
    {
        std::clog << "compute_concordance v" << version << "\n\n";

        std::clog << "Options:     Input VCF File A     " << input_vcf_files[0] << "\n";
        std::clog << "             Input VCF File B     " << input_vcf_files[1] << "\n";
        std::clog << "             sample concordance   " << variant_concordance_txt_file << "\n";
        std::clog << "             variant concordance  " << sample_concordance_txt_file << "\n";
        
        print_str_op("         [f] filter               ", fexp);
        print_int_op("         [i] intervals            ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {


    };

    ~Igor()
    {
    };

    private:
};

}

void compute_concordance(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.compute_concordance();
    igor.print_stats();
}
