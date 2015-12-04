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

#include "merge_candidate_variants2.h"

namespace
{

class Evidence
{
    public:
    uint32_t i, m;
    uint32_t* e;
    uint32_t* n;
    kstring_t samples;
    int32_t esum, nsum;
    double af;
    double lr;
    bcf1_t *v;

    Evidence(uint32_t m)
    {
        this->m = m;
        i = 0;
        e = (uint32_t*) malloc(m*sizeof(uint32_t));
        n = (uint32_t*) malloc(m*sizeof(uint32_t));
        samples = {0,0,0};
        esum = 0;
        nsum = 0;
        af = 0;
        lr = 0;
        v = NULL;
    };

    ~Evidence()
    {
        i = 0;
        free(e);
        free(n);
        if (samples.m) free(samples.s);
        v = NULL;
    };

    void clear()
    {
        i = 0;
        samples.l = 0;
        esum = 0;
        nsum = 0;
        af = 0;
        lr = 0;
        v = NULL;
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
    std::string input_vcf_file_list;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;
    float snp_variant_score_cutoff;
    float indel_variant_score_cutoff;

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
    uint32_t no_samples;
    uint32_t no_candidate_snps;
    uint32_t no_candidate_indels;

    /////////
    //tools//
    /////////
    LogTool *lt;
    VariantManip * vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc =
"Merge candidate variants across samples.\n\
Each VCF file is required to have the FORMAT flags E and N and should have exactly one sample.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", true, "", "str", cmd);
            TCLAP::ValueArg<float> arg_snp_variant_score_cutoff("c", "c", "SNP variant score cutoff [30]", false, 30, "float", cmd);
            TCLAP::ValueArg<float> arg_indel_variant_score_cutoff("d", "d", "Indel variant score cutoff [30]", false, 30, "float", cmd);

            cmd.parse(argc, argv);

            input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            snp_variant_score_cutoff = arg_snp_variant_score_cutoff.getValue();
            indel_variant_score_cutoff = arg_indel_variant_score_cutoff.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());

            ///////////////////////
            //parse input VCF files
            ///////////////////////
            htsFile *file = hts_open(input_vcf_file_list.c_str(), "r");
            if (file==NULL)
            {
                std::cerr << "cannot open " << input_vcf_file_list.c_str() << "\n";
                exit(1);
            }
            kstring_t *s = &file->line;
            while (hts_getline(file, KS_SEP_LINE, s) >= 0)
            {
                if (s->s[0]!='#')
                {
                    input_vcf_files.push_back(std::string(s->s));
                }
            }
            hts_close(file);
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
        bcf_hdr_append(odw->hdr, "##fileformat=VCFv4.2");
        bcf_hdr_transfer_contigs(sr->hdrs[0], odw->hdr);
        bcf_hdr_append(odw->hdr, "##QUAL=Maximum variant score of the alternative allele likelihood ratio: -10 * log10 [P(Non variant)/P(Variant)] amongst all individuals.");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description=\"Number of samples.\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"Samples with evidence. (up to first 10 samples)\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=E,Number=.,Type=Integer,Description=\"Evidence read counts for each sample\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=N,Number=.,Type=Integer,Description=\"Read counts for each sample\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=ESUM,Number=1,Type=Integer,Description=\"Total evidence read count\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NSUM,Number=1,Type=Integer,Description=\"Total read count\">");

        odw->write_hdr();

        ///////////////
        //general use//
        ///////////////
        variant = {0,0,0};
        no_samples = sr->nfiles;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_candidate_snps = 0;
        no_candidate_indels = 0;
     
        /////////
        //tools//
        /////////
        lt = new LogTool();
        vm = new VariantManip();
    }

    void merge_candidate_variants2()
    {
        int32_t *E = (int32_t*) malloc(1*sizeof(int32_t));
        int32_t *N = (int32_t*) malloc(1*sizeof(int32_t));
        int32_t no_E = 1, no_N = 1;

        uint32_t nfiles = sr->get_nfiles();
        int32_t e[nfiles];
        int32_t n[nfiles];
        int32_t esum = 0;
        int32_t nsum = 0;
        kstring_t sample_names = {0,0,0};
        float af = 0;
        uint32_t no = 0;

        //obtain sample names
        char* index2sample[nfiles];
        for (uint32_t i=0; i<nfiles; ++i)
        {
            index2sample[i] = bcf_hdr_get_sample_name(sr->hdrs[i], 0);
        }

        bcf1_t* nv = bcf_init();
        Variant var;
        std::vector<bcfptr*> current_recs;
        while(sr->read_next_position(current_recs))
        {
            //aggregate statistics
            no = 0;
            esum = 0;
            nsum = 0;
            sample_names.l = 0;
            bool max_variant_score_gt_cutoff = false;
            float max_variant_score = 0;
            int32_t vtype;

            for (uint32_t i=0; i<current_recs.size(); ++i)
            {
                int32_t file_index = current_recs[i]->file_index;
                bcf1_t *v = current_recs[i]->v;
                bcf_hdr_t *h = current_recs[i]->h;

                if (i==0)
                {
                    //update variant information
                    bcf_clear(nv);
                    bcf_set_chrom(odw->hdr, nv, bcf_get_chrom(h, v));
                    bcf_set_pos1(nv, bcf_get_pos1(v));
                    bcf_update_alleles(odw->hdr, nv, const_cast<const char**>(bcf_get_allele(v)), bcf_get_n_allele(v));
                    vtype = vm->classify_variant(odw->hdr, nv, var);
                }

                if (bcf_get_format_int32(h, v, "E", &E, &no_E) < 0 ||
                    bcf_get_format_int32(h, v, "N", &N, &no_N) < 0)
                {
                    fprintf(stderr, "[E:%s:%d %s] cannot get format values E or N from %s\n", __FILE__, __LINE__, __FUNCTION__, sr->file_names[i].c_str());
                    exit(1);
                }

                float variant_score = bcf_get_qual(v);

                if (bcf_float_is_missing(variant_score))
                {
                    variant_score = 0;
                }

                if ((vtype == VT_SNP && variant_score >= snp_variant_score_cutoff) ||
                    (vtype == VT_INDEL && variant_score >= indel_variant_score_cutoff))
                {
                    max_variant_score_gt_cutoff = true;

                    if (max_variant_score < variant_score)
                    {
                        max_variant_score = variant_score;
                    }
                }
                
                e[i] = E[0];
                n[i] = N[0];
                esum += E[0];
                nsum += N[0];
                //just output first 10 samples
                if (no<=9)
                {    
                    if (i) kputc(',', &sample_names);
                    kputs(index2sample[file_index], &sample_names);
                }
                ++no;
            }

            if (max_variant_score_gt_cutoff)
            {
                bcf_update_info_int32(odw->hdr, nv, "NSAMPLES", &no, 1);
                bcf_update_info_string(odw->hdr, nv, "SAMPLES", sample_names.s);
                bcf_update_info_int32(odw->hdr, nv, "E", &e, no);
                bcf_update_info_int32(odw->hdr, nv, "N", &n, no);
                bcf_update_info_int32(odw->hdr, nv, "ESUM", &esum, 1);
                bcf_update_info_int32(odw->hdr, nv, "NSUM", &nsum, 1);
                bcf_set_qual(nv, max_variant_score);

                odw->write(nv);

                if (vtype == VT_SNP)
                {
                    ++no_candidate_snps;
                }
                else if (vtype == VT_INDEL)
                {
                    ++no_candidate_indels;
                }
            }
        }

        sr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "merge_candidate_variants v" << version << "\n\n";
        std::clog << "options: [L] input VCF file list         " << input_vcf_file_list << " (" << input_vcf_files.size() << " files)\n";
        std::clog << "         [o] output VCF file             " << output_vcf_file << "\n";
        std::clog << "         [c] SNP variant score cutoff    " << snp_variant_score_cutoff << "\n";
        std::clog << "         [d] Indel variant score cutoff  " << indel_variant_score_cutoff << "\n";
        print_int_op("         [i] intervals                   ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: Total Number of Candidate SNPs                 " << no_candidate_snps << "\n";
        std::clog << "       Total Number of Candidate Indels               " << no_candidate_indels << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {
    };

    private:
};

}

void merge_candidate_variants2(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.merge_candidate_variants2();
    igor.print_stats();
}

