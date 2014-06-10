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

#include "genotype.h"
#include "lhmm_genotyping_record.h"
#include "chmm.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string sample_id;
    std::string input_vcf_file;
    std::string input_sam_file;
    std::string output_vcf_file;
    std::string ref_fasta_file;
    std::string model;
    std::vector<GenomeInterval> intervals;
    bool iterate_by_site;
    bool debug;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *vodr;
    BAMOrderedReader *sodr;
    BCFOrderedWriter *vodw;
    bcf1_t *v;
    bam1_t *s;

    /////////
    //stats//
    /////////
    int32_t no_snps_genotyped;
    int32_t no_indels_genotyped;

    /////////
    //tools//
    /////////
    faidx_t *fai;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Genotypes SNPs and Indels for a sample\n";

            version = "0.57";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_sam_file("b", "b", "input BAM file", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file", false, "-", "file", cmd);
            TCLAP::ValueArg<std::string> arg_sample_id("s", "s", "sample ID", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_model("m", "m", "model [lhmm]", true, "", "str", cmd);
            TCLAP::SwitchArg arg_iterate_by_site("c", "c", "iterate by candidate sites", cmd, false);
            TCLAP::SwitchArg arg_debug("d", "d", "debug alignments", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            input_sam_file = arg_input_sam_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
            sample_id = arg_sample_id.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            iterate_by_site = arg_iterate_by_site.getValue();
            debug = arg_debug.getValue();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    ~Igor()
    {
    };

    void print_options()
    {
        std::clog << "genotype v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [b] input BAM file        " << input_sam_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "         [c] iterate by site       " << (iterate_by_site ? "yes" : "no") << "\n";
        std::clog << "         [s] sample ID             " << sample_id << "\n";
        print_ref_op("         [r] reference fasta file  ", ref_fasta_file);
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void initialize()
    {
        if (model!="lhmm" && model!="chmm")
        {
            fprintf(stderr, "[E:%s] Probe information appears to be missing, cannot proceed unless reference FASTA file is available\n", __FUNCTION__);
            exit(1);
        }

        //////////////////////
        //i/o initialization//
        //////////////////////
        vodr = new BCFOrderedReader(input_vcf_file, intervals);
        sodr = new BAMOrderedReader(input_sam_file, intervals);
        vodw = new BCFOrderedWriter(output_vcf_file, 0);
        bcf_hdr_add_sample(vodw->hdr, sample_id.c_str());

        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=CT,Number=1,Type=String,Description=\"Count Genotype, number of repeat units as the genotype\">");
        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=QT,Number=1,Type=String,Description=\"Qualitative Genotype, e for exact sequence, i for inexact sequence\">");
        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=GQ,Number=1,Type=String,Description=\"Genotype Quality\">");
        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");
        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">");
        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allele Depth\">");
        bcf_hdr_append(vodw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED scaled genotype likelihood\">");

        //////////////////////
        //tool initialization//
        //////////////////////

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_snps_genotyped = 0;
        no_indels_genotyped = 0;
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "Stats: SNPs genotyped     " << no_snps_genotyped << "\n";
        std::clog << "       Indels genotyped   " << no_indels_genotyped << "\n";
        std::clog << "\n";
    }

    void genotype()
    {
        if (iterate_by_site)
        {
            GenotypingRecord *record;
            if (model=="lhmm")
            {
                record = new LHMMGenotypingRecord(fai);
            }
            else 
            {
            }
                            
            bcf1_t *v = vodw->get_bcf1_from_pool();
            bam1_t *s = bam_init1();
            int32_t no_read = 0;
            GenomeInterval interval;
            while (vodr->read(v))
            {
                record->set(v);
                std::string chrom(bcf_get_chrom(vodr->hdr, v));
                interval.set(chrom, bcf_get_pos1(v), bcf_get_pos1(v));
                sodr->jump_to_interval(interval);
                while(sodr->read(s))
                {
                    record->genotype(s);
                }

//                record.print(vodw);
                ++no_read;
            }
        }
        else //iterate by reads
        {
            //get reads from bam
                    //let VCFPool process reads

//            //pick up chromosomes from bam and vcf
//            //hash table for sequences
//            khash_t(sdict) *h = kh_init(sdict);
//            khiter_t k;
//            int32_t success;
//            if (intervals.empty())
//            {
//                kstring_t s = {0,0,0};
//                char** seqs = bam_hdr_get_target_name(sodr->hdr);
//                for (uint32_t i=0; i<bam_hdr_get_n_targets(sodr->hdr); ++i)
//                {
//                    if (kh_get(sdict, h, seqs[i])==kh_end(h))
//                    {
//                        char* seq = strdup(seqs[i]);
//                        k = kh_put(sdict, h, seq, &success);
//                        if (!success)
//                        {
//                            intervals.push_back(GenomeInterval(std::string(seq)));
//                        }
//                        else
//                        {
//                            free(seq);
//                        }
//                    }
//                }
//                free(seqs);
//
//                int32_t nseqs;
//                const char** seqs1 = bcf_hdr_seqnames(vodr->hdr, &nseqs);
//                for (uint32_t i=0; i<nseqs; ++i)
//                {
//                    if (kh_get(sdict, h, seqs1[i])==kh_end(h))
//                    {
//                        char* seq = strdup(seqs1[i]);
//                        k = kh_put(sdict, h, seq, &success);
//                        if (!success)
//                        {
//                            intervals.push_back(GenomeInterval(std::string(seq)));
//                        }
//                        else
//                        {
//                            free(seq);
//                        }
//                    }
//                }
//                free(seqs1);
//
//                //clean up
//                if (s.m) free(s.s);
//                for (k = kh_begin(h); k != kh_end(h); ++k)
//                {
//                    if (kh_exist(h, k))
//                    {
//                        free((char*)kh_key(h, k));
//                    }
//                }
//                kh_clear(sdict, h);
//            }
//
//            VCFGenotypingPool pool;
//
//            //iterate through chromosomes
//            for (uint32_t i=0; i<intervals.size(); ++i)
//            {
//                vodr->jump_to_interval(intervals[i]);
//                sodr->jump_to_interval(intervals[i]);
//
//                bcf1_t *v = vodw->get_bcf1_from_pool();
//                while (vodr->read(v))
//                {
//                    int32_t count = 0;
//                    while(sodr->read(s))
//                    {
//                        //update VCF records that overlap with read
//                        pool.process_read(s);
//                        ++count;
//                    }
//
//                    //remove records that occur before this read
//                    pool.flush();
//                }
//            }
//
//            vodw->close();
        }
    }

    private:

    void swap(double& a, double& b)
    {
        b = (a=a+b) - b;
        a -= b;
    }
};

}

void genotype(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.initialize();
    igor.print_options();
    igor.genotype();
    igor.print_stats();
}
