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

#include "discover2.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::vector<GenomeInterval> intervals;
    std::string output_vcf_file;
    std::string input_bam_file;
    std::string ref_fasta_file;
    std::string sample_id;
    uint32_t mapq_cutoff;
    uint32_t baseq_cutoff;
    //takes on snps, mnps, indels
    std::string variant_type;
    uint32_t evidence_allele_count_cutoff;
    double fractional_evidence_allele_count_cutoff;

    uint16_t exclude_flag;

    ///////
    //i/o//
    ///////
    BAMOrderedReader *odr;
    bam1_t *s;

    BCFOrderedWriter *odw;
    bcf1_t *v;

    /////////
    //stats//
    /////////
    uint32_t no_reads;
    uint32_t no_overlapping_reads;
    uint32_t no_passed_reads;
    uint32_t no_exclude_flag_reads;
    uint32_t no_low_mapq_reads;

    /////////
    //tools//
    /////////
    VariantHunter *variantHunter;
    BCFOrderedWriter *odw;
    
    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Discovers variants from reads in a BAM file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_sample_id("s", "s", "sample ID", true, "", "str", cmd);
            TCLAP::ValueArg<uint32_t> arg_mapq_cutoff("m", "m", "MAPQ cutoff for alignments [20]", false, 20, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_baseq_cutoff("q", "q", "base quality cutoff for bases [13]", false, 13, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_evidence_allele_count_cutoff("e", "e", "evidence count cutoff for candidate allele [2]", false, 2, "int", cmd);
            TCLAP::ValueArg<double> arg_fractional_evidence_allele_count_cutoff("f", "f", "fractional evidence cutoff for candidate allele [0.1]", false, 0.1, "float", cmd);
            TCLAP::ValueArg<std::string> arg_variant_type("v", "v", "variant types [snps,mnps,indels]", false, "snps,mnps,indels", "str", cmd);
            TCLAP::ValueArg<std::string> arg_input_bam_file("b", "b", "input BAM file", true, "", "string", cmd);

            cmd.parse(argc, argv);

            input_bam_file = arg_input_bam_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            sample_id = arg_sample_id.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
            mapq_cutoff = arg_mapq_cutoff.getValue();
            baseq_cutoff = arg_baseq_cutoff.getValue();
            variant_type = arg_variant_type.getValue();
            evidence_allele_count_cutoff = arg_evidence_allele_count_cutoff.getValue();
            fractional_evidence_allele_count_cutoff = arg_fractional_evidence_allele_count_cutoff.getValue();
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
        exclude_flag = 0x0704;

        odr = new BAMOrderedReader(input_bam_file, intervals);
        s = bam_init1();

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bam_hdr_transfer_contigs_to_bcf_hdr(odr->hdr, odw->hdr);
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=E,Number=1,Type=Integer,Description=\"Number of reads containing evidence of the alternate allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Total number of reads at a candidate locus with reads that contain evidence of the alternate allele\">");
        bcf_hdr_add_sample(odw->hdr, sample_id.c_str());
        bcf_hdr_add_sample(odw->hdr, NULL);
        //bcf_hdr_set_n_sample(odw->hdr, 1);
                                
        v = NULL;

        std::vector<std::string> variant_types;
        split(variant_types, ",", variant_type);
        uint32_t vtype = 0;
        for (uint32_t i = 0; i<variant_types.size(); ++i)
        {
            if (variant_types[i] == "snps")
            {
                vtype |= SNP;
            }
            else if (variant_types[i] == "mnps")
            {
                vtype |= MNP;
            }
            else if (variant_types[i] == "indels")
            {
                vtype |= INDEL;
            }
            else if (variant_types[i] == "all")
            {
                vtype = SNP|MNP|INDEL;
            }
        }

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_reads = 0;
        no_overlapping_reads = 0;
        no_passed_reads = 0;
        no_exclude_flag_reads = 0;
        no_low_mapq_reads = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        faidx_t *fai = fai_load(ref_fasta_file.c_str());
        if (fai==NULL) 
        {
            fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
            exit(1);
        }
        variantHunter = new VariantHunter(vtype,
                                    evidence_allele_count_cutoff,
                                    fractional_evidence_allele_count_cutoff,
                                    baseq_cutoff,
                                    fai,
                                    odw);
    }

    void discover()
    {
        odw->write_hdr();

        //for tracking overlapping reads
        khash_t(rdict) *reads = kh_init(rdict);
        khiter_t k;
        int32_t ret;

        while (odr->read(s))
        {   
            ++no_reads;

            //this read is the first of the pair
            if (bam_get_mpos1(s) && (bam_get_tid(s)==bam_get_mtid(s)))
            {
                //first mate
                if (bam_get_mpos1(s)>bam_get_pos1(s))
                {
                    //overlapping
                    if (bam_get_mpos1(s)<=(bam_get_pos1(s) + bam_get_l_qseq(s) - 1))
                    {
                        //add read that has overlapping
                        //duplicate the record and perform the stitching later
                        char* qname = strdup(bam_get_qname(s));
                        k = kh_put(rdict, reads, qname, &ret);
                        if (!ret)
                        {
                            //already present
                            free(qname);
                        }
                        kh_val(reads, k) = {bam_get_pos1(s), bam_get_pos1(s)+bam_get_l_qseq(s)-1};
                    }
                }
                else
                {
                    //check overlap
                    //todo: perform stitching in future
                    if((k = kh_get(rdict, reads, bam_get_qname(s)))!=kh_end(reads))
                    {
                        if (kh_exist(reads, k))
                        {
                            free((char*)kh_key(reads, k));
                            kh_del(rdict, reads, k);
                            ++no_overlapping_reads;
                        }
                        //continue;
                    }
                }
            }

            if(bam_get_flag(s) & exclude_flag)
            {
                //1. unmapped
                //2. secondary alignment
                //3. not passing QC
                //4. PCR or optical duplicate
                ++no_exclude_flag_reads;
                continue;
            }

            if (bam_get_mapq(s) < mapq_cutoff)
            {
                //filter short aligments and those with too many indels (?)
                ++no_low_mapq_reads;
                continue;
            }

            if (0)
            {
               bam_print(s);
            }

            variantHunter->process_read(odr->hdr, s);



            variantHunter->output_read(odr->hdr, s);

            ++no_passed_reads;
        }

        odw->close();

    };

    void bam_print(bam1_t *s)
    {
        const char* chrom = bam_get_chrom(odr->hdr, s);
        uint32_t pos1 = bam_get_pos1(s);
        kstring_t seq = {0,0,0};
        bam_get_seq_string(s, &seq);
        uint32_t len = bam_get_l_qseq(s);
        kstring_t qual = {0,0,0};
        bam_get_qual_string(s, &qual);
        kstring_t cigar_string = {0,0,0};
        bam_get_cigar_string(s, &cigar_string);
        kstring_t cigar_expanded_string = {0,0,0};
        bam_get_cigar_expanded_string(s, &cigar_expanded_string);
        uint16_t flag = bam_get_flag(s);
        uint32_t mapq = bam_get_mapq(s);

        std::cerr << "##################" << "\n";
        std::cerr << "read no  : " << no_reads << "\n";
        std::cerr << "chrom-pos: " << chrom << "-" << pos1 << "\n";
        std::cerr << "read     : " << seq.s << "\n";
        std::cerr << "qual     : " << qual.s << "\n";
        std::cerr << "cigar_str: " << cigar_string.s << "\n";
        std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
        std::cerr << "len      : " << len << "\n";
        std::cerr << "mapq     : " << mapq << "\n";
        std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
        std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";

        if (seq.m) free(seq.s);
        if (qual.m) free(qual.s);
        if (cigar_string.m) free(cigar_string.s);
        if (cigar_expanded_string.m) free(cigar_expanded_string.s);
    }

    void print_options()
    {
        std::clog << "discover v" << version << "\n\n";

        std::clog << "options: [b] input BAM File               " << input_bam_file << "\n";
        std::clog << "         [o] output VCF File              " << output_vcf_file << "\n";
        std::clog << "         [s] sample ID                    " << sample_id << "\n";
        std::clog << "         [r] reference FASTA File         " << ref_fasta_file << "\n";
        std::clog << "         [m] MAPQ cutoff                  " << mapq_cutoff << "\n";
        std::clog << "         [q] base quality cutoff          " << baseq_cutoff << "\n";
        std::clog << "         [v] variant type(s)              " << variant_type << "\n";
        std::clog << "         [e] evidence cutoff              " << evidence_allele_count_cutoff << "\n";
        std::clog << "         [f] fractional evidence cutoff   " << fractional_evidence_allele_count_cutoff<< "\n";
        print_int_op("         [i] intervals                    ", intervals);
        std::clog << "\n";

    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. reads              : " << no_reads << "\n";
        std::clog << "       no. overlapping reads  : " << no_overlapping_reads << "\n";
        std::clog << "       no. low mapq reads     : " << no_low_mapq_reads << "\n";
        std::clog << "       no. passed reads       : " << no_passed_reads << "\n";
        std::clog << "       no. exclude flag reads : " << no_exclude_flag_reads << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void discover2(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.discover();
    igor.print_stats();
};
