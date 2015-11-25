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

#include "discover.h"

namespace
{

/**
 * For detecting overlapping reads.
 */
typedef struct
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

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
    bool ignore_md;
    int32_t debug;

    //options for selecting reads
    khash_t(rdict) *reads;

    //sample properties
    uint32_t ploidy;

    //read filters
    uint32_t read_mapq_cutoff;
    uint16_t read_exclude_flag;
    bool ignore_overlapping_read;

    //variables for keeping track of chromosome
    std::string chrom; //current chromosome
    int32_t tid;       // current sequence id in bam
    int32_t rid;       // current sequence id in bcf

    ///////
    //i/o//
    ///////
    BAMOrderedReader *odr;
    bam1_t *s;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    /////////////////////////
    //useful shared variables
    /////////////////////////
    int32_t *gt;

    /////////
    //stats//
    /////////
    uint32_t no_reads;
    uint32_t no_overlapping_reads;
    uint32_t no_passed_reads;
    uint32_t no_exclude_flag_reads;
    uint32_t no_low_mapq_reads;
    uint32_t no_unaligned_cigars;
    uint32_t no_malformed_del_cigars;
    uint32_t no_malformed_ins_cigars;
    uint32_t no_salvageable_ins_cigars;

    uint32_t no_snps;
    uint32_t no_ts;
    uint32_t no_tv;
    uint32_t no_insertions;
    uint32_t no_deletions;
    uint32_t no_left_soft_clips;
    uint32_t no_right_soft_clips;

    /////////
    //tools//
    /////////
    Pileup pileup;
    VariantFilter vf;
    LogTool lt;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Discovers variants from reads in a SAM/BAM/CRAM file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<uint32_t> arg_debug("d", "d", "debug [0]", false, 0, "int", cmd);
            TCLAP::SwitchArg arg_ignore_md("z", "z", "ignore MD tags [0]", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_sample_id("s", "s", "sample ID []", true, "", "str", cmd);

            TCLAP::ValueArg<uint32_t> arg_ploidy("p", "p", "ploidy [2]", false, 2, "int", cmd);

            //Reads
            TCLAP::ValueArg<uint32_t> arg_read_mapq_cutoff("t", "t", "MAPQ cutoff for alignments (>) [20]", false, 20, "int", cmd);
            TCLAP::SwitchArg arg_ignore_overlapping_read("l", "l", "ignore overlapping reads [false]", cmd, false);
            TCLAP::ValueArg<uint32_t> arg_read_exclude_flag("a", "a", "read exclude flag [0x0704]", false, 0x0704, "int", cmd);

            //Reference Bias
            TCLAP::ValueArg<float> arg_reference_bias("B", "B", "reference bias [0.1]", false, 0.1, "float", cmd);
            TCLAP::ValueArg<float> arg_lr_cutoff("C", "C", "likelihood ratio cutoff [0]", false, -1, "float", cmd);

            //SNP
            TCLAP::ValueArg<uint32_t> arg_snp_baseq_cutoff("q", "q", "base quality cutoff for bases [13]", false, 13, "int", cmd);
            TCLAP::ValueArg<uint32_t> arg_snp_e_cutoff("e", "e", "snp evidence count cutoff [2]", false, 2, "int", cmd);
            TCLAP::ValueArg<float> arg_snp_f_cutoff("f", "f", "snp fractional evidence cutoff [0]", false, 0, "float", cmd);
            TCLAP::ValueArg<float> arg_snp_desired_type_I_error("j", "j", "snp desired type I error [0.0]", false, 0, "float", cmd);
            TCLAP::ValueArg<float> arg_snp_desired_type_II_error("k", "k", "snp desired type II error [0.00001]", false, 0.00001, "float", cmd);

            //DEL
            TCLAP::ValueArg<uint32_t> arg_deletion_e_cutoff("u", "u", "deletion count cutoff [2]", false, 2, "int", cmd);
            TCLAP::ValueArg<float> arg_deletion_f_cutoff("v", "v", "deletion fractional evidence cutoff [0]", false, 0, "float", cmd);
            TCLAP::ValueArg<float> arg_deletion_desired_type_I_error("m", "m", "deletion desired type I error [0.0]", false, 0, "float", cmd);
            TCLAP::ValueArg<float> arg_deletion_desired_type_II_error("n", "n", "deletion desired type II error [0.00001]", false, 0.00001, "float", cmd);

            //INS
            TCLAP::ValueArg<uint32_t> arg_insertion_e_cutoff("g", "g", "insertion count cutoff [2]", false, 2, "int", cmd);
            TCLAP::ValueArg<float> arg_insertion_f_cutoff("h", "h", "insertion fractional evidence cutoff [0]", false, 0, "float", cmd);
            TCLAP::ValueArg<float> arg_insertion_desired_type_I_error("c", "c", "insertion desired type I error [0.0]", false, 0, "float", cmd);
            TCLAP::ValueArg<float> arg_insertion_desired_type_II_error("w", "w", "insertion desired type II error [0.00001]", false, 0.00001, "float", cmd);

            //SOFT CLIPS
            TCLAP::ValueArg<float> arg_sclip_mq_cutoff("x", "x", "soft clipped mean quality cutoff [0]", false, 0, "float", cmd);
            TCLAP::ValueArg<uint32_t> arg_sclip_u_cutoff("y", "y", "soft clipped unique sequences cutoff [0]", false, 1, "float", cmd);

            TCLAP::ValueArg<std::string> arg_input_bam_file("b", "b", "input SAM/BAM/CRAM file []", true, "", "string", cmd);

            cmd.parse(argc, argv);

            debug = arg_debug.getValue();
            ignore_md = arg_ignore_md.getValue();
            input_bam_file = arg_input_bam_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            sample_id = arg_sample_id.getValue();
            ploidy = arg_ploidy.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
            read_mapq_cutoff = arg_read_mapq_cutoff.getValue();
            ignore_overlapping_read = arg_ignore_overlapping_read.getValue();
            read_exclude_flag = arg_read_exclude_flag.getValue();

            vf.set_reference_bias(arg_reference_bias.getValue());
            vf.set_lr_cutoff(arg_lr_cutoff.getValue());

            vf.set_snp_baseq_cutoff(arg_snp_baseq_cutoff.getValue());
            vf.set_snp_e_cutoff(arg_snp_e_cutoff.getValue());
            vf.set_snp_f_cutoff(arg_snp_f_cutoff.getValue());
            vf.set_snp_desired_type_I_error(arg_snp_desired_type_I_error.getValue());
            vf.set_snp_desired_type_II_error(arg_snp_desired_type_II_error.getValue());

            vf.set_deletion_e_cutoff(arg_deletion_e_cutoff.getValue());
            vf.set_deletion_f_cutoff(arg_deletion_f_cutoff.getValue());
            vf.set_deletion_desired_type_I_error(arg_deletion_desired_type_I_error.getValue());
            vf.set_deletion_desired_type_II_error(arg_deletion_desired_type_II_error.getValue());

            vf.set_insertion_e_cutoff(arg_insertion_e_cutoff.getValue());
            vf.set_insertion_f_cutoff(arg_insertion_f_cutoff.getValue());
            vf.set_insertion_desired_type_I_error(arg_insertion_desired_type_I_error.getValue());
            vf.set_insertion_desired_type_II_error(arg_insertion_desired_type_II_error.getValue());

            vf.set_sclip_mq_cutoff(arg_sclip_mq_cutoff.getValue());
            vf.set_sclip_u_cutoff(arg_sclip_u_cutoff.getValue());

            vf.sync();
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

        //fails the following types
        //1. unmapped reads
        //2. secondary reads alignments
        //3. failed QC filter
        //4. duplicate
        //read_exclude_flag = 0x0704;

        odr = new BAMOrderedReader(input_bam_file, intervals, ref_fasta_file);
        s = bam_init1();

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bam_hdr_transfer_contigs_to_bcf_hdr(odr->hdr, odw->hdr);
        bcf_hdr_append(odw->hdr, "##QUAL=Variant score of the alternative allele likelihood ratio: -10 * log10 [P(Non variant)/P(Variant)].");
        bcf_hdr_append(odw->hdr, "##ALT=<ID=RSC,Description=\"Right Soft Clip\">");
        bcf_hdr_append(odw->hdr, "##ALT=<ID=LSC,Description=\"Left Soft Clip\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Soft clipped Sequence\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=E,Number=1,Type=Integer,Description=\"Number of reads containing evidence of the alternate allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Total number of reads at a candidate locus with reads that contain evidence of the alternate allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=MQS,Number=.,Type=Float,Description=\"Mean qualities of soft clipped bases.\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=STR,Number=.,Type=String,Description=\"Strands of soft clipped sequences.\">");

        bcf_hdr_add_sample(odw->hdr, sample_id.c_str());
        bcf_hdr_add_sample(odw->hdr, NULL);
        v = bcf_init();

        //for tracking overlapping reads
        reads = kh_init(rdict);

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_reads = 0;
        no_overlapping_reads = 0;
        no_passed_reads = 0;
        no_exclude_flag_reads = 0;
        no_low_mapq_reads = 0;
        no_unaligned_cigars = 0;
        no_malformed_del_cigars = 0;
        no_malformed_ins_cigars = 0;
        no_salvageable_ins_cigars = 0;

        no_snps = 0;
        no_ts = 0;
        no_tv = 0;
        no_insertions = 0;
        no_deletions = 0;
        no_left_soft_clips = 0;
        no_right_soft_clips = 0;

        //////////////////////////////////////
        //discovery variables initialization//
        //////////////////////////////////////
        chrom = "";
        tid = -1;
        rid = -1;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        pileup.set_reference(ref_fasta_file);
        pileup.set_debug(debug);
    }

    /**
     * Print BAM for debugging purposes.
     */
    void bam_print_key_values(bam_hdr_t *h, bam1_t *s)
    {
        const char* chrom = bam_get_chrom(h, s);
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

        uint8_t *aux;
        char* md = NULL;
        (aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(aux));

        std::cerr << "##################" << "\n";
        std::cerr << "chrom:pos: " << chrom << ":" << pos1 << "\n";
        std::cerr << "read     : " << seq.s << "\n";
        std::cerr << "qual     : " << qual.s << "\n";
        std::cerr << "cigar_str: " << cigar_string.s << "\n";
        std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
        std::cerr << "len      : " << len << "\n";
        std::cerr << "mapq     : " << mapq << "\n";
        std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
        std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";
        std::cerr << "md       : " << (aux?md:"") << "\n";
        std::cerr << "##################" << "\n";

        if (seq.m) free(seq.s);
        if (qual.m) free(qual.s);
        if (cigar_string.m) free(cigar_string.s);
        if (cigar_expanded_string.m) free(cigar_expanded_string.s);
    }

    /**
     * Filter reads.
     *
     * Returns true if read is failed.
     */
    bool filter_read(bam1_t *s)
    {
        khiter_t k;
        int32_t ret;

        if (ignore_overlapping_read)
        {
            //this read is part of a mate pair on the same contig
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
                    if((k = kh_get(rdict, reads, bam_get_qname(s)))!=kh_end(reads))
                    {
                        if (kh_exist(reads, k))
                        {
                            free((char*)kh_key(reads, k));
                            kh_del(rdict, reads, k);
                            ++no_overlapping_reads;
                        }
                        //set this on to remove overlapping reads.
                        return false;
                    }
                }
            }
        }
        //should we join them up?

        if(bam_get_flag(s) & read_exclude_flag)
        {
            //1. unmapped
            //2. secondary alignment
            //3. not passing QC
            //4. PCR or optical duplicate
            ++no_exclude_flag_reads;
            return false;
        }

        if (bam_get_mapq(s) < read_mapq_cutoff)
        {
            //filter short aligments and those with too many indels (?)
            ++no_low_mapq_reads;
            return false;
        }

        //*****************************************************************
        //should we have an assertion on the correctness of the bam record?
        //Is, Ds not sandwiched in M
        //leading and trailing Is - convert to S
        //no Ms!!!!!
        //*****************************************************************
        int32_t n_cigar_op = bam_get_n_cigar_op(s);
        if (n_cigar_op)
        {
            uint32_t *cigar = bam_get_cigar(s);
            bool seenM = false;
            int32_t last_opchr = '^';

            for (int32_t i = 0; i < n_cigar_op; ++i)
            {
                int32_t opchr = bam_cigar_opchr(cigar[i]);
                int32_t oplen = bam_cigar_oplen(cigar[i]);
                if (opchr=='S')
                {
                    if (i!=0 && i!=n_cigar_op-1)
                    {
                        std::cerr << "S issue\n";
                        bam_print_key_values(odr->hdr, s);
                        //++malformed_cigar;
                    }
                }
                else if (opchr=='M')
                {
                    seenM = true;
                }
                else if (opchr=='D')
                {
                    if (last_opchr!='M' || (i<=n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M'))
                    {
                        std::cerr << "D issue\n";
                        ++no_malformed_del_cigars;
                        bam_print_key_values(odr->hdr, s);
                    }
                }
                else if (opchr=='I')
                {
                    if (last_opchr!='M' || (i<n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M'))
                    {
                        if (last_opchr!='M')
                        {
                            if (last_opchr!='^' && last_opchr!='S')
                            {
                                std::cerr << "leading I issue\n";
                                bam_print_key_values(odr->hdr, s);
                                ++no_malformed_ins_cigars;
                            }
                            else
                            {
                                ++no_salvageable_ins_cigars;
                            }
                        }
                        else if (i==n_cigar_op-1)
                        {

                            ++no_salvageable_ins_cigars;
                        }
                        else if (i==n_cigar_op-2 && (bam_cigar_opchr(cigar[i+1])=='S'))
                        {
                            ++no_salvageable_ins_cigars;
                        }
                        else
                        {
                            std::cerr << "trailing I issue\n";
                            bam_print_key_values(odr->hdr, s);
                            ++no_malformed_ins_cigars;
                        }

                    }
                }

                last_opchr = opchr;
            }

            if (!seenM)
            {
                std::cerr << "NO! M issue\n";
                bam_print_key_values(odr->hdr, s);
                ++no_unaligned_cigars;
            }
        }

        //check to see that hash should be cleared when encountering new contig.
        //some bams may not be properly formed and contain orphaned sequences that
        //can retained in the hash
        if (bam_get_tid(s)!=tid)
        {
            for (k = kh_begin(reads); k != kh_end(reads); ++k)
            {
                if (kh_exist(reads, k))
                {
                    free((char*)kh_key(reads, k));
                    kh_del(rdict, reads, k);
                }
            }

            //tid is not updated here, it is handled by flush()
            //kh_destroy(rdict, h);
        }

        return true;
    }

    /**
     * Compute -10 * log10 likelihood ratio P(Non Variant)/P(Variant) for an alternative SNP allele.
     */
    float compute_snp_variant_score(std::vector<uint32_t>& REF_Q, std::vector<uint32_t>& ALT_Q)
    {
//        double pRR = 1;
//        double pRA = 1;
//        double pAA = 1;
//        double p;
//        double theta = 0.001;
//
//        for (uint32_t i=0; i<REF_Q.size(); ++i)
//        {
//            p = lt.pl2prob(REF_Q[i]);
//            pRR *= 1-p;
//            pRA *= 0.5;
//            pAA *= p;
//        }
//
//        for (uint32_t i=0; i<ALT_Q.size(); ++i)
//        {
//            p = lt.pl2prob(ALT_Q[i]);
//            pRR *= p;
//            pRA *= 0.5;
//            pAA *= 1-p;
//        }
//
//        double ln_lr = log10(pRR/((1-theta)*pRR+0.33*theta*pRA+0.67*theta*pAA));
//        ln_lr = ln_lr>0 ? 0 : ln_lr;
//
//        return (float) (-10 * ln_lr);
        
        ///////
        float lg_theta = -3; // theta = 0.001;
        float lg_one_minus_theta = -0.0004345118; // 1-theta = 0.999;
        float lg_0_5 = -0.30103;
        float lg_one_third = -0.4771213;
        float lg_two_thirds = -0.1760913;

        float lg_pRR = 0;
        float lg_pRA = 0;
        float lg_pAA = 0;

        for (uint32_t i=0; i<REF_Q.size(); ++i)
        {
            lg_pRR += lt.pl2pl_one_minus_p(REF_Q[i])/-10.0;
            lg_pRA += lg_0_5;
            lg_pAA += REF_Q[i]/-10.0;
        }
        
        for (uint32_t i=0; i<ALT_Q.size(); ++i)
        {
            lg_pRR += ALT_Q[i]/-10.0;
            lg_pRA += lg_0_5;
            lg_pAA += lt.pl2pl_one_minus_p(ALT_Q[i])/-10.0;
        }
        
        float lg_lr = lg_one_minus_theta + lg_pRR;
        lg_lr = lt.log10sum(lg_lr, lg_one_third+lg_theta+lg_pRA);
        lg_lr = lt.log10sum(lg_lr, lg_two_thirds+lg_theta+lg_pAA);
        lg_lr = lg_pRR - lg_lr;

        if (lg_lr>0)
        {
            return 0;
        }    
        else
        {
            return -10 * lg_lr;
        }
    }

    /**
     * Compute -10 * log10 likelihood ratio P(Non Variant)/P(Variant) for an alternative Indel allele.
     */
    float compute_indel_variant_score(uint32_t e, uint32_t n)
    {
        float ln_theta = -6.907755; // theta = 0.001;
        float ln_one_minus_theta = -0.0010005; // 1-theta = 0.999;
        float ln_one_third = -1.098612;
        float ln_two_thirds = -0.4054651;
        float ln_0_001 = -6.907755;   //equivalent to QUAL=30
        float ln_0_999 = -0.0010005; 
        float ln_0_5 = -0.6931472;

        float ln_pRR = (n-e) * ln_0_999 + e * ln_0_001;
        float ln_pRA = n * ln_0_5;
        float ln_pAA = e * ln_0_999 + (n-e) * ln_0_001;

        float ln_lr = ln_one_minus_theta + ln_pRR;
        ln_lr = logspace_add(ln_lr, ln_one_third+ln_theta+ln_pRA); //call from Rmath logspace arithmetic.
        ln_lr = logspace_add(ln_lr, ln_two_thirds+ln_theta+ln_pAA);
        ln_lr = ln_pRR - ln_lr;

        if (ln_lr>0)
        {
            return 0;
        }    
        else
        {
            return -10 * ln_lr * 0.4342945;
        }
    }

    /**
     * Write out pileupPosition as a VCF entry if it contains a variant.
     */
    void write_to_vcf(uint32_t rid, uint32_t gpos1, PileupPosition& p)
    {
        int32_t gts[2] = {0x0002,0x0004};

        if (p.R=='N' || p.R=='X')
        {
            return;
        }

        if (p.X[1]+p.X[2]+p.X[4]+p.X[8]+p.X[15]>p.E+p.N ||
           p.D.size()+p.I.size()>p.N)
        {
            std::cerr << "*******************\n";
            std::cerr << "Evidence count ISSUES\n";
            std::cerr << "#snp alts : " << (p.X[1]+p.X[2]+p.X[4]+p.X[8]+p.X[15]) << "\n";
            std::cerr << "#del alts : " << p.D.size() << "\n";
            std::cerr << "#ins alts : " << p.I.size() << "\n";
            std::cerr << "#reads (N): " << p.N << "\n";
            std::cerr << "#reads (E): " << p.E << "\n";
            p.print(gpos1);
            std::cerr << "*******************\n";
        }

        std::string alleles = "";
        int32_t N = 0;
        int32_t E = 0;
        float variant_score = 0;
        kstring_t new_alleles = {0,0,0};

        p.E = 0;
        p.F = 0;

        if (p.X[1])
        {
            if (vf.filter_snp(p.X[1], p.N+p.E-p.F))
            {
                bcf_clear(v);
                bcf_set_rid(v, rid);
                bcf_set_pos1(v, gpos1);
                bcf_set_n_sample(v, 1);
                alleles.clear();
                alleles.append(1, p.R);
                alleles.append(1, ',');
                alleles.append(1, 'A');
                bcf_update_alleles_str(odw->hdr, v, alleles.c_str());
                E = p.X[1];
                bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                N = p.N+p.E;
                bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                variant_score = compute_snp_variant_score(p.REF_Q, p.ALT_Q);
                bcf_set_qual(v, variant_score);
                odw->write(v);

                ++no_snps;

                if (p.R=='G')
                {
                    ++no_ts;
                }
                else
                {
                    ++no_tv;
                }
            }
        }

        if (p.X[2])
        {
            if (vf.filter_snp(p.X[2], p.N+p.E-p.F))
            {
                bcf_clear(v);
                bcf_set_rid(v, rid);
                bcf_set_pos1(v, gpos1);
                bcf_set_n_sample(v, 1);
                v->unpacked = 1;
                alleles.clear();
                alleles.append(1, p.R);
                alleles.append(1, ',');
                alleles.append(1, 'C');
                bcf_update_alleles_str(odw->hdr, v, alleles.c_str());
                E = p.X[2];
                bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                N = p.N+p.E;
                bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                variant_score = compute_snp_variant_score(p.REF_Q, p.ALT_Q);
                bcf_set_qual(v, variant_score);
                odw->write(v);

                ++no_snps;

                if (p.R=='T')
                {
                    ++no_ts;
                }
                else
                {
                    ++no_tv;
                }
            }
        }

        if (p.X[4])
        {
            if (vf.filter_snp(p.X[4], p.N+p.E-p.F))
//            if (vf.filter_snp(p.X[4], p.N+p.E))
            {
                bcf_clear(v);
                bcf_set_rid(v, rid);
                bcf_set_pos1(v, gpos1);
                bcf_set_n_sample(v, 1);
                alleles.clear();
                alleles.append(1, p.R);
                alleles.append(1, ',');
                alleles.append(1, 'G');
                bcf_update_alleles_str(odw->hdr, v, alleles.c_str());
                E = p.X[4];
                bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                N = p.N+p.E;
                bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                variant_score = compute_snp_variant_score(p.REF_Q, p.ALT_Q);
                bcf_set_qual(v, variant_score);
                odw->write(v);

                ++no_snps;

                if (p.R=='A')
                {
                    ++no_ts;
                }
                else
                {
                    ++no_tv;
                }
            }
        }

        if (p.X[8])
        {
            if (vf.filter_snp(p.X[8], p.N+p.E-p.F))
//            if (vf.filter_snp(p.X[8], p.N+p.E))
            {
                bcf_clear(v);
                bcf_set_rid(v, rid);
                bcf_set_pos1(v, gpos1);
                bcf_set_n_sample(v, 1);
                alleles.clear();
                alleles.append(1, p.R);
                alleles.append(1, ',');
                alleles.append(1, 'T');
                bcf_update_alleles_str(odw->hdr, v, alleles.c_str());
                E = p.X[8];
                bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                N = p.N+p.E;
                bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                variant_score = compute_snp_variant_score(p.REF_Q, p.ALT_Q);
                bcf_set_qual(v, variant_score);
                odw->write(v);

                ++no_snps;

                if (p.R=='C')
                {
                    ++no_ts;
                }
                else
                {
                    ++no_tv;
                }
            }
        }

        if (p.D.size()!=0)
        {
            for (std::map<std::string, uint32_t>::iterator i = p.D.begin(); i!=p.D.end(); ++i)
            {
                const std::string& del = i->first;
                if (del.find_first_of("Nn")==std::string::npos)
                {
                    E = i->second;
                    N = p.N;

                    if (vf.filter_del(E, N))
                    {
                        bcf_clear(v);
                        bcf_set_rid(v, rid);
                        bcf_set_pos1(v, gpos1);
                        bcf_set_n_sample(v, 1);
                        alleles.clear();
                        alleles.append(1, p.R);
                        alleles.append(del);
                        alleles.append(1, ',');
                        alleles.append(1, p.R);
                        bcf_update_alleles_str(odw->hdr, v, alleles.c_str());
                        bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                        bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                        variant_score = compute_indel_variant_score(E, N);
                        bcf_set_qual(v, variant_score);
                        odw->write(v);

                        ++no_deletions;
                    }
                }
            }
        }

        if (p.I.size()!=0)
        {
            for (std::map<std::string, uint32_t>::iterator i = p.I.begin(); i!=p.I.end(); ++i)
            {
                const std::string& ins = i->first;
                if (ins.find_first_of("Nn")==std::string::npos)
                {
                    E = i->second;
                    N = p.N;

                    if (vf.filter_ins(E, N))
                    {
                        bcf_clear(v);
                        bcf_set_rid(v, rid);
                        bcf_set_pos1(v, gpos1);
                        bcf_set_n_sample(v, 1);
                        alleles.clear();
                        alleles.append(1, p.R);
                        alleles.append(1, ',');
                        alleles.append(1, p.R);
                        alleles.append(ins);
                        bcf_update_alleles_str(odw->hdr, v, alleles.c_str());
                        bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                        bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                        variant_score = compute_indel_variant_score(E, N);
                        bcf_set_qual(v, variant_score);
                        odw->write(v);

                        ++no_insertions;
                    }
                }
            }
        }

        if (p.J.size()!=0)
        {
            if (false && p.J.size()>=vf.get_sclip_u_cutoff())
            {
                for (std::map<std::string, SoftClipInfo>::iterator i = p.J.begin(); i!=p.J.end(); ++i)
                {
                    bcf_clear(v);
                    bcf_set_rid(v, rid);
                    bcf_set_pos1(v, gpos1);
                    bcf_set_n_sample(v, 1);

                    new_alleles.l = 0;
                    kputc(p.R, &new_alleles);
                    kputc(',', &new_alleles);
                    kputs("<LSC:", &new_alleles);
                    kputw(i->first.size(), &new_alleles);
                    kputc('>', &new_alleles);
                    bcf_update_alleles_str(odw->hdr, v, new_alleles.s);

                    bcf_update_info_string(odw->hdr, v, "SEQ", i->first.c_str());

                    bcf_update_genotypes(odw->hdr, v, &gts, ploidy);
                    SoftClipInfo& info = i->second;
                    uint32_t no = info.no;
                    E = no;
                    bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                    N = p.N+p.E;
                    bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                    bcf_update_format_float(odw->hdr, v, "MQS", &info.mean_quals[0], no);
                    bcf_update_format_char(odw->hdr, v, "STR", &info.strands[0], no);
                    odw->write(v);

                    ++no_left_soft_clips;
                }
            }
        }

        if (p.K.size()!=0)
        {
            if (false && p.K.size()>=vf.get_sclip_u_cutoff())
            {
                for (std::map<std::string, SoftClipInfo>::iterator i = p.K.begin(); i!=p.K.end(); ++i)
                {
                    bcf_clear(v);
                    bcf_set_rid(v, rid);
                    bcf_set_pos1(v, gpos1);
                    bcf_set_n_sample(v, 1);

                    new_alleles.l = 0;
                    kputc(p.R, &new_alleles);
                    kputc(',', &new_alleles);
                    kputs("<RSC:", &new_alleles);
                    kputw(i->first.size(), &new_alleles);
                    kputc('>', &new_alleles);
                    bcf_update_alleles_str(odw->hdr, v, new_alleles.s);

                    bcf_update_info_string(odw->hdr, v, "SEQ", i->first.c_str());

                    bcf_update_genotypes(odw->hdr, v, &gts, ploidy);
                    SoftClipInfo& info = i->second;
                    uint32_t no = info.no;
                    E = no;
                    bcf_update_format_int32(odw->hdr, v, "E", &E, 1);
                    N = p.N+p.E;
                    bcf_update_format_int32(odw->hdr, v, "N", &N, 1);
                    float* mean_quals = info.mean_quals.data();
                    bcf_update_format_float(odw->hdr, v, "MQS", mean_quals, no);
                    char* s = info.strands.data();
                    bcf_update_format_char(odw->hdr, v, "STR", s, no);
                    odw->write(v);

                    ++no_right_soft_clips;
                }
            }
        }

        if (new_alleles.m) free(new_alleles.s);
    }

    /**
     * Check if flushable.
     *
     * returns
     *    0 - not flushable
     *    1 - flushable
     *   -1 - flushable, must update with new sequence
     */
    int32_t flushable(bam1_t* s)
    {
        //different sequence, flush everything.
        if (pileup.get_tid()!=bam_get_tid(s))
        {
            return -1;
        }

        //position is smaller than window size
        if (pileup.get_window_size()>=bam_get_pos1(s))
        {
            return 0;
        }

        //some leeway for flushing the records
        if (bam_get_pos1(s)-pileup.get_window_size()>pileup.get_gbeg1())
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }

    /**
     * Flush records out before gpos1.
     */
    void flush(bam1_t* s)
    {
        int32_t ret = 0;
        if ((ret=flushable(s)))
        {
            if (ret==1)
            {
                if (debug>=3) std::cerr << "FLUSHING " << pileup.get_gbeg1() << " to " << bam_get_pos1(s)-pileup.get_window_size() << "\n";

                uint32_t cpos1 = pileup.get_gbeg1();
                uint32_t gpos1 = bam_get_pos1(s)-pileup.get_window_size();
                uint32_t lend0 = pileup.get_gend1()<gpos1 ? pileup.end() : pileup.g2i(gpos1);

//                uint32_t j;
//                for (j=pileup.end(); j!=pileup.begin(); j=pileup.inc(j,1))
//                {
//                    if(!pileup[j].is_cleared())
//                    {
//                        std::cerr << "NOT VALID!!!!\n";
//                        pileup[j].print();
//                    }
//                }

                uint32_t i;
                for (i=pileup.begin(); i!=lend0; i=pileup.inc(i,1))
                {
                    write_to_vcf(rid, cpos1, pileup[i]);
                    pileup[i].clear();
                    ++cpos1;
                }

                pileup.set_gbeg1(cpos1);
                pileup.set_beg0(i);
            }
            else //ret == -1
            {
                if (debug>=3) std::cerr << "FLUSHING " << pileup.get_gbeg1() << " to " << pileup.get_gend1() << "\n";

//                uint32_t j;
//                pileup.print_state();
//                for (j=pileup.end(); j!=pileup.begin(); j=pileup.inc(j,1))
//                {
//                    if(!pileup[j].is_cleared())
//                    {
//                        std::cerr << "NOT VALID!!!!\n";
//                        pileup[j].print();
//                    }
//                }

                uint32_t cpos1 = pileup.get_gbeg1();
                uint32_t i;
                for (i=pileup.begin(); i!=pileup.end(); i=pileup.inc(i,1))
                {
                    write_to_vcf(rid, cpos1, pileup[i]);
                    pileup[i].clear();
                    ++cpos1;
                }

                tid = bam_get_tid(s);
                chrom.assign(bam_get_chrom(odr->hdr, s));
                rid = bcf_hdr_name2id(odw->hdr, chrom.c_str());
                pileup.set_tid(tid);
                pileup.set_chrom(chrom);
                pileup.set_gbeg1(0);
                pileup.set_beg0(i);
            }
        }
    }

    /**
     * Flush records till pileup is empty.
     */
    void flush()
    {
        uint32_t cpos1 = pileup.get_gbeg1();
        uint32_t i;

        for (i=pileup.begin(); i!=pileup.end(); i=pileup.inc(i,1))
        {
            write_to_vcf(rid, cpos1, pileup[i]);
            pileup[i].clear();
            ++cpos1;
        }

        tid = -1;
        rid = -1;
        pileup.set_tid(-1);
        pileup.set_gbeg1(0);
        pileup.set_beg0(0);
        pileup.set_end0(0);
    }

    /**
     * Processes cigar string and MD tag to minimize access of the reference file.
     * Aggregates observed read information in pileup data structure.
     * Outputs to vcf.
     */
    void process_read(bam1_t *s)
    {
        if (debug>=1) bam_print_key_values(odr->hdr, s);

        flush(s);

        uint32_t tid = bam_get_tid(s);
        uint32_t pos1 = bam_get_pos1(s);
        uint8_t* seq = bam_get_seq(s);
        uint8_t* qual = bam_get_qual(s);
        int32_t l_qseq = bam_get_l_qseq(s);
        uint32_t* cigar = bam_get_cigar(s);
        char strand = bam_is_rev(s) ? '-' : '+';

        uint8_t *md_aux;
        char* md;

        if (!ignore_md && (md_aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(md_aux)))
        {
            //iterate cigar
            uint32_t n_cigar_op = bam_get_n_cigar_op(s);

            char* mdp = md;
            uint32_t cpos1 = pos1; //current 1 based genome position
            uint32_t spos0 = 0;    //current position in read sequence

            //variables for I's embedded in Matches in the MD tag
            uint32_t md_mlen_left = 0;

            if (n_cigar_op)
            {
                uint32_t *cigar = bam_get_cigar(s);
                bool seenM = false;

                if (debug>=3) pileup.print_state();
                for (uint32_t i = 0; i < n_cigar_op; ++i)
                {
                    uint32_t oplen = bam_cigar_oplen(cigar[i]);
                    char opchar = bam_cigar_opchr(cigar[i]);

                    if (debug) std::cerr << "CIGAR: " << oplen << " " << opchar << "\n";

                    if (opchar=='S')
                    {
                        if (i==n_cigar_op-1 || (i==0 && n_cigar_op>=2 && bam_cigar_opchr(cigar[1])=='M'))
                        {
                            //add to S evidence
                            std::string ins = "";
                            float mean_qual = 0;
                            for (size_t j=0; j<oplen ; ++j)
                            {
                                ins += bam_base2char(bam_seqi(seq, spos0+j));
                                mean_qual += qual[spos0+j];
                            }
                            mean_qual /= oplen;

                            if (cpos1==pos1)
                            {
                                if (debug) std::cerr << "\t\t\tadding LSCLIP: " << cpos1 << "\t" << ins << " {" << mean_qual << "}\n";
                                pileup.add_lsclip(cpos1, ins, mean_qual, strand);
                            }
                            else if (seenM)
                            {
                                if (debug) std::cerr << "\t\t\tadding RSCLIP: " << (cpos1-1) << "\t" << ins << " {" << mean_qual << "}\n";
                                pileup.add_rsclip(cpos1-1, ins, mean_qual, strand);
                            }
                        }

                        spos0 += oplen;
                    }
                    else if (opchar=='M')
                    {
                        uint32_t lpos1 = cpos1; // we need this because M contains matches and mismatches
                        uint32_t sspos0 = spos0; // we need this because M contains matches and mismatches
                        uint32_t mlen = oplen;
                        seenM = true;

                        if (debug) std::cerr << "\t\tmd len left : " << md_mlen_left << "\n";
                        if (debug) std::cerr << "\t\tmlen : " << mlen << "\n";
                        if (debug) std::cerr << "\t\tmdp : " << mdp << "\n";

                        //left over MD matches to handle.
                        if (md_mlen_left)
                        {
                            uint32_t ilen = md_mlen_left<=mlen ? md_mlen_left : mlen;
                            pileup.add_ref(lpos1, sspos0, ilen, seq);

                            if (debug)
                            {
                                uint32_t gbeg1 = lpos1;
                                uint32_t gend1 = lpos1+ilen-1;

                                std::cerr << "\t\t\tadding REF: " << gbeg1 << "-" << gend1 << ":";
                                for (size_t j=sspos0; j<=(sspos0+ilen-1); ++j)
                                {
                                    std::cerr << (bam_base2char(bam_seqi(seq, j)));
                                }
                                std::cerr << " (" << gend1-gbeg1+1 << ") [" <<  mlen-ilen << "]\n";
                            }

                            lpos1 += ilen;
                            sspos0 += ilen;

                            if (md_mlen_left>=mlen)
                            //yet another insertion
                            {
                                md_mlen_left -= ilen;
                                cpos1 += ilen;
                                spos0 += ilen;
                                continue;
                            }
                            //a snp next
                            else
                            {
                                md_mlen_left = 0;
                                mlen -= ilen;
                                //go to loop in the next section
                            }
                        }

                        while (*mdp)
                        {
                            if (isalpha(*mdp)) //SNPs
                            {
                                char ref = toupper(*mdp);
                                char alt = (bam_base2char(bam_seqi(seq, spos0+(lpos1-cpos1))));
                                if (debug) std::cerr << "\tMD: Mismatch " << ref << "\n";
                                if (debug) std::cerr << "\t\t\tadding SNP: " << lpos1 << ":" << ref << "/" << alt << " [" << (mlen-1)<< "]\n";
//                                if (qual[sspos0]>vf.get_snp_baseq_cutoff())
//                                {
                                    pileup.add_snp(lpos1, ref, alt, qual[sspos0], vf.get_snp_baseq_cutoff());
//                                }
                                ++lpos1;
                                ++mdp;
                                ++sspos0;
                                --mlen;
                            }
                            else if (isdigit(*mdp)) //matches
                            {
                                char* end = 0;
                                int32_t len = std::strtol(mdp, &end, 10);
                                mdp = end;

                                if (debug) std::cerr << "\tMD: Match " << len << "\n";

                                if (len)
                                {
                                    uint32_t ilen = len<=mlen ? len : mlen;

                                    if (debug)
                                    {
                                        uint32_t gbeg1 = lpos1;
                                        uint32_t gend1 = lpos1+ilen-1;

                                        //std::cerr << "\t\t\tadding REF: " << gbeg1 << "-" << gend1 << ":";
                                        for (size_t j=sspos0; j<=(sspos0+ilen-1); ++j)
                                        {
                                            std::cerr << (bam_base2char(bam_seqi(seq, j)));
                                        }
                                        std::cerr << " (" << gend1-gbeg1+1 << ") [" <<  mlen-ilen << "]\n";
                                    }

                                    pileup.add_ref(lpos1, sspos0, ilen, seq);

                                    lpos1 += ilen;
                                    sspos0 += ilen;

                                    //next up an insertion
                                    if (len>mlen)
                                    {
                                        md_mlen_left = len - mlen;
                                        break;
                                    }
                                    else
                                    {
                                        mlen -= ilen;
                                    }
                                }
                            }
                            else // deletion
                            {
                                break;
                            }

                            if (mlen==0)
                            {
                                break;
                            }
                        }

                        //note that only insertions, matches and mismatches can only occur here.

                        cpos1 += oplen;
                        spos0 += oplen;
                    }
                    else if (opchar=='D')
                    {
                        bool is_del = false;

                        if (*mdp=='0') ++mdp;

                        if (*mdp!='^')
                        {
                            bam_print_key_values(odr->hdr, s);
                            std::cerr << "mdp: " << mdp << "\n";
                            std::cerr << "inconsistent MD and cigar, deletion does not occur at the right place.\n";
                            exit(1);
                        }
                        else
                        {
                            ++mdp;
                            std::string del = "";
                            while (isalpha(*mdp))
                            {
                                del += toupper(*mdp);
                                ++mdp;
                            }

                            if (debug) std::cerr << "\t\t\tadding DEL: " << (cpos1-1) << " " << del << "\n";
                            pileup.add_del((cpos1-1), del);

                            cpos1 += oplen;
                        }
                    }
                    else if (opchar=='I')
                    {
                        //leading Is
                        if (!seenM)
                        {
//                            //add to S evidence
//                            std::string ins = "";
//                            float mean_qual = 0;
//                            for (size_t j=0; j<oplen ; ++j)
//                            {
//                                ins += bam_base2char(bam_seqi(seq, spos0+j));
//                                mean_qual += qual[spos0+j];
//                            }
//                            mean_qual /= oplen;
//
//                            if (mean_qual>sclip_mq_cutoff)
//                            {
//                                if (debug) std::cerr << "\t\t\tadding LSCLIP: " << cpos1 << "\t" << ins << " {" << mean_qual << "}\n";
//                                pileup.add_lsclip(cpos1, ins, mean_qual, strand);
//                            }

                            spos0 += oplen;
                        }
                        //trailing Is
                        else if (i==n_cigar_op-1 || (i+2==n_cigar_op && bam_cigar_opchr(cigar[n_cigar_op-1])=='S'))
                        {
                            //bam_print_key_values(odr->hdr, s);
                            spos0 += oplen;
                        }
                        else
                        {
                            //insertions are not present in MD tags
                            //may be handled independently of future matches
                            std::string ins = "";
                            for (size_t i=0; i<oplen ; ++i)
                            {
                                ins += bam_base2char(bam_seqi(seq, spos0+i));
                            }

                            if (debug) std::cerr << "\t\t\tadding INS: " << (cpos1-1) << " " << ins  << "\n";
                            pileup.add_ins((cpos1-1), ins, pos1);

                            spos0 += oplen;
                        }
                    }
                    else
                    {
                        std::cerr << "never seen before state " << opchar << "\n";
                    }

                    if (debug>=2) pileup.print_state_extent();
                }

                //pileup.print_state();
                //update last matching base
                if (seenM)
                {
                    pileup.update_read_end(cpos1-1);
                }
            }
        }
        else //don't use MD tag
        {
            //iterate cigar
            uint32_t n_cigar_op = bam_get_n_cigar_op(s);
            uint32_t cpos1 = pos1; //current 1 based genome position
            uint32_t spos0 = 0;    //current position in read sequence
            bool seenM = false;

            if (n_cigar_op)
            {
                uint32_t *cigar = bam_get_cigar(s);

                if (debug>=3) pileup.print_state();

                for (uint32_t i = 0; i < n_cigar_op; ++i)
                {
                    uint32_t oplen = bam_cigar_oplen(cigar[i]);
                    char opchar = bam_cigar_opchr(cigar[i]);

                    if (debug) std::cerr << "CIGAR: " << oplen << " " << opchar << "\n";

                    if (opchar=='S')
                    {
                        if (i==n_cigar_op-1 || (i==0 && n_cigar_op>=2 && bam_cigar_opchr(cigar[1])=='M'))
                        {
                            std::string ins = "";
                            float mean_qual = 0;
                            for (size_t j=0; j<oplen ; ++j)
                            {
                                ins += bam_base2char(bam_seqi(seq, spos0+j));
                                mean_qual += qual[spos0+j];
                            }
                            mean_qual /= oplen;

                            if (cpos1==pos1)
                            {
                                if (debug) std::cerr << "\t\t\tadding LSC: " << cpos1 << "\t" << ins << " {" << mean_qual << "}\n";
                                pileup.add_LSC(cpos1, ins, mean_qual, strand);

                            }
                            else
                            {
                                if (debug) std::cerr << "\t\t\tadding RSC: " << (cpos1-1) << "\t" << ins << " {" << mean_qual << "}\n";
                                pileup.add_RSC(cpos1-1, ins, mean_qual, strand);
                            }
                        }

                        spos0 += oplen;
                    }
                    else if (opchar=='M')
                    {
                        seenM = true;
                        if (debug) std::cerr << "\t\t\tadding matches: " << cpos1 << "\t" << spos0 << " (" << oplen << ")\n";
                        pileup.add_M(cpos1, spos0, oplen, seq, qual, vf.get_snp_baseq_cutoff());

                        cpos1 += oplen;
                        spos0 += oplen;
                    }
                    else if (opchar=='D')
                    {
                        std::string del = "";
                        for (size_t i=0; i<oplen ; ++i)
                        {
                            del += bam_base2char(bam_seqi(seq, spos0+i));
                        }

                        if (debug) std::cerr << "\t\t\tadding DEL: " << (cpos1-1) << " " << del << "\n";
                        pileup.add_D((cpos1-1), oplen);

                        cpos1 += oplen;
                    }
                    else if (opchar=='I')
                    {
                        //leading Is
                        if (!seenM)
                        {
//                            //add to S evidence
//                            std::string ins = "";
//                            float mean_qual = 0;
//                            for (size_t j=0; j<oplen ; ++j)
//                            {
//                                ins += bam_base2char(bam_seqi(seq, spos0+j));
//                                mean_qual += qual[spos0+j];
//                            }
//                            mean_qual /= oplen;
//
//                            if (mean_qual>sclip_mq_cutoff)
//                            {
//                                if (debug) std::cerr << "\t\t\tadding LSCLIP: " << cpos1 << "\t" << ins << " {" << mean_qual << "}\n";
//                                pileup.add_LSC(cpos1, ins, mean_qual, strand);
//                            }

                            spos0 += oplen;
                        }
                        //trailing Is
                        else if (i==n_cigar_op-1 || (i+2==n_cigar_op && bam_cigar_opchr(cigar[n_cigar_op-1])=='S'))
                        {
                            //bam_print_key_values(odr->hdr, s);
                            spos0 += oplen;
                        }
                        else
                        {
                            std::string ins = "";
                            for (size_t i=0; i<oplen ; ++i)
                            {
                                ins += bam_base2char(bam_seqi(seq, spos0+i));
                            }

                            if (debug) std::cerr << "\t\t\tadding INS: " << (cpos1-1) << " " << ins  << "\n";
                            pileup.add_I((cpos1-1), ins, pos1);

                            spos0 += oplen;
                        }
                    }
                    else
                    {
                        std::cerr << "never seen before state " << opchar << "\n";
                    }

                    if (debug>=2)  pileup.print_state_extent();
                }

                //pileup.print_state();
                //update last matching base
                if (seenM)
                {
                    pileup.update_read_end(cpos1-1);
                }
            }
        }
    }

    /**
     * Discover variants.
     */
    void discover()
    {
        odw->write_hdr();
        while (odr->read(s))
        {
            ++no_reads;

            if (!filter_read(s))
            {
                continue;
            }

            process_read(s);
            if (debug>=3) pileup.print_state();
            ++no_passed_reads;

            //if (no_passed_reads==1) break;

            if ((no_reads & 0x0000FFFF) == 0)
            {
                std::cerr << pileup.get_chrom() << ":" << pileup.get_gbeg1() << "\n";
            }
        }
        flush();
        odw->close();
    };

    void print_options()
    {
        std::clog << "discover v" << version << "\n\n";

        std::clog << "options: [b] input BAM File                       " << input_bam_file << "\n";
        std::clog << "         [o] output VCF File                      " << output_vcf_file << "\n";
        std::clog << "         [s] sample ID                            " << sample_id << "\n";
        std::clog << "         [r] reference FASTA File                 " << ref_fasta_file << "\n";
        std::clog << "         [p] ploidy                               " << ploidy << "\n";
        std::clog << "         [z] ignore MD tags                       " << (ignore_md ? "true": "false") << "\n";
        print_int_op("         [i] intervals                            ", intervals);
        std::clog << "\n";
        std::clog << "         [B] reference bias                       " << vf.get_reference_bias() << "\n";
        std::clog << "         [C] likelihood ratio cutoff              " << vf.get_lr_cutoff() << "\n";
        std::clog << "\n";
        std::clog << "         [t] read mapping quality cutoff          " << read_mapq_cutoff << "\n";
        std::clog << "         [l] ignore overlapping read              " << (ignore_overlapping_read ? "true" : "false") << "\n";
        std::clog << "         [a] read flag filter                     " << std::showbase << std::hex << read_exclude_flag << std::dec << "\n";
        std::clog << "\n";
        std::clog << "         [q] snp base quality cutoff              " << vf.get_snp_baseq_cutoff() << "\n";
        std::clog << "         [e] snp evidence cutoff                  " << vf.get_snp_e_cutoff() << "\n";
        std::clog << "         [f] snp fractional evidence cutoff       " << vf.get_snp_f_cutoff() << "\n";
        std::clog << "         [j] snp desired_type_I_error             " << vf.get_snp_desired_type_I_error() << "\n";
        std::clog << "         [k] snp desired_type_II_error            " << vf.get_snp_desired_type_II_error() << "\n";
        std::clog << "\n";
        std::clog << "         [u] deletion evidence cutoff             " << vf.get_deletion_e_cutoff() << "\n";
        std::clog << "         [v] deletion fractional evidence cutoff  " << vf.get_deletion_f_cutoff() << "\n";
        std::clog << "         [m] deletion desired_type_I_error        " << vf.get_deletion_desired_type_I_error() << "\n";
        std::clog << "         [n] deletion desired_type_II_error       " << vf.get_deletion_desired_type_II_error() << "\n";
        std::clog << "\n";
        std::clog << "         [g] insertion evidence cutoff            " << vf.get_insertion_e_cutoff() << "\n";
        std::clog << "         [h] insertion fractional evidence cutoff " << vf.get_insertion_f_cutoff() << "\n";
        std::clog << "         [c] insertion desired_type_I_error       " << vf.get_insertion_desired_type_I_error() << "\n";
        std::clog << "         [w] insertion desired_type_II_error      " << vf.get_insertion_desired_type_II_error() << "\n";
        std::clog << "\n";
        std::clog << "         [x] soft clip mean quality cutoff        " << vf.get_sclip_mq_cutoff() << "\n";
        std::clog << "         [y] soft clip unique events cutoff       " << vf.get_sclip_u_cutoff() << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. reads                    : " << no_reads << "\n";
        std::clog << "       no. overlapping reads        : " << no_overlapping_reads << "\n";
        std::clog << "       no. low mapq reads           : " << no_low_mapq_reads << "\n";
        std::clog << "       no. passed reads             : " << no_passed_reads << "\n";
        std::clog << "       no. exclude flag reads       : " << no_exclude_flag_reads << "\n";
        std::clog << "\n";
        std::clog << "       no. unaligned cigars         : " << no_unaligned_cigars << "\n";
        std::clog << "       no. malformed del cigars     : " << no_malformed_del_cigars << "\n";
        std::clog << "       no. malformed ins cigars     : " << no_malformed_ins_cigars << "\n";
        std::clog << "       no. salvageable ins cigars   : " << no_salvageable_ins_cigars << "\n";
        std::clog << "\n";
        std::clog << "       no. variants                 : " << (no_snps+no_insertions+no_deletions+no_left_soft_clips+no_right_soft_clips) << "\n";
        std::clog << "           no. snps (ts/tv)         : " << no_snps << " (" << std::fixed << std::setprecision(2) << (float)no_ts/no_tv << ")\n";
        std::clog << "               no. transitions      : " << no_ts << "\n";
        std::clog << "               no. transversions    : " << no_tv << "\n";
        std::clog << "           no. indels (ins/del)     : " << (no_insertions + no_deletions) << std::fixed << std::setprecision(2) << " (" << (float)no_insertions/no_deletions << ")\n";
        std::clog << "               no. insertions       : " << no_insertions << "\n";
        std::clog << "               no. deletions        : " << no_deletions << "\n";
        std::clog << "       no. soft clips               : " << (no_left_soft_clips+no_right_soft_clips) << "\n";
        std::clog << "               no. left soft clips  : " << no_left_soft_clips << "\n";
        std::clog << "               no. right soft clips : " << no_right_soft_clips << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {
        //clear hash of reads
        for (khiter_t k = kh_begin(reads); k != kh_end(reads); ++k)
        {
            if (kh_exist(reads, k))
            {
                free((char*)kh_key(reads, k));
                kh_del(rdict, reads, k);
            }
        }
    };

    private:
};

}

void discover(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.discover();
    igor.print_stats();
};
