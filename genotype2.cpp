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

#include "genotype2.h"

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

    std::string version;

    ///////////
    //options//
    ///////////
    std::string sample_id;
    std::string input_vcf_file;
    std::string input_sam_file;
    std::string output_vcf_file;
    std::string ref_fasta_file;
    std::string mode;
    bool ignore_md;
    int32_t debug;

    //variables for keeping track of chromosome
    std::string chrom; //current chromosome
    int32_t tid; // current sequence id in bam
    int32_t rid; // current sequence id in bcf

    //read filters
    uint32_t read_mapq_cutoff;
    uint16_t read_exclude_flag;
    bool ignore_overlapping_read;

    ///////
    //i/o//
    ///////
    BAMOrderedReader *odr;
    BCFGenotypingBufferedReader *gbr;
    BCFOrderedWriter *odw;

    std::vector<GenomeInterval> intervals;

    //options for selecting reads
    khash_t(rdict) *reads;

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

    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;
    uint32_t no_vntrs_genotyped;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Genotypes SNPs, Indels, VNTRs for each sample.\n";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            //Reads
            TCLAP::ValueArg<uint32_t> arg_read_mapq_cutoff("t", "t", "MAPQ cutoff for alignments (>=) [0]", false, 0, "int", cmd);
            TCLAP::SwitchArg arg_ignore_overlapping_read("l", "l", "ignore overlapping reads [false]", cmd, false);
            TCLAP::ValueArg<uint32_t> arg_read_exclude_flag("a", "a", "read exclude flag [0x0704]", false, 0x0704, "int", cmd);


            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_input_sam_file("b", "b", "input SAM/BAM/CRAM file []", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file", false, "-", "string", cmd);
            TCLAP::ValueArg<std::string> arg_sample_id("s", "s", "sample ID []", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> arg_mode("m", "m", "mode [d]\n"
                 "              d : iterate by read for dense genotyping.\n"
                 "                 (e.g. 50m variants close to one another).\n"
                 "              s : iterate by sites for sparse genotyping.\n"
                 "                 (e.g. 100 variants scattered over the genome).",
                 false, "d", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference FASTA file []", true, "", "string", cmd);
            TCLAP::ValueArg<uint32_t> arg_debug("d", "d", "debug [0]", false, 0, "int", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            mode = arg_mode.getValue();
            input_vcf_file = arg_input_vcf_file.getValue();
            input_sam_file = arg_input_sam_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            sample_id = arg_sample_id.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            ref_fasta_file = arg_ref_fasta_file.getValue();
            debug = arg_debug.getValue();

            read_mapq_cutoff = arg_read_mapq_cutoff.getValue();
            ignore_overlapping_read = arg_ignore_overlapping_read.getValue();
            read_exclude_flag = arg_read_exclude_flag.getValue();

        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }

        //////////////////////
        //i/o initialization//
        //////////////////////

        //fails the following types
        //1. unmapped reads
        //2. secondary reads alignments
        //3. failed QC filter
        //4. duplicate
        //read_exclude_flag = 0x0704;

        //input sam
        odr = new BAMOrderedReader(input_sam_file, intervals);

        //input vcf
        gbr = new BCFGenotypingBufferedReader(input_vcf_file, intervals);

        //output vcf
        odw = new BCFOrderedWriter(output_vcf_file);
        bcf_hdr_transfer_contigs(gbr->odr->hdr, odw->hdr);
        bcf_hdr_add_sample(odw->hdr, strdup(sample_id.c_str()));
        bcf_hdr_add_sample(odw->hdr, NULL);

        //INFO fields
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "MOTIF", "1", "String", "Canonical motif in an VNTR or homopolymer", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU", "1", "String", "Repeat unit in a VNTR or homopolymer", true);

        bcf_hdr_append_info_with_backup_naming(odw->hdr, "RL", "1", "Float", "Reference repeat unit length", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "LL", "1", "Float", "Longest repeat unit length", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "CONCORDANCE", "1", "Float", "Concordance of repeat unit.", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU_COUNTS", "2", "Integer", "Number of exact repeat units and total number of repeat units.", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "FLANKS", "2", "Integer", "Left and right flank positions of the Indel, left/right alignment invariant, not necessarily equal to POS.", true);

        bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_RL", "1", "Float", "Fuzzy reference repeat unit length", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_LL", "1", "Float", "Fuzzy longest repeat unit length", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_CONCORDANCE", "1", "Float", "Fuzzy concordance of repeat unit.", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_RU_COUNTS", "2", "Integer", "Fuzzy number of exact repeat units and total number of repeat units.", true);
        bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_FLANKS", "2", "Integer", "Fuzzy left and right flank positions of the Indel, left/right alignment invariant, not necessarily equal to POS.", true);

        bcf_hdr_append_info_with_backup_naming(odw->hdr, "TR", "1", "String", "Tandem repeat associated with this indel.", true);

        bcf_hdr_append(odw->hdr, "##INFO=<ID=LARGE_REPEAT_REGION,Number=0,Type=Flag,Description=\"Very large repeat region, vt only detects up to 1000bp long regions.\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=FLANKSEQ,Number=1,Type=String,Description=\"Flanking sequence 10bp on either side of detected repeat region.\">");
        
        //FILTERS
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_snp,Description=\"Overlaps with snp\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_indel,Description=\"Overlaps with indel\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_vntr,Description=\"Overlaps with VNTR\">");

        //COMMON
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED scaled genotype likelihoods\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depth\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=ADF,Number=A,Type=Integer,Description=\"Allele Depth (Forward strand)\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=ADR,Number=A,Type=Integer,Description=\"Allele Depth (Reverse strand)\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");

        //NONREF SNP
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=BQ,Number=.,Type=Integer,Description=\"Base Qualities\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=MQ,Number=.,Type=Integer,Description=\"Phred-scaled Map Qualities\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=CY,Number=.,Type=Integer,Description=\"Cycle of base\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=ST,Number=1,Type=String,Description=\"Strand of allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AL,Number=.,Type=Integer,Description=\"Alleles - 0,1,2,... for reference and alternate alleles. "
                                 "-1 : other allele, -2 : deletion, -3 : insertion\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=NM,Number=.,Type=Integer,Description=\"Number of mismatches per read\">");

        //NONREF INDEL
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=RQ,Number=.,Type=Integer,Description=\"Phred-scaled Reference Allele Qualities for Indels\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AQ,Number=.,Type=Integer,Description=\"Phred-scaled Alternative Allele Qualities for Indels, the number of entries is ploidy*no_alleles\">");

        //VNTR
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=CG,Number=.,Type=Float,Description=\"Repeat count genotype\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=CT,Number=.,Type=Float,Description=\"Repeat counts\">");

        //REF
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=BQSUM,Number=1,Type=Integer,Description=\"Sum of Base Qualities\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AQSUM,Number=A,Type=Integer,Description=\"Sum of Allele Likelihoods\">");

        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Depth of forward reference alleles\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DPR,Number=1,Type=Integer,Description=\"Depth of reverse reference alleles\">");

        odw->write_hdr();


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

        no_snps_genotyped = 0;
        no_indels_genotyped = 0;
        no_vntrs_genotyped = 0;

        //for tracking overlapping reads
        reads = kh_init(rdict);

        //////////////////////////////////////
        //discovery variables initialization//
        //////////////////////////////////////
        chrom = "";
        tid = -1;
        rid = -1;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip(ref_fasta_file);
    };

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
        //can be retained in the hash
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

            tid = bam_get_tid(s);
        }

        return true;
    }

    void genotype2()
    {
        if (mode=="d")
        {
            //iterate sam
            bam_hdr_t *h = odr->hdr;
            bam1_t * s = bam_init1();
            while (odr->read(s))
            {
                ++no_reads;

                if (!filter_read(s))
                {
                    continue;
                }

                gbr->flush(odw, h, s);
                gbr->process_read(h, s);

                ++no_passed_reads;
                if ((no_reads & 0x0000FFFF) == 0)
                {
                    std::cerr << bam_get_chrom(h,s) << ":" << bam_get_pos1(s) << " ("  << gbr->buffer.size() << ")\n";
                }
            }

            no_snps_genotyped = gbr->no_snps_genotyped;
            no_indels_genotyped = gbr->no_indels_genotyped;
            no_vntrs_genotyped = gbr->no_vntrs_genotyped;

            gbr->flush(odw, h, s, true);
            odw->close();
        }
        else if (mode=="s")
        {
            //iterate VCF file
                //random access per site
        }
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

    void print_options()
    {
        std::clog << "genotype2 v" << version << "\n\n";

        std::clog << "options:     input VCF File                       " << input_vcf_file << "\n";
        std::clog << "         [b] input BAM File                       " << input_sam_file << "\n";
        std::clog << "         [o] output VCF File                      " << output_vcf_file << "\n";
        std::clog << "         [s] sample ID                            " << sample_id << "\n";
        std::clog << "         [r] reference FASTA File                 " << ref_fasta_file << "\n";
        std::clog << "         [z] ignore MD tags                       " << (ignore_md ? "true": "false") << "\n";
        std::clog << "         [m] mode of genotyping                   " << mode << "\n";
        print_int_op("         [i] intervals                            ", intervals);
        std::clog << "\n";
        std::clog << "         [t] read mapping quality cutoff          " << read_mapq_cutoff << "\n";
        std::clog << "         [l] ignore overlapping read              " << (ignore_overlapping_read ? "true" : "false") << "\n";
        std::clog << "         [a] read flag filter                     " << std::showbase << std::hex << read_exclude_flag << std::dec << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "genotype2 v" << version << "\n\n";

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
        std::clog << "       no. SNPs genotyped           : " << no_snps_genotyped<< "\n";
        std::clog << "       no. Indels genotyped         : " << no_indels_genotyped << "\n";
        std::clog << "       no. VNTRs genotyped          : " << no_vntrs_genotyped << "\n";
        std::clog << "\n";
    }

    ~Igor()
    {
        kh_destroy(rdict, reads);
    };

    private:
};

}

void genotype2(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.genotype2();
    igor.print_stats();
}

