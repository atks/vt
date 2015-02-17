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

#include "discover2.h"

namespace
{

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
    uint32_t mapq_cutoff;
    uint32_t baseq_cutoff;

    uint32_t evidence_allele_count_cutoff;
    double fractional_evidence_allele_count_cutoff;

    uint16_t exclude_flag;

    khash_t(rdict) *reads;
    khiter_t k;
        int32_t ret;
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
    Pileup pileup;

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
            TCLAP::ValueArg<std::string> arg_input_bam_file("b", "b", "input BAM file", true, "", "string", cmd);

            cmd.parse(argc, argv);

            input_bam_file = arg_input_bam_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            output_vcf_file = arg_output_vcf_file.getValue();
            sample_id = arg_sample_id.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();
            mapq_cutoff = arg_mapq_cutoff.getValue();
            baseq_cutoff = arg_baseq_cutoff.getValue();
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

        odr = new BAMOrderedReader(input_bam_file, intervals, ref_fasta_file);
        s = bam_init1();

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bam_hdr_transfer_contigs_to_bcf_hdr(odr->hdr, odw->hdr);
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=E,Number=1,Type=Integer,Description=\"Number of reads containing evidence of the alternate allele\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Total number of reads at a candidate locus with reads that contain evidence of the alternate allele\">");
        bcf_hdr_add_sample(odw->hdr, sample_id.c_str());
        bcf_hdr_add_sample(odw->hdr, NULL);
        v = NULL;

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
        
    }

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
        char* md;
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
        std::cerr << "md       : " << md << "\n";


        if (seq.m) free(seq.s);
        if (qual.m) free(qual.s);
        if (cigar_string.m) free(cigar_string.s);
        if (cigar_expanded_string.m) free(cigar_expanded_string.s);
    }

    bool filter_read(bam1_t *s)
    {
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
            return false;
        }

        if (bam_get_mapq(s) < mapq_cutoff)
        {
            //filter short aligments and those with too many indels (?)
            ++no_low_mapq_reads;
            return false;
        }
        
        return true;
    }

    /**
     * Flush records out before gpos1.
     */
    void flush(bam1_t* s)
    {
        uint32_t gbeg1 = pileup.get_gbeg1();
        if (pileup.flushable(bam_get_tid(s), bam_get_pos1(s)))
        {
            for (uint32_t i=pileup.begin(); i!=pileup.end(); i=pileup.inc(i,1))
            {
                std::cerr << "\t\t\t" << gbeg1++ << "\n";
                pileup[i].print();
                pileup[i].clear();
            }
        }
    }

    /**
     * Flush records out before gpos1.
     */
    void flush()
    {
        uint32_t gbeg1 = pileup.get_gbeg1();
        for (uint32_t i=pileup.begin(); i!=pileup.end(); i=pileup.inc(i,1))
        {
            std::cerr << "\t\t\t" << gbeg1++ << "\n";
            pileup[i].print();
            pileup[i].clear();
        }
    }
    
    void process_read(bam1_t *s)
    {
        bam_print_key_values(odr->hdr, s);

        uint32_t tid = bam_get_tid(s);
        uint32_t pos1 = bam_get_pos1(s);
        uint8_t* seq = bam_get_seq(s);
        uint8_t* qual = bam_get_qual(s);
        int32_t l_qseq = bam_get_l_qseq(s);
        uint32_t* cigar = bam_get_cigar(s);

        uint8_t *aux;
        char* md;
        (aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(aux));

        //iterate cigar
        uint32_t n_cigar_op = bam_get_n_cigar_op(s);

        char* mdp = md;
        uint32_t cpos1 = pos1; //current 1 based genome position
        uint32_t spos0 = 0;    //current position in read sequence

        flush(s);

        if (n_cigar_op)
        {
            uint32_t *cigar = bam_get_cigar(s);
            for (uint32_t i = 0; i < n_cigar_op; ++i)
            {
                uint32_t oplen = bam_cigar_oplen(cigar[i]);
                char opchar = bam_cigar_opchr(cigar[i]);

                std::cerr << oplen << " " << opchar << "\n";

                if (opchar=='S')
                {
                    //add to S evidence
                    //do nothing

                    std::string ins = "";
                    float mean_qual = 0;
                    for (size_t i=0; i<oplen ; ++i)
                    {
                        ins += bam_base2char(bam_seqi(seq, spos0+i));
                        mean_qual += qual[spos0+i];
                    }

                    mean_qual /= oplen;

                    if (cpos1==pos1)
                    {
                        std::cerr << "\t\t\tadding LSCLIP: " << cpos1 << "\t" << ins << " (" << mean_qual << ")\n";
                        pileup.add_lsclip(cpos1, ins);    
                    }
                    else
                    {
                        std::cerr << "\t\t\tadding RSCLIP: " << cpos1 << "\t" << ins << " (" << mean_qual << ")\n";
                        pileup.add_rsclip(cpos1, ins); 
                    }

                    spos0 += oplen;
                }
                else if (opchar=='M')
                {
                    uint32_t lpos1 = cpos1; // we need this because M contains matches and mismatches
                    uint32_t sspos0 = spos0; // we need this because M contains matches and mismatches
                    uint32_t i = 0;
                    while (*mdp)
                    {
                        if (isdigit(*mdp)) //matches
                        {
                            char* end = 0;
                            int32_t len = std::strtol(mdp, &end, 10);
                            mdp = end;

                            std::cerr << "\tMatch " << len << "\n";

                            if (len)
                            {
                                uint32_t gbeg1 = lpos1;
                                uint32_t gend1 = lpos1+len-1;
                                std::cerr << "\t\t\tadding REF: " << gbeg1 << "-" << gend1 << ":";
                                for (size_t i=sspos0; i<=(sspos0+len-1); ++i)
                                {
                                    std::cerr << (bam_base2char(bam_seqi(seq, i)));
                                }
                                std::cerr << "\n";
                                std::cerr << "\t\t\t\t sspos0: " << sspos0  << " " << (sspos0+len-1) << "\n";
                                
                                pileup.add_ref(gbeg1, sspos0, len, seq, false);
                                
                                lpos1 += len;
                                sspos0 += len;
                            }
                        }
                        else if (isalpha(*mdp)) //SNPs
                        {
                            char ref = *mdp;
                            char alt = (bam_base2char(bam_seqi(seq, spos0+(lpos1-cpos1))));
                            std::cerr << "\tMismatch " << ref << "\n";
                            std::cerr << "\t\t\tadding SNP: " << lpos1 << ":" << ref << "/" << alt << "\n";
                            pileup.add_snp(lpos1, ref, alt, false);

                            ++lpos1;
                            ++mdp;
                            ++sspos0;
                        }
                        else
                        {
                            //unexpected MD tag value
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

                    if (*mdp!='^')
                    {
                        std::cerr << "inconsistent MD and cigar\n";
                        exit(1);
                    }
                    else
                    {
                        ++mdp;
                        std::string del = "";
                        while (isalpha(*mdp))
                        {
                            del += *mdp;
                            ++mdp;
                        }

                        std::cerr << "\t\t\tadding DEL: " << cpos1 << " " << del << "\n";
                        pileup.add_ins(cpos1, del);

                        cpos1 += oplen;
                    }
                }
                else if (opchar=='I')
                {
                    //insertions are not present in MD tags
                    //may be handled independently of future matches
                    std::string ins = "";
                    for (size_t i=0; i<oplen ; ++i)
                    {
                        ins += bam_base2char(bam_seqi(seq, spos0+i));
                    }

                    std::cerr << "\t\t\tadding INS: " << cpos1 << " " << ins  << "\n";
                    pileup.add_ins(cpos1, ins);
                    
                    spos0 += oplen;
                }
                else
                {
                    std::cerr << "never seen before state " << opchar << "\n";
                }
            }
        }
    }

    void discover()
    {
        odw->write_hdr();

        //for tracking overlapping reads
        reads = kh_init(rdict);
        
        while (odr->read(s))
        {
            ++no_reads;

            if (!filter_read(s)) 
            {
                continue;
            }

     
            process_read(s);

            ++no_passed_reads;

            //break;
        }

        //flush();

        odw->close();
    };

    void print_options()
    {
        std::clog << "discover v" << version << "\n\n";

        std::clog << "options: [b] input BAM File               " << input_bam_file << "\n";
        std::clog << "         [o] output VCF File              " << output_vcf_file << "\n";
        std::clog << "         [s] sample ID                    " << sample_id << "\n";
        std::clog << "         [r] reference FASTA File         " << ref_fasta_file << "\n";
        std::clog << "         [m] MAPQ cutoff                  " << mapq_cutoff << "\n";
        std::clog << "         [q] base quality cutoff          " << baseq_cutoff << "\n";
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
