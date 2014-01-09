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

namespace
{

typedef struct
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

/**
 * Maintains read information and allows for additional reads
 * till VCF record can be printed out.
 */
class VCFAuxRecord
{
    public:
    bcf1_t *v;
    int32_t vtype;
    uint32_t read_no;

    khash_t(rdict) *reads;
    khiter_t k;
    int32_t ret;

    char* ref_probe;
    char* alt_probe;
    uint32_t plen;

    kstring_t rqs;
    kstring_t aqs;

    VCFAuxRecord(bcf_hdr_t *h, bcf1_t *v)
    {
        reads = kh_init(rdict);
        
        this->v = v;
        
        ref_probe = (char*) malloc(1);
        
        int32_t len;
        int32_t ret = bcf_get_info_string(h, v, "REFPROBE", &ref_probe, &len);
        
        std::cerr << ref_probe << " " << len << " " << ret << "\n";
        bcf_print(h,v);
        
        
        exit(1);
        
//        int32_t probeLength = strlen(refProbe);
//        int32_t variantLengthDifference = (int32_t)strlen(alt)-(int32_t)strlen(ref);
        
    };

    ~VCFAuxRecord(){};

    /**
     * Prints records
     */
    void genotype_read(bam1_t *b)
    {
        if (vtype == VT_SNP)
        {
            //genotype SNP
        }
        else if (vtype & VT_INDEL)
        {
            //genotype Indel
        }
    }

    void genotype_indel(bam1_t* s)
    {
        //maximum depth cap
        if (read_no>255)
        {
            return;
        }

        //has secondary alignment or fail QC or is duplicate or is unmapped
        if (bam_get_flag(s) & 0x0704)
        {
            return;
        }

        //ignore poor quality mappings
        if (bam_get_mapq(s)<13)
        {
            return;
        }

        //this read is the first of the pair
        if (bam_get_mpos1(s) && (bam_get_tid(s)==bam_get_mtid(s)))
        {
            //first mate 
            if (bam_get_mpos1(s)>bam_get_pos1(s))
            {
                //overlapping (only insert a paired end if you know it will overlap)
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
                    }
                }
                return; // for the time being, just drop this read.
            }
        }




//
//        //read name
//
//        //sequence
//        readseq.l=0;
//        bam_get_seq_string(s, &readseq);
//
//        //qual
//        bam_get_qual_string(s, &readqual, readseq.s);
//
//        //map qual
//        uint32_t mapQual = bam_get_mapq(s);
//
//        std::string refCigar = "";
//        std::string altCigar = "";
//
//        double refllk = 0, altllk = 0;
//
//        lhmm_ref.align(refllk, refProbe, readseq.s, readqual.s);
//        lhmm_ref.computeLogLikelihood(refllk, lhmm_ref.getPath(), readqual.s);
//
//        lhmm_alt.align(altllk, altProbe, readseq.s, readqual.s);
//        lhmm_alt.computeLogLikelihood(altllk, lhmm_alt.getPath(), readqual.s);
//
//        std::string pad = "\t";
//        if (debug)
//        {
//            prec.log << pad << "==================\n";
//            prec.log << pad <<  (prec.read_no+1) << ") " << chrom << ":" << (pos0+1) << ":" << ref << ":" << alt << ":" << variantLengthDifference << "\n";
//            prec.log << pad <<  read_name << "\n";
//            prec.log << pad << "==================\n";
//            prec.log << pad << "ref probe     " << refProbe << " (" << plen << "/" << probeLength << ")\n";
//            prec.log << pad << "read sequence " << readseq.s  << "\n";
//            prec.log << pad << "==================\n";
//            lhmm_ref.printAlignment(pad, prec.log);
//            prec.log << pad << "ref llk: " << refllk << "\n";
//            prec.log << pad << "expected indel location: " << plen+1 << "\n";
//            prec.log << pad << "==================\n";
//            prec.log << pad << "==================\n";
//            prec.log << pad << "alt probe     " << altProbe << " (" << plen << ")\n";
//            prec.log << pad << "read sequence " << readseq.s  << "\n";
//            prec.log << pad << "==================\n";
//            lhmm_alt.printAlignment(pad, prec.log);
//            prec.log << pad << "alt llk: " << altllk << "\n";
//            prec.log << pad << "expected indel location: " << plen+1 << "\n";
//            prec.log << pad << "==================\n\n";
//        }
//
//        //check if the Insertions and Deletions are at the expected places.
//        //deletion
//        uint32_t ref_rpos=0, alt_rpos=0;
//        if (variantLengthDifference<0)
//        {
//            //deletion
//            if (!deletion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                 insertion_start_exists(lhmm_alt, plen+1, alt_rpos))
//            {
//                if (refllk<altllk)
//                {
//                    swap(refllk, altllk);
//                }
//
//                //ref allele
//                if (debug)
//                    std::cerr << pad << "DEL: REF ALLELE\n";
//            }
//            else if (deletion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                    !insertion_start_exists(lhmm_alt, plen+1, alt_rpos))
//            {
//                if (refllk>altllk)
//                {
//                    swap(refllk, altllk);
//                }
//
//                //alt allele
//                if (debug)
//                    std::cerr << pad << "DEL: ALT ALLELE\n";
//            }
//            else
//            {
//                refllk = 0;
//                altllk = 0;
//                if (debug)
//                    std::cerr << pad << "Deletion not at expected location, set as ambiguous\n";
//            }
//        }
//        else if (variantLengthDifference>0)
//        {
//            //insertion
//            if (!insertion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                 deletion_start_exists(lhmm_alt, plen+1, alt_rpos))
//            {
//                if (refllk<altllk)
//                {
//                    swap(refllk, altllk);
//                }
//
//                //ref allele
//                if (debug)
//                    std::cerr << pad << "INS: REF ALLELE\n";
//            }
//            else if (insertion_start_exists(lhmm_ref, plen+1, ref_rpos) &&
//                    !deletion_start_exists(lhmm_alt, plen+1, alt_rpos))
//            {
//                if (refllk>altllk)
//                {
//                    swap(refllk, altllk);
//                }
//
//                //alt allele
//                if (debug)
//                    std::cerr << pad << "INS: ALT ALLELE\n";
//            }
//            else
//            {
//                refllk = 0;
//                altllk = 0;
//                if (debug)
//                    std::cerr << pad << "Insertion not at expected location "  << plen+1 << ", set as ambiguous\n";
//            }
//        }
//
//        if (debug)
//        {
//          std::cerr << pad << "++++++++++++++++++\n";
//            std::cerr << pad << "reflk " << refllk << "\n";
//          std::cerr << pad << "altlk " << altllk << "\n";
//        }
//
//        uint32_t baseqr = 0, baseqa = 0;
//        //char allele = 'N';
//        uint32_t rpos = 0;
//        if (refllk>altllk)
//        {
//          //allele = 'R';
//            rpos = alt_rpos;
//            baseqa = round(-10*(altllk-refllk));
//        }
//        else if (refllk<altllk)
//        {
//          //allele = 'A';
//              rpos = ref_rpos;
//              baseqr = round(-10*(refllk-altllk));
//        }
//
//        prec.rqs << (uint8_t)(baseqr>93? 126 : baseqr+33);
//        prec.aqs << (uint8_t)(baseqa>93? 126 : baseqa+33);
//        prec.mqs << (uint8_t)(mapQual>93? 126 : mapQual+33);
//        int32_t tp =  pack_tp(!bam_is_rev(s), bam_is_fread1(s), rpos, bam_get_l_qseq(s));
//        prec.tps << std::hex << (tp<16?"0":"") << tp << std::dec ;
//
//        ++prec.read_no;
//
//        free(refProbe);
//        free(altProbe);
//        free(plenstring);
    }

    /**
     * Prints records
     */
    void print(BCFOrderedWriter *odw)
    {

    }

    /**
     * Prints records
     */
    void clear()
    {
    }
};

class VCFGenotypingPool
{
    public:
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    std::list<VCFAuxRecord*> buffer;
    std::list<VCFAuxRecord*> pool;

    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;

    uint32_t bref, vref;
    uint32_t bstart, bend, vpos;

    VCFGenotypingPool()
    {
        no_snps_genotyped = 0;
        no_indels_genotyped = 0;
    }
 
    void flush()
    {
        //flush out remaining records
        std::list<VCFAuxRecord*>::iterator prec_iter = buffer.begin();
        while (prec_iter!=buffer.end())
        {
            VCFAuxRecord* prec = *prec_iter;
            bcf1_t *v = prec->v;

            if (bcf_get_var_type(v) == VCF_SNP)
            {
                //prec->print(OUT_GLF);
                ++no_snps_genotyped;
            }
            else if (bcf_get_var_type(v) == VCF_INDEL || bcf_get_var_type(v) == VCF_OTHER)
            {
                //prec->print(OUT_GLF);
                ++no_indels_genotyped;
            }

            //move record to pool
            prec->clear();
            pool.push_front(prec);
            prec_iter = buffer.erase(prec_iter);
        }
    }

    /**
     * Adds pileup records till the first record after epos1.
     * Coordinates all 1 based.
     */
    bool add_prec(int32_t epos1)
    {
        //add records only if the new record overlaps with read
        if (buffer.size()!=0)
        {
            bcf1_t *v = buffer.back()->v;
            int32_t vpos1 = bcf_get_pos1(v);

            if (vpos1>epos1)
            {
                return false;
            }
        }

        bool added_record = false;
        bcf1_t *v = odw->get_bcf1_from_pool();
        while (odr->read(v))
        {
            VCFAuxRecord *p = NULL;
            if (pool.size()!=0)
            {
                p = pool.front();
                pool.pop_front();
            }
            else
            {
                p = new VCFAuxRecord(odr->hdr, v);
            }

            p->clear();

            //bcf_set_variant_types(p->v);

            buffer.push_back(p);
            added_record = true;

            if ((bcf_get_pos1(p->v))>epos1)
                break;
        }

        return added_record;
    }

    int32_t read_is_before_first_vcf_record(bam1_t *s)
    {
        bcf1_t *v = (*(buffer.begin()))->v;
        int32_t vpos1 = bcf_get_pos1(v);
        int32_t epos1 = bam_get_end_pos1(s);

        return (epos1+100000)<vpos1 ? vpos1 : 0;
    }

//    void genotype_snp(VCFAuxRecord &prec, bam1_t* s)
//    {
//        //maximum depth cap
//        if (prec.read_no>255)
//        {
//            return;
//        }
//
//        //has secondary alignment or fail QC or is duplicate or is unmapped
//        if (bam_get_flag(s) & 0x0704)
//        {
//            return;
//        }
//
//        //make this setable
//        //ignore poor quality mappings
//        if(bam_get_mapq(s)<20)
//        {
//            return;
//        }
//
//        bcf1_t *v = prec.v;
//        const char* chrom = bcf_get_chrom(ivcf_hdr, v);
//        uint32_t pos0 = bcf_get_pos0(v);
//        char ref = bcf_get_snp_ref(v);
//        char alt = bcf_get_snp_alt(v);
//        char* read_name = bam_get_qname(s);
//
//        //handle mate pairs
//      if(prec.read_ids.find(read_name)==prec.read_ids.end())
//      {
//          prec.read_ids[read_name] = 1;
//      }
//      else
//      {
//          prec.read_ids.erase(read_name);
//          return;
//      }
//
//      //get base and qual
//      char base, qual; int32_t rpos;
//        bam_get_base_and_qual(s, v->pos, base, qual, rpos);
//
//        //strand
//        char strand = bam_is_rev(s) ? 'R' : 'F';
//
//        //map qual
//        uint32_t mapQual = bam_get_mapq(s);
//
//        //fail to find the expected mapped base on the read
//        if (rpos==BAM_READ_INDEX_NA)
//        {
//            return;
//        }
//
//        //ignore ambiguous bases
//        if (base=='N')
//        {
//            return;
//        }
//
//      //////////////////////////////////////////////////
//      //perform GL computation here, uses map alignments
//      //////////////////////////////////////////////////
//        std::string refCigar = "";
//        std::string altCigar = "";
//
//        //compute genotype likelihood using alignment coordinates
//        double refllk = log10Emission(ref, base, qual);
//        double altllk = log10Emission(alt, base, qual);
//
//        ////////////////////////////////
//        //aggregate genotype likelihoods
//        ////////////////////////////////
//        uint32_t baseqr = 0, baseqa = 0;
//
//        if (refllk>altllk)
//        {
//            baseqa = round(-10*(altllk-refllk));
//        }
//        else if (refllk<altllk)
//        {
//            baseqr = round(-10*(refllk-altllk));
//        }
//
//        if (debug)
//        {
//            kstring_t cigar;
//            cigar.l = cigar.m = 0; cigar.s = 0;
//            bam_get_cigar_string(s, &cigar);
//
//            kstring_t aligned_read;
//          aligned_read.l = aligned_read.m = 0;
//          aligned_read.s = 0;
//
//          kstring_t aligned_qual;
//          aligned_qual.l = aligned_qual.m = 0;
//          aligned_qual.s = 0;
//
//          kstring_t expanded_cigar;
//          expanded_cigar.l = expanded_cigar.m = 0;
//          expanded_cigar.s = 0;
//
//          kstring_t annotations;
//          annotations.l = annotations.m = 0;
//          annotations.s = 0;
//
//            bam_get_seq_string(s, &readseq);
//            bam_get_qual_string(s, &readqual, readseq.s);
//
//            char allele = strand=='F' ? 'N' : 'n';
//            if (refllk>altllk)
//            {
//              allele = strand=='F' ? 'R' : 'r';
//            }
//            else if (refllk<altllk)
//            {
//              allele = strand=='F' ? 'A' : 'a';
//            }
//
//            std::string pad = "\t";
//          prec.log << pad << "==================\n";
//          prec.log << pad <<  (prec.read_no+1) << ") " << chrom << ":" << (pos0+1) << ":" << ref << ":" << alt << "\n";
//            prec.log << pad << bam_get_qname(s) << "\n";
//            prec.log << pad << "==================\n";
//            print_read_alignment(readseq, readqual, cigar, aligned_read, aligned_qual, expanded_cigar, annotations, rpos);
//
//            prec.log << pad << "read  " << aligned_read.s  << "\n";
//            prec.log << pad << "qual  " << aligned_qual.s  << "\n";
//            prec.log << pad << "cigar " << expanded_cigar.s  << "\n";
//            prec.log << pad << "anno  " << annotations.s  << "\n";
//            prec.log << pad << "==================\n";
//            prec.log << pad << "base   " << base  << "\n";
//            prec.log << pad << "qual   " << (int32_t)(qual-33)  << "\n";
//            prec.log << pad << "rpos   " << rpos  << "\n";
//          prec.log << pad << "refllk " << refllk << "\n";
//          prec.log << pad << "altllk " << altllk << "\n";
//          prec.log << pad << "allele " << allele << "\n\n";
//      }
//
//        prec.rqs << (uint8_t) (baseqr>93 ? 126 : baseqr+33);
//      prec.aqs << (uint8_t) (baseqa>93 ? 126 : baseqa+33);
//      prec.mqs << (uint8_t) (mapQual>93 ? 126 : mapQual+33);
//      int32_t tp =  pack_tp(!bam_is_rev(s), bam_is_fread1(s), rpos, bam_get_l_qseq(s));
//      prec.tps << std::hex << (tp<16?"0":"") << tp << std::dec ;
//
//      ++prec.read_no;
//    }



//    double log10Emission(char readBase, char probeBase, uint32_t qual)
//    {
//        double e = lt.pl2prob(qual);
//
//        if (readBase=='N' || probeBase=='N')
//        {
//            return 0;
//        }
//
//        return readBase!=probeBase ? log10(e/3) : log10(1-e);
//    };

//    void generateCigarString(kstring_t* expanded_cigar, kstring_t* cigar)
//    {
//        expanded_cigar->l = 0;
//        int32_t i=0, lastIndex = cigar->l-1;
//        std::stringstream token;
//
//        if (lastIndex<0)
//        {
//            return;
//        }
//        char c;
//        bool seenM = false;
//
//        while (i<=lastIndex)
//        {
//            c = cigar->s[i];
//
//            //captures the count
//            if (c<'A')
//            {
//                token << c;
//            }
//
//            if (c>'A' ||
//                i==lastIndex)
//            {
//                uint32_t count;
//                std::string s = token.str();
//                str2uint32(s, count);
//
//                //it is possible for I's to be observed before the first M's in the cigar string
//                //in this case, we treat them as 'S'
//                if (!seenM)
//                {
//                    if (c=='I')
//                    {
//                        c = 'S';
//                    }
//                    else if (c=='M')
//                    {
//                        seenM = true;
//                    }
//                }
//
//                for (uint32_t j=0; j<count; ++j)
//                    kputc_(c, expanded_cigar);
//                token.str("");
//            }
//
//            ++i;
//        }
//
//        kputc_(0, expanded_cigar);
//    };

//    void print_read_alignment(kstring_t& read, kstring_t& qual, kstring_t& cigar,
//                              kstring_t& aligned_read, kstring_t& aligned_qual, kstring_t& expanded_cigar, kstring_t& annotations,
//                              uint32_t rpos)
//    {
//        aligned_read.l = 0;
//        aligned_qual.l = 0;
//        expanded_cigar.l = 0;
//        generateCigarString(&expanded_cigar, &cigar);
//        annotations.l = 0;
//
//        uint32_t j=0;
//        for (uint32_t i=0; i<expanded_cigar.l; ++i)
//        {
//            char state = expanded_cigar.s[i];
//            if (state=='M' || state=='S' || state=='I')
//            {
//                kputc(read.s[j], &aligned_read);
//                kputc(qual.s[j], &aligned_qual);
//
//                if (j==rpos)
//                {
//                    kputc('^', &annotations);
//                }
//                else
//                {
//                    kputc(' ', &annotations);
//                }
//                ++j;
//            }
//            else if (state=='D')
//            {
//                kputc(' ', &aligned_read);
//                kputc(' ', &aligned_qual);
//            }
//        }
//    }

};

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
    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;

    /////////
    //tools//
    /////////
    LogTool lt;
    LHMM lhmm_ref, lhmm_alt;

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
            TCLAP::SwitchArg arg_iterate_by_site("c", "c", "iterate by candidate sites", cmd, false);
            TCLAP::SwitchArg arg_debug("d", "d", "debug alignments", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            input_sam_file = arg_input_sam_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
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
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////
        vodr = new BCFOrderedReader(input_vcf_file, intervals);
        sodr = new BAMOrderedReader(input_sam_file, intervals);
        vodw = new BCFOrderedWriter(output_vcf_file, 0);

        

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
        //works well for sparse number of candidate variants
        if (iterate_by_site)
        {
            
            
            bcf1_t *v = vodw->get_bcf1_from_pool();
            while (vodr->read(v))
            {
                VCFAuxRecord record(vodr->hdr, v);
                
                //construct interval from vcf record
                std::string chrom(bcf_get_chrom(vodr->hdr, v));
                GenomeInterval interval(chrom, bcf_get_pos1(v), bcf_get_pos1(v));
                sodr->jump_to_interval(interval);
                while(sodr->read(s))
                {
                    
                }
            }
        }
        //works well for dense set of candidate variants
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
