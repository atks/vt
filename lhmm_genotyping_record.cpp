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

#include "lhmm_genotyping_record.h"

typedef struct
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

/**
 * Maintains read information and allows for additional reads
 * till VCF record can be printed out.
 */
class LHMMGenotypingRecord : GenotypingRecord
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

    LHMMGenotypingRecord(bcf_hdr_t *h, bcf1_t *v)
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

    ~LHMMGenotypingRecord(){};

    /**
     * Prints records
     */
    void genotype(bam1_t *b)
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