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

LHMMGenotypingRecord::LHMMGenotypingRecord(faidx_t *fai):GenotypingRecord(fai)
{
    this->fai = fai;
};

LHMMGenotypingRecord::LHMMGenotypingRecord(bcf_hdr_t *h, bcf1_t *v, faidx_t *fai):GenotypingRecord(h,v,fai)
{
    initialize(h, v, fai);
};
    
LHMMGenotypingRecord::~LHMMGenotypingRecord(){};

/**
 * Initializes a candidate variant for genotyping.
 */
bool LHMMGenotypingRecord::initialize(bcf_hdr_t *h, bcf1_t *v, faidx_t *fai)
{
    this->h = h;
    this->v = v;
    this->fai = fai;
    
    return set(v);
}
    
/**
 * Initializes a candidate VCF record. Returns false if failure.
 */
bool LHMMGenotypingRecord::set(bcf1_t *v)
{
    reads = kh_init(rdict);
    
    std::cerr << "working???\n";
    
    int32_t ret1,ret2,ret3,len;            
    ret1 = bcf_get_info_string(h, v, "REFPROBE", &ref_probe, &len);
    ret2 = bcf_get_info_string(h, v, "ALTPROBE", &alt_probe, &len);
    ret3 = bcf_get_info_int32(h, v, "PLEN", &plen, &len);
    
    if (ret1<0 || ret2<0 || ret3<0)
    {
        if (fai)
        {
             //ignore alleles with N
            if (strchr(bcf_get_alt(v, 0), 'N') || strchr(bcf_get_alt(v, 1), 'N'))
            {
                error_msg = "Alleles have N.";
                return false;
            }

            std::vector<std::string> probes;
            std::vector<std::string> alleles;

            int32_t preambleLength = 0;
            
            //populate alleles
            for (int32_t i=0; i<bcf_get_n_allele(v); ++i)
            {
                alleles.push_back(std::string(bcf_get_alt(v,i)));
            }
            
            generate_probes(bcf_get_chrom(h,v), bcf_get_pos1(v), 1, alleles, probes, 20, preambleLength);

            //remove ill defined probes
            bool skip = false;
            for (size_t i=1; i<probes.size()-1; ++i)
            {
                if(strchr(probes[i].c_str(), 'N'))
                {
                    skip = true;
                }
            }

            if (skip)
            {
                error_msg = "Probes contain N.";
                return false;
            }    
            
            return true;
        }
        else
        {
            fprintf(stderr, "[E %s:%d %s] Probe information appears to be missing, cannot proceed unless reference FASTA file is available\n", __FILE__, __LINE__, __FUNCTION__);
            exit(1);
        }
    }            
    return true;
}
    
/**
 * Prints records
 */
void LHMMGenotypingRecord::genotype(bam1_t *b)
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

void LHMMGenotypingRecord::genotype_indel(bam1_t* s)
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
void LHMMGenotypingRecord::print(BCFOrderedWriter *odw)
{

}

/**
 * Prints records
 */
void LHMMGenotypingRecord::clear()
{
}

/**
 * Generates a probing haplotype with flanks around the variant of interest.
 * Flanks are equal length
 */
void LHMMGenotypingRecord::generate_probes(const char* chrom,
                        int32_t pos1, 
                        uint32_t probeDiff, //
                        std::vector<std::string>& alleles, //store alleles
                        std::vector<std::string>& probes, //store probes
                        uint32_t min_flank_length,
                        int32_t& preambleLength) //store preamble length
{
    //map to usable number
    probes.resize(alleles.size(), "");

    //check allele lengths
    std::map<uint32_t, uint32_t> alleleLengths;
    for (uint32_t i=0; i<alleles.size(); ++i)
    {
       alleleLengths[alleles[i].size()]=1;
    }

    //for SNPs and MNPs and block substitutions
    if (alleleLengths.size()==1)
    {
        //just get flanking sequences
        //append preamble
        std::map<char, uint32_t> bases;
        std::string preamble;
        std::string postamble;
        char *base;
        uint32_t i = 1;
        int32_t ref_len;
        while (bases.size()<4 || preamble.size()<min_flank_length)
        {
            base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), pos1-1, pos1-1, &ref_len);
            preamble.append(1,base[0]);
            bases[base[0]] = 1;
            if (ref_len>0) free(base);
            ++i;
        }

        bases.clear();
        i=0;
        uint32_t alleleLength = alleles[0].size();
        while (bases.size()<4 || postamble.size()<min_flank_length)
        {
            base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), pos1+alleleLength+i, pos1+alleleLength+i, &ref_len);
            postamble.append(1,base[0]);
            bases[base[0]] = 1;
            if (ref_len>0) free(base);
            ++i;
        }

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            probes[i] = std::string(preamble.rbegin(), preamble.rend()).append(alleles[i]);
            probes[i] = probes[i].append(postamble);
        }

        preambleLength = preamble.size();
    }
    //for Indels and Complex Substitutions
    else
    {
        //find gald
        uint32_t min_len = alleles[0].size();
        uint32_t max_len = alleles[0].size();
        for (uint32_t i=1; i<alleles.size(); ++i)
        {
            if (alleles[i].size()<min_len)
                min_len = alleles[i].size();
            if (alleles[i].size()>max_len)
                max_len = alleles[i].size();
        }
        uint32_t gald = max_len-min_len;

        uint32_t currentDiff = 0;
        //current length of probe
        uint32_t length = 0;
        //number of point differences for each probe wrt the reference
        std::vector<uint32_t> diff(alleles.size(), 0);
        probes.resize(alleles.size(), "");

        generate_probes(chrom, pos1, min_flank_length, currentDiff, length, gald, diff, alleles, probes);

        //append preamble
        std::map<char, uint32_t> bases;
        std::string preamble;
        char* base;
        uint32_t i = 1;
        int32_t ref_len = 0;
        while (bases.size()<4 && preamble.size()<min_flank_length)
        {
            base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), pos1+i-1, pos1+i-1, &ref_len);
            preamble.append(1,base[0]);
            bases[base[0]] = 1;
            ++i;
            if (base[0]=='N')
            {
                break;
            }
            if (ref_len>0) free(base);
        }

        preambleLength = preamble.size();

        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            probes[i] = std::string(preamble.rbegin(), preamble.rend()).append(probes[i]);
        }
    }
}

/**
 * Iteratively called function for generating a haplotype.
 */
void LHMMGenotypingRecord::generate_probes(const char* chrom,
                        int32_t pos1,
                        uint32_t flankLength,
                        uint32_t& currentDiff,
                        uint32_t& length,
                        uint32_t gald,
                        std::vector<uint32_t>& diff,
                        std::vector<std::string>& alleles,
                        std::vector<std::string>& probes)
{
    if (currentDiff<alleles.size() || length<=2*gald+flankLength)
    {
        std::map<std::string, uint32_t> probeHash;
        //extend probes
        for (uint32_t i=0; i<alleles.size(); ++i)
        {
            {
                //copy from allele
                if (length<alleles[i].size())
                {
                    probes[i].append(1,alleles[i].at(length));
                }
                else//copy from reference
                {
                    int32_t start1 = (pos1+length-alleles[i].size()+alleles[0].size()-1);
                    int32_t ref_len;
                    char* base = faidx_fetch_uc_seq(fai, const_cast<char*>(chrom), start1 , start1, &ref_len);
                    probes[i].append(1, base[0]);
                    if (ref_len>0) free(base);
                }
            }
            probeHash[probes[i]] = 1;
        }

        currentDiff = probeHash.size();
        ++length;
        //std::cerr << probes[0] << "\n" << probes[1] << "\n";
        generate_probes(chrom, pos1, flankLength, currentDiff, length, gald, diff, alleles, probes);
    }
}
