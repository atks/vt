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

#include "bam_variant_extractor.h"

/**
 * Constructor
 * baseq_cutoff - q value cutoff to select candidate SNPs
 */
BAMVariantExtractor::BAMVariantExtractor(size_t evidence_allele_count_cutoff,
              double fractional_evidence_allele_count_cutoff,
              size_t baseq_cutoff,
              std::string& ref_fasta_file)
: evidence_allele_count_cutoff(evidence_allele_count_cutoff),
  fractional_evidence_allele_count_cutoff(fractional_evidence_allele_count_cutoff),
  baseq_cutoff(baseq_cutoff)
{
    bcf1_t* v = bcf_init();
    
    alleles = {0,0,0};
    read_seq = {0,0,0};
    qual = {0,0,0};
    cigar = {0,0,0};
};

/**
 * Transfer read into a buffer for processing later
 */
void BAMVariantExtractor::process_read(bam_hdr_t *h, bam1_t *s)
{
    //extract relevant information from sam record
    char* chrom = bam_get_chrom(h, s);
    int32_t tid = bam_get_tid(s);
    int32_t pos0 = bam_get_pos0(s);
    
    //flush variants
    if (this->tid!=tid)
    {
        flush_variant_buffer();
    }
    else
    {
        this->chrom.assign(chrom);
        this->tid = tid;
        flush_variant_buffer();
        //reset buffer
    }

    int32_t ref_len;
    char* genome_seq = faidx_fetch_seq(fai, chrom, pos0-1, pos0+cigar.l, &ref_len);

    uint8_t* sq = bam_get_seq(s);
    bam_get_seq_string(s, &read_seq);
    
    size_t genome_pos0 = 0;
    
    const bam1_core_t *c = &s->core;
    

    flush_variant_buffer();
    if (ref_len>0) free(genome_seq);
};

/**
 * Empty variant buffer records that are completed.
 */
bool BAMVariantExtractor::flush_variant_buffer()
{
    return true;
};

/**
 * Processes buffer to pick up variant
 */
bool BAMVariantExtractor::next_variant(bcf1_t* v)
{
    return true;
};


/**
 * Checks if a variant is normalized.
 */
bool BAMVariantExtractor::is_biallelic_normalized(std::string& ref, std::string& alt)
{
    bool last_base_same = ref.at(ref.size()-1) == alt.at(alt.size()-1);
    bool exists_len_one_allele = ref.size()==1 || alt.size()==1;
    bool first_base_same = ref.at(0) == alt.at(0);

    if (last_base_same || (!exists_len_one_allele && first_base_same))
    {
        return false;
    }
    else
    {
        return true;
    }
}

/**
 * Normalize a biallelic variant.
 */
void BAMVariantExtractor::normalize_biallelic(size_t pos0, std::string& ref, std::string& alt)
{
    if (!is_biallelic_normalized(ref, alt))
    {
        bool to_right_trim = true;
        bool to_left_extend = false;

        while (to_right_trim || to_left_extend)
        {
            //checks if right trimmable or left extendable
            if (!ref.empty() && !alt.empty())
            {
                if (ref.at(ref.size()-1) == alt.at(alt.size()-1))
                {
                    to_right_trim = true;
                    to_left_extend = false;
                }

                if (pos0==0 && (ref.size()==1||alt.size()==1))
                {
                    to_right_trim = false;
                    to_left_extend = false;
                }
            }
            else
            {
                to_right_trim = false;
                to_left_extend = true;
            }

            if (to_right_trim)
            {
                ref.erase(ref.size()-1);
                alt.erase(alt.size()-1);
            }

            if (to_left_extend)
            {
                --pos0;
                int ref_len = 0;
                char *refseq = faidx_fetch_uc_seq(fai, chrom.c_str(), pos0, pos0, &ref_len);
                if (!refseq)
                {
                    fprintf(stderr, "[%s:%d %s] failure to extrac base from fasta file: %s:%zu: >\n", __FILE__, __LINE__, __FUNCTION__, chrom.c_str(), pos0);
                    exit(1);
                }
                char base = refseq[0];
                free(refseq);

                ref.insert(0, 1, base);
                alt.insert(0, 1, base);
            }
        }

        bool to_left_trim =  true;

        while (to_left_trim)
        {
            //checks if left trimmable.
            if (ref.size()==1 || ref.at(0)!=alt.at(0))
            {
                to_left_trim = false;
            }

            if (to_left_trim)
            {
                ref.erase(0, 1);
                alt.erase(0, 1);
                ++pos0;
            }
        }
    }
};