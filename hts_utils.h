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

#ifndef HTS_UTILS_H
#define HTS_UTILS_H

#include <string>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <map>
#include <queue>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "utils.h"

/**********
 *FAI UTILS
 **********/
typedef struct {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
};

/**
 * An alternate sequence fetcher for upper case sequence.
 */
char *faidx_fetch_uc_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);

/**************
 *BAM HDR UTILS
 **************/

/**
 * Copies contigs found in bam header to bcf header.
 */
void bam_hdr_transfer_contigs_to_bcf_hdr(const bam_hdr_t *sh, bcf_hdr_t *vh);

/**
 * Get number of sequences.
 */
#define bam_hdr_get_n_targets(h) ((h)->n_targets)

/**
 * Get list of sequence names.
 */
#define bam_hdr_get_target_name(h) ((h)->target_name)

/**
 * Get list of sequence lenghts.
 */
#define bam_hdr_get_target_len(h) ((h)->target_len)

/**********
 *BAM UTILS
 **********/

//used when a base on a read is not aligned to the human genome reference
//in particular for soft clips and insertions
#define BAM_READ_INDEX_NA -1

/**
 * Gets the start position of the first mapped base in the sequence.
 */
#define bam_get_chrom(h, s) ((h)->target_name[(s)->core.tid])

/**
 * Gets the 1 based start position of the first mapped base in the read.
 */
#define bam_get_pos1(s) ((s)->core.pos + 1)

/**
 * Gets the 0 based start position of the first mapped base in the read.
 */
#define bam_get_pos0(s) ((s)->core.pos)

/**
 * Gets the end position of the last mapped base in the read.
 */
int32_t bam_get_end_pos1(bam1_t *srec);

/**
 * Gets the template ID of the paired read.
 */
#define bam_get_mtid(s) ((s)->core.mtid)

/**
 * Gets the start position of the first mapped base in the read.
 */
#define bam_get_mpos1(s) ((s)->core.mpos)

/**
 * Gets the read sequence from a bam record
 */
void bam_get_seq_string(bam1_t *s, kstring_t *seq);

/**
 * Gets the base qualities from a bam record
 */
void bam_get_qual_string(bam1_t *s, kstring_t *qual);

/**
 * Gets the cigar sequence from a bam record
 */
#define bam_get_n_cigar_op(b) ((b)->core.n_cigar)

/**
 * Gets the cigar from a BAM record
 */
void bam_get_cigar_string(bam1_t *s, kstring_t *cigar);

/**
 * Gets the cigar string from a bam record
 */
void bam_get_cigar_expanded_string(bam1_t *s, kstring_t *cigar_string);

/**
 * Is this sequence the first read?
 */
#define bam_is_fread1(s) (((s)->core.flag&BAM_FREAD1) != 0)

/**
 * Is this sequence the second read?
 */
#define bam_is_fread2(s) (((s)->core.flag&BAM_FREAD2) != 0)

/**
 * Gets the base in the read that is mapped to a genomic position.
 * Returns -1 if it does not exists.
 */
void bam_get_base_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos);

/**
 * Gets the base in the read that is mapped to a genomic position.
 * Returns true if successful, false if the location is a deletion or not on the read.
 */
void bam_get_base_and_qual_and_read_and_qual(bam1_t *s, uint32_t pos, char& base, char& qual, int32_t& rpos, kstring_t* readseq, kstring_t* readqual);

/**
 * Gets flag
 */
#define bam_get_flag(s) ((s)->core.flag)

/**
 * Get map quality
 */
#define bam_get_mapq(s) ((s)->core.qual)

/**
 *Get tid - e.g. chromosome id
 */
#define bam_get_tid(s) ((s)->core.tid)

/**
 * Get read length
 */
#define bam_get_l_qseq(s) ((s)->core.l_qseq)

/**************
 *BCF HDR UTILS
 **************/

/**
 * Copies contigs found in bcf header to another bcf header.
 */
void bcf_hdr_transfer_contigs(const bcf_hdr_t *sh, bcf_hdr_t *vh);

/**
 * Extracts sequence length by rid.
 */
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq);

/**
 * Get samples from bcf header
 */
#define bcf_hdr_get_samples(h) ((h)->samples)

/**
 * Get ith sample name from bcf header
 */
#define bcf_hdr_get_sample_name(h, i) ((h)->samples[i])

/**
 * Get number of samples in bcf header
 */
int32_t bcf_hdr_get_n_sample(bcf_hdr_t *h);

/**
 * Checks if a info header exists.
 */
bool bcf_hdr_info_exists(bcf_hdr_t *h, const char* key);

/**
 * Reads header of a VCF file and returns the bcf header object.
 * This wraps around vcf_hdr_read from the original htslib to
 * allow for an alternative header file to be read in.
 *
 * this searches for the alternative header saved as <filename>.hdr
 */
bcf_hdr_t *bcf_alt_hdr_read(htsFile *fp);

/**
 * Subsets a record by samples.
 * @imap - indices the subsetted samples
 */
int bcf_hdr_subset_samples(const bcf_hdr_t *h, bcf1_t *v, std::vector<int32_t>& imap);

/**********
 *BCF UTILS
 **********/

/**
 * Gets number of expected genotypes from number of allelles for a ploidy of 2.
 */
#define bcf_an2gn(n) (((n+1)*n)>>1)

/**
 * Gets a string representation of a variant.
 */
void bcf_variant2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Gets a sorted string representation of a variant.
 */
void bcf_variant2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Gets a string representation of the alleles of a variant.
 */
void bcf_alleles2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Gets a sorted string representation of the alleles of a variant.
 */
void bcf_alleles2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * Prints a VCF record to STDERR.
 */
void bcf_print(bcf_hdr_t *h, bcf1_t *v);

/**
 * Prints a VCF record in compact string representation to STDERR.
 */
void bcf_print_lite(bcf_hdr_t *h, bcf1_t *v);

/**
 * Prints a VCF record in compact string representation to STDERR with alleles sorted.
 */
void bcf_print_lite_sorted(bcf_hdr_t *h, bcf1_t *v);

/**
 * Prints a VCF record in compact string representation to STDERR with a new line.
 */
void bcf_print_liten(bcf_hdr_t *h, bcf1_t *v);

/**
 * Get number of chromosomes
 */
#define bcf_get_n_chrom(h) ((h)->n[BCF_DT_CTG])

/**
 * Get chromosome name
 */
#define bcf_get_chrom(h, v) ((h)->id[BCF_DT_CTG][(v)->rid].key)

/**
 * Set chromosome name
 */
void bcf_set_chrom(bcf_hdr_t *h, bcf1_t *v, const char* chrom);

/**
 * Get RID
 */
#define bcf_get_rid(v) ((v)->rid)

/**
 * Set RID
 */
#define bcf_set_rid(v, c) ((v)->rid=(c))

/**
 * Check if variant is passed
 */
bool bcf_is_passed(bcf_hdr_t *h, bcf1_t *v);

/**
 * Get 0-based position
 */
#define bcf_get_pos0(v) ((v)->pos)

/**
 * Set 0-based position
 */
#define bcf_set_pos0(v, p) ((v)->pos = (p))

/**
 * Get 1-based position
 */
#define bcf_get_pos1(v) ((v)->pos+1)

/**
 * Get 1-based end position
 */
#define bcf_get_end_pos1(v) ((v)->pos + strlen((v)->d.allele[0]))

/**
 * Set 1-based position
 */
#define bcf_set_pos1(v, p) ((v)->pos = (p)-1)

/**
 * Get id
 */
#define bcf_get_id(v) ((v)->d.id)

/**
 * Set id.
 */
void bcf_set_id(bcf1_t *v, char* id);

/**
 * Get allele
 */
#define bcf_get_allele(v) ((v)->d.allele)

/**
 * Get n_alleles of a bcf record
 */
#define bcf_get_n_allele(v) ((v)->n_allele)

/**
 * Get reference base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_ref(v) ((v)->d.allele[0][0])

/**
 * Get alternative base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_alt(v) ((v)->d.allele[1][0])

/**
 * Get reference allele
 */
#define bcf_get_ref(v) ((v)->d.allele[0])

/**
 * Get alternative allele
 */
#define bcf_get_alt(v, i) ((v)->d.allele[i])

/**
 * Get variant type
 */
#define bcf_get_var_type(v) ((v)->d.var_type)

/**
 * Get number of samples in bcf record
 */
#define bcf_get_n_sample(v) ((v)->n_sample)

/**
 * Set number of samples in bcf record
 */
#define bcf_set_n_sample(v, n) ((v)->n_sample = (n))

/**
 * Get qual
 */
#define bcf_get_qual(v) ((v)->qual)

/**
 * Set qual
 */
#define bcf_set_qual(v, q) ((v)->qual = (q))

#endif
