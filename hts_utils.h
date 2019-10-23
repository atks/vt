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

#include <sys/stat.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "htslib/hfile.h"
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

/**
 * An alternate sequence fetcher for upper case sequence.
 */
char *faidx_fetch_uc_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);

/**********
 *HTS UTILS
 **********/

/**
 * Checks file extension for use in writing files.
 */
bool str_ends_with(std::string const & value, std::string const & ending);

/**
 * Checks file extension based on filename.
 */
int32_t hts_filename_type(std::string const& value);

/**************
 *BAM HDR UTILS
 **************/

/**
 * Copies contigs found in bam header to bcf header.
 */
void bam_hdr_transfer_contigs_to_bcf_hdr(const bam_hdr_t *sh, bcf_hdr_t *vh);

/**
 * Gets sample names from bam header.
 */
std::string bam_hdr_get_sample_name(bam_hdr_t* hdr);

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
 * Gets the chromosome name of the tid.
 */
#define bam_get_chrom(h, s) ((h)->target_name[(s)->core.tid])

/**
 * Gets the 1 based start position of the first mapped base in the read.
 */
#define bam_get_pos1(s) ((int32_t)(s)->core.pos + 1)

/**
 * Gets the 0 based start position of the first mapped base in the read.
 */
#define bam_get_pos0(s) ((int32_t)(s)->core.pos)

/**
 * Gets the end position of the last mapped base in the read.
 */
int32_t bam_get_end_pos1(bam1_t *srec);

/**
 * Gets the template ID of the read.
 */
#define bam_get_tid(s) ((s)->core.tid)

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
 * Gets the base qualities from a bam record.
 */
void bam_get_qual_string(bam1_t *s, kstring_t *qual);

/**
 * Gets the cigar sequence from a bam record.
 */
#define bam_get_n_cigar_op(b) ((b)->core.n_cigar)

/**
 * Gets the cigar from a BAM record.
 */
void bam_get_cigar_string(bam1_t *s, kstring_t *cigar);

/**
 * Gets the cigar string from a bam record.
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
 * Converts base encoding to character.
 */
#define bam_base2char(b) ("XACXGXXXTXXXXXXN"[(b)])

/**
 * Gets flag.
 */
#define bam_get_flag(s) ((s)->core.flag)

/**
 * Get map quality.
 */
#define bam_get_mapq(s) ((s)->core.qual)

/**
 * Get tid - e.g. chromosome id.
 */
#define bam_get_tid(s) ((s)->core.tid)

/**
 * Get read length.
 */
#define bam_get_l_qseq(s) ((s)->core.l_qseq)

/**
 * Prints a bam.
 */
void bam_print(bam_hdr_t *h, bam1_t *s);

/**************
 *BCF HDR UTILS
 **************/

/**
 * Prints header.
 */
void bcf_hdr_print(bcf_hdr_t *h);

/**
 * Copies contigs found in bcf header to another bcf header.
 */
void bcf_hdr_transfer_contigs(const bcf_hdr_t *sh, bcf_hdr_t *vh);

/**
 * Checks if a particular header type exists
 * @hdr  - header
 * @type - BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT, BCF_HL_CTG
 * @key  - the key name
 */
bool bcf_hdr_exists(bcf_hdr_t *hdr, int32_t type, const char *key);

/**
 * Extracts sequence length by rid.
 */
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq);

/**
 * Copies an info fields from one record to another
 * @hsrc  - source header
 * @vsrc  - source bcf1_t
 * @hdest - destination header
 * @vdest - destination bcf1_t
 * @key   - the key name
 * @type  - BCF_HT_FLAG, BCF_HT_INT, BCF_HT_REAL, BCF_HT_STR
 */
void bcf_copy_info_field(bcf_hdr_t *hsrc, bcf1_t* vsrc, bcf_hdr_t *hdest, bcf1_t* vdest, const char* key, int32_t type);

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
 * Reads header of a VCF file and returns the bcf header object.
 * This wraps around vcf_hdr_read from the original htslib to
 * allow for an alternative header file to be read in.
 *
 * this searches for the alternative header saved as <filename>.hdr
 * If the VCF files is BCF, any alternative header is ignored.
 */
bcf_hdr_t *bcf_alt_hdr_read(htsFile *fp);

/**
 * Subsets a record by samples.
 * @imap - indices the subsetted samples
 */
int bcf_hdr_subset_samples(const bcf_hdr_t *h, bcf1_t *v, std::vector<int32_t>& imap);

/**
 * Help function for adding a header with a backup tag name.
 * If the <tag> is already present, a new tag is attempted
 * in the format <tag>_1 to <tag>_9.  If <tag>_9 failed,
 * the function will not add any new tag and will return
 * an empty string.
 *
 * Returns the tag that was inserted or updated.
 */
std::string bcf_hdr_append_info_with_backup_naming(bcf_hdr_t *h, std::string tag, std::string number, std::string type, std::string description, bool rename);

/**
 * Translates BCF_VL types to string.
 */
std::string bcf_hdr_vl2str(int32_t id);

/**
 * Translates BCF_BT types to string.
 */
std::string bcf_hdr_ht2str(int32_t id);

/**********
 *BCF UTILS
 **********/

/**
 * Gets number of expected genotypes from number of alleles for a ploidy of 2.
 */
#define bcf_an2gn(n) (((n+1)*n)>>1)

/**
 * Gets number of genotypes from number of alleles and ploidy.
 *
 * Returns 0 if number of alleles and genotypes are not consistent.
 */
uint32_t bcf_ap2g(uint32_t no_allele, uint32_t no_ploidy);

/**
 * Gets alleles from number of ploidy and genotype index.
 */
void bcf_pg2a(uint32_t no_ploidy, uint32_t genotype_index, std::vector<int32_t>& alleles);

/**
 * Gets number of ploidy from number of alleles and genotypes.
 */
uint32_t bcf_ag2p(uint32_t no_alleles, uint32_t no_genotypes);

/**
 * Gets genotype from genotype index and ploidy.
 */
std::vector<int32_t> bcf_ip2g(int32_t genotype_index, uint32_t no_ploidy);

/**
 * Gets index of a genotype of p ploidy.
 */
uint32_t bcf_g2i(int32_t* g, uint32_t p);

/**
 * Gets index of a genotype of n ploidy.
 */
uint32_t bcf_g2i(std::vector<int32_t>& g);

/**
 * Gets index of a diploid genotype.
 */
#define bcf_2g2c(i,j) ((i)+((((j)+1)*(j))>>1))

///**
// * Gets number of genotypes from number of alleles and ploidy.
// */
//uint32_t bcf_g2i(std::string genotype);

/**
 * Returns true if a is before b, false otherwise.
 */
bool bcf_is_in_order(bcf1_t *a, bcf1_t *b);

/**
 * Returns a copy v that only has the chr:pos1:ref:alt information.
 */
bcf1_t* bcf_copy_variant(bcf_hdr_t *h, bcf1_t *v);

/**
 * Gets a string representation of a variant.
 */
std::string bcf_variant2string(bcf_hdr_t *h, bcf1_t *v);

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
const char* bcf_get_chrom(bcf_hdr_t *h, bcf1_t *v);

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
#define bcf_get_pos0(v) ((int32_t)(v)->pos)

/**
 * Set 0-based position
 */
#define bcf_set_pos0(v, p) ((v)->pos = (p))

/**
 * Get 1-based position
 */
#define bcf_get_pos1(v) ((int32_t)(v)->pos+1)

/**
 * Get 1-based end position
 */
#define bcf_get_end1(v) ((int32_t)(v)->pos + strlen((v)->d.allele[0]))

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
 * Gets an info flag.
 */
bool bcf_get_info_flg(bcf_hdr_t *h, bcf1_t *v, const char* tag);

/**
 * Sets an info flag.
 */
void bcf_set_info_flg(bcf_hdr_t *h, bcf1_t *v, const char* tag, bool value);

/**
 * Gets an info integer.
 */
int32_t bcf_get_info_int(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t default_value = 0);

/**
 * Sets an info integer.
 */
void bcf_set_info_int(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t value);

/**
 * Gets an info integer vector.
 */
std::vector<int32_t> bcf_get_info_int_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t default_size=0, int32_t default_value=0);

/**
 * Sets an info integer vector.
 */
void bcf_set_info_int_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::vector<int32_t>& values);

/**
 * Gets an info float.
 */
float bcf_get_info_flt(bcf_hdr_t *h, bcf1_t *v, const char* tag, float default_value = 0);

/**
 * Sets an info float.
 */
void bcf_set_info_flt(bcf_hdr_t *h, bcf1_t *v, const char* tag, float value);

/**
 * Gets an info integer vector.
 */
std::vector<float> bcf_get_info_flt_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t default_size=0, float default_value=0);

/**
 * Sets an info float vector.
 */
void bcf_set_info_flt_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::vector<float>& values);

/**
 * Gets an info string.
 */
std::string bcf_get_info_str(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::string default_value = ".");

/**
 * Sets an info string.
 */
void bcf_set_info_str(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::string default_value = ".");

/**
 * Gets an info string vector.
 */
std::vector<std::string> bcf_get_info_str_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::string default_value = ".");

/**
 * Sets an info string vector.
 */
void bcf_set_info_str_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::vector<std::string> values);

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
 * Get qual
 */
#define bcf_get_qual(v) ((v)->qual)

/**
 * Get n_flt of a bcf record
 */
#define bcf_get_n_filter(v) ((v)->d.n_flt)

/**
 * Get ith format name
 */
#define bcf_get_format(h, v, i) (h)->id[BCF_DT_ID][(v->d.fmt)[i].id].key
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

/**
 * Creates a dummy header with hs37d5 contigs for testing purposes.
 */
bcf_hdr_t* bcf_create_dummy_hdr();

/**
 * Creates a dummy bcf record representing the variant for testing purposes.
 *
 * @variant - 1:123:ACT:AC/ACCCC
 */
bcf1_t* bcf_create_dummy_record(bcf_hdr_t* h, std::string& variant);

/**********
 *LOG UTILS
 **********/

/**
 * Prints a message to STDERR and abort.
 */
void error(const char * msg, ...);

/**
 * Gives a notice to STDERR.
 */
void notice(const char * msg, ...);

#endif
