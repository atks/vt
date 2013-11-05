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
//
//#ifndef FT_VCF 
//#define FT_VCF 0 
//#endif
//
//#ifndef FT_VCF_GZ 
//#define FT_VCF 1 
//#endif
//
//#ifndef FT_BCF 
//#define FT_BCF 2 
//#endif

#ifndef HTS_UTILS_H
#define HTS_UTILS_H

#include <iostream>
#include <sstream>
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
#include "htslib/vcfutils.h"
#include "utils.h"

/*******
 *COMMON
 *******/

/**********
 *BAM UTILS
 **********/

//used when a base on a read is not aligned to the human genome reference
//in particular for soft clips and insertions
#define BAM_READ_INDEX_NA -1

/**
 *Gets the start position of the first mapped base in the sequence.
 */
#define bam_get_pos1(b) ((b)->core.pos)

/**
 *Gets the end position of the last mapped base in the sequence.
 */
int32_t bam_get_end_pos1(bam1_t *srec);

/**
 *Gets the read sequence from a bam record
 */
void bam_get_seq_string(bam1_t *srec, kstring_t *seq);

/**
 *Gets the base qualities from a bam record, when N is observed, a placeholder value of 0 is entered
 */
void bam_get_qual_string(bam1_t *srec, kstring_t *qual, char* seq);

/**
 *Gets the cigar sequence from a bam record
 */
#define bam_get_n_cigar_op(b) ((b)->core.n_cigar)

/**
 *Gets the cigar sequence from a bam record
 */
void bam_get_cigar_string(bam1_t *srec, kstring_t *str);

/**
 * Is this sequence the first read?
 */
#define bam_is_fread1(b) (((b)->core.flag&BAM_FREAD1) != 0)

/**
 * Is this sequence the second read?
 */
#define bam_is_fread2(b) (((b)->core.flag&BAM_FREAD2) != 0)

/**
 *Gets the base in the read that is mapped to a genomic position.
 *Returns -1 if it does not exists.
 */
void bam_get_base_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos);

/**
 *Gets the base in the read that is mapped to a genomic position.
 *Returns true if successful, false if the location is a deletion or not on the read.
 */
void bam_get_base_and_qual_and_read_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos, kstring_t* readseq, kstring_t* readqual);

/**
 *Gets flag
 */
#define bam_get_flag(b) ((b)->core.flag)

/**
 *Get map quality
 */
#define bam_get_mapq(b) ((b)->core.qual)

/**
 *Get tid - e.g. chromosome id
 */
#define bam_get_tid(b) ((b)->core.tid)

/**
 *Get read length
 */
#define bam_get_l_qseq(b) ((b)->core.l_qseq)

/**
 *Adds the contigs required for build version hs37d5
 */
void bcf_add_hs37d5_contig_headers(bcf_hdr_t *hdr);

/**
 * Gets a string representation of a variant.
 */
void bcf_get_variant(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);

/**
 * bcf_get_fmt_ptr() - returns pointer to a FORMAT field
 * @header: for access to BCF_DT_ID dictionary
 * @line:   VCF line obtained from vcf_parse1
 * @fmt:    one of GT,PL,...
 *
 * Returns bcf_fmt_t* if the call succeeded, or returns NULL when the field
 * is not available.
 */
bcf_fmt_t *bcf_get_fmt_ptr(const bcf_hdr_t *header, bcf1_t *line, char *tag);

/**********
 *BCF UTILS
 **********/  
/**
 *Get number of chromosomes
 */
#define bcf_get_n_chrom(h) ((h)->n[BCF_DT_CTG])

/**
 * Subsets a record by samples.
 * @imap - indices the subsetted samples
 */
int bcf_subset_samples(const bcf_hdr_t *h, bcf1_t *v, std::vector<int32_t>& imap);

/**
 *Get chromosome name
 */
#define bcf_get_chrom(h, b) ((h)->id[BCF_DT_CTG][(b)->rid].key)

/**
 *Set chromosome name
 */
void bcf_set_chrom(bcf_hdr_t *h, bcf1_t *v, char* chrom);

/**
 *Get RID
 */
#define bcf_get_rid(b) ((b)->rid)

/**
 *Check if variant is passed
 */
bool bcf_is_passed(bcf_hdr_t *h, bcf1_t *v);

/**
 *Get RID
 */
//void bcf_set_rid(b) ((b)->rid)

/**
 *Get chromosome names defined in vcf.h
 */
//bcf_seqnames(bcf_hdr_t* h, int* nseqs)

/**
 *Get 0-based position
 */
#define bcf_get_pos0(b) ((b)->pos)

/**
 *Set 0-based position
 */
#define bcf_set_pos0(v, p) ((v)->pos = (p))

/**
 *Get 1-based position
 */
#define bcf_get_pos1(b) ((b)->pos+1)

/**
 *Set 1-based position
 */
#define bcf_set_pos1(v, p) ((v)->pos = (p)-1);

/**
 *Format variant, returns in chr:pos:ref:alts(,) form
 */
void bcf_format_variant(bcf_hdr_t *h, bcf1_t *v, kstring_t *s);

/**
 *Set allele
 */
void bcf_set_allele(bcf1_t *v, std::vector<std::string> alleles);

/**
 *Get allele
 */
#define bcf_get_allele(b) ((b)->d.allele)

/**
 *Get n_alleles of a bcf record
 */
#define bcf_get_n_allele(b) ((b)->n_allele)

/**
 *Get reference base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_ref(b) ((b)->d.allele[0][0])

/**
 *Get alternative base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_alt(b) ((b)->d.allele[1][0])

/**
 *Get reference allele for a Indel, assumes the record is a biallelic Indel
 */
#define bcf_get_indel_ref(b) ((b)->d.allele[0])

/**
 *Get alternative allele for a Indel, assumes the record is a biallelic Indel
 */
#define bcf_get_indel_alt(b) ((b)->d.allele[1])

/**
 *Get reference allele
 */
#define bcf_get_ref(b) ((b)->d.allele[0])

/**
 *Get alternative allele
 */
#define bcf_get_alt(b, i) ((b)->d.allele[i])

/**
 *Get variant type
 */
#define bcf_get_var_type(b) ((b)->d.var_type)

/**
 *Get samples from bcf header
 */
#define bcf_hdr_get_samples(b) ((b)->samples)

/**
 *Get number of samples in bcf header
 */
int32_t bcf_hdr_get_n_sample(bcf_hdr_t *h);

/**
 *Add single sample to BCF header
 */
void bcf_hdr_add_sample(bcf_hdr_t *h, char *s);

/**
 *Gets sequence names and lengths
 */
void bcf_hdr_get_seqs_and_lens(const bcf_hdr_t *h, const char**& seqs, int32_t*& lens, int *n);

const char *hts_parse_reg1(const char *s, int *beg, int *end);

/**
 *synchronize dictionaries
 */
int bcf_hdr_sync(bcf_hdr_t *h);

/**
 * Formats format and genotypes for individuals
 */
int vcf_format1_format_genotypes(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);

/**
 * Reads header of a VCF file and returns the bcf header object.
 * This wraps around vcf_hdr_read from the original htslib to
 * allow for an alternative header file to be read in.
 *
 * this searches for the alternative header saved as <filename>.hdr
 */
bcf_hdr_t *bcf_alt_hdr_read(htsFile *fp);

/**
 * Formats genotypes for individuals
 */
int vcf_format1_genotypes(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);

/**
 * Formats genotypes for individuals
 */
int vcf_format1_shared_to_missing(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);

/**
 * Returns a bcf_fmt_t pointer associated with tag
 */
bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *h, bcf1_t *v, const char *tag);

/**
 * Get the jth integer value for individual i.
 * @i ith individual
 * @j the jth value
 *
 * returns false when the value is missing.
 */
bool bcf_get_fmt_int(int32_t* val, int32_t i, int32_t j, bcf_fmt_t* fmt);

/**
 * Get string value for individual i.
 *
 * returns false when the value is missing.
 */
bool bcf_get_fmt_str(kstring_t* s, int32_t i, bcf_fmt_t* fmt);

/**
 * Get genotype for individual i.
 * @i ith individual
 *
 * returns false when the value is missing.
 */
bool bcf_get_fmt_gt(kstring_t* s, int32_t i, bcf_fmt_t* f);

/**
 * Returns c string of a single values info tag
 */
char* bcf_get_info1(const bcf_hdr_t *h, bcf1_t *v, const char *tag);

#endif
