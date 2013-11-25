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

//#include <iostream>
//#include <sstream>
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
 * Gets the start position of the first mapped base in the sequence.
 */
#define bam_get_chrom(h, b) ((h)->target_name[(b)->core.tid])

/**
 * Gets the start position of the first mapped base in the sequence.
 */
#define bam_get_pos1(b) ((b)->core.pos)

/**
 * Gets the end position of the last mapped base in the sequence.
 */
int32_t bam_get_end_pos1(bam1_t *srec);

/**
 * Gets the read sequence from a bam record
 */
void bam_get_seq_string(bam1_t *srec, kstring_t *seq);

/**
 * Gets the base qualities from a bam record, when N is observed, a placeholder value of 0 is entered
 */
void bam_get_qual_string(bam1_t *srec, kstring_t *qual, const char* seq);

/**
 * Gets the cigar sequence from a bam record
 */
#define bam_get_n_cigar_op(b) ((b)->core.n_cigar)

/**
 * Gets the cigar sequence from a bam record
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
 * Gets the base in the read that is mapped to a genomic position.
 * Returns -1 if it does not exists.
 */
void bam_get_base_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos);

/**
 * Gets the base in the read that is mapped to a genomic position.
 * Returns true if successful, false if the location is a deletion or not on the read.
 */
void bam_get_base_and_qual_and_read_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos, kstring_t* readseq, kstring_t* readqual);

/**
 * Gets flag
 */
#define bam_get_flag(b) ((b)->core.flag)

/**
 * Get map quality
 */
#define bam_get_mapq(b) ((b)->core.qual)

/**
 *Get tid - e.g. chromosome id
 */
#define bam_get_tid(b) ((b)->core.tid)

/**
 * Get read length
 */
#define bam_get_l_qseq(b) ((b)->core.l_qseq)

/**************
 *BCF HDR UTILS
 **************/

/**
 * Get samples from bcf header
 */
#define bcf_hdr_get_samples(b) ((b)->samples)

/**
 * Get number of samples in bcf header
 */
int32_t bcf_hdr_get_n_sample(bcf_hdr_t *h);

/**
 * Gets sequence names and lengths
 */
void bcf_hdr_get_seqs_and_lens(const bcf_hdr_t *h, const char**& seqs, int32_t*& lens, int *n);

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
 * Gets a string representation of a variant.
 */
void bcf_variant2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var);
 
/**
 * Get number of chromosomes
 */
#define bcf_get_n_chrom(h) ((h)->n[BCF_DT_CTG])

/**
 * Get chromosome name
 */
#define bcf_get_chrom(h, b) ((h)->id[BCF_DT_CTG][(b)->rid].key)

/**
 * Set chromosome name
 */
void bcf_set_chrom(bcf_hdr_t *h, bcf1_t *v, const char* chrom);

/**
 * Get RID
 */
#define bcf_get_rid(b) ((b)->rid)

/**
 * Check if variant is passed
 */
bool bcf_is_passed(bcf_hdr_t *h, bcf1_t *v);

/**
 * Get 0-based position
 */
#define bcf_get_pos0(b) ((b)->pos)

/**
 * Set 0-based position
 */
#define bcf_set_pos0(v, p) ((v)->pos = (p))

/**
 * Get 1-based position
 */
#define bcf_get_pos1(b) ((b)->pos+1)

/**
 * Set 1-based position
 */
#define bcf_set_pos1(v, p) ((v)->pos = (p)-1);

/**
 * Set allele
 */
void bcf_set_allele(bcf1_t *v, std::vector<std::string> alleles);

/**
 * Get allele
 */
#define bcf_get_allele(b) ((b)->d.allele)

/**
 * Get n_alleles of a bcf record
 */
#define bcf_get_n_allele(b) ((b)->n_allele)

/**
 * Get reference base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_ref(b) ((b)->d.allele[0][0])

/**
 * Get alternative base for a SNP, assumes the record is a biallelic SNP
 */
#define bcf_get_snp_alt(b) ((b)->d.allele[1][0])

/**
 * Get reference allele
 */
#define bcf_get_ref(b) ((b)->d.allele[0])

/**
 * Get alternative allele
 */
#define bcf_get_alt(b, i) ((b)->d.allele[i])

/**
 * Get variant type
 */
#define bcf_get_var_type(b) ((b)->d.var_type)

#endif
