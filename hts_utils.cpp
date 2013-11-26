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

#include "hts_utils.h"
 
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t) 
typedef khash_t(vdict) vdict_t;

/**********
 *BAM UTILS
 **********/
 
/**
 * Gets the read sequence from a bam record
 */
void bam_get_seq_string(bam1_t *srec, kstring_t *seq)
{
    seq->l=0;
    uint8_t* s = bam_get_seq(srec);
    for (uint16_t i = 0; i < bam_get_l_qseq(srec); ++i)
    {
        kputc("=ACMGRSVTWYHKDBN"[bam_seqi(s, i)], seq);
    }
};

/**
 * Gets the base qualities from a bam record, when N is observed, a placeholder value of 0(!, 33 adjusted) is entered
 */
void bam_get_qual_string(bam1_t *srec, kstring_t *qual, const char* seq)
{
    qual->l=0;
    uint32_t offset = 0;
    uint8_t* s = bam_get_qual(srec);
    for (int32_t i = 0; i < bam_get_l_qseq(srec); ++i) 
    {
        kputc(s[i-offset] + 33, qual);
    }
};

/**
 * Gets the cigar string from a bam record
 */
void bam_get_cigar_string(bam1_t *srec, kstring_t *str)
{
    str->l=0;
    int32_t n_cigar_op = bam_get_n_cigar_op(srec);
    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(srec);
        for (int32_t i = 0; i < n_cigar_op; ++i)
        {
            kputw(bam_cigar_oplen(cigar[i]), str);
            kputc(bam_cigar_opchr(cigar[i]), str);
        }
    }
}

/**
 * Gets the base in the read that is mapped to a genomic position.
 * Extracts the read sequence and aualities too.
 */
void bam_get_base_and_qual_and_read_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos, kstring_t* readseq, kstring_t* readqual)
{
    bam1_core_t *c = &srec->core;
    int32_t rlen = c->l_qseq;
    uint32_t cpos = c->pos; //reference coordinates of the first mapped base
    rpos = 0; //read coordinates

    kstring_t str;
    str.l = str.m = 0, str.s = 0;
    base = 'N';
    qual = 0;

    if (c->n_cigar)
    {
        uint32_t *cigar = bam_get_cigar(srec);
        for (uint32_t i = 0; i < c->n_cigar; ++i)
        {
            char op = bam_cigar_opchr(cigar[i]);
            str.l = 0;
            kputw(bam_cigar_oplen(cigar[i]), &str);
            char* stop;
            uint32_t len = strtol(str.s, &stop, 10);
            assert(stop);

            if (op=='M')
            {
                if (pos>=cpos && pos<=cpos+len-1)
                {
                    rpos += pos-cpos;
                    break;
                }

                cpos += len;
                rpos += len;
            }
            else if (op=='D')
            {
                if (pos>=cpos && pos<=cpos+len-1)
                {
                    rpos = -1;
                    break;
                }

                cpos += len;
            }
            else if (op=='S' || op=='I')
            {
                rpos += len;
            }
        }

        //std::cout << "bpos " << bpos << "\n";
        if (rpos>=0 && rpos<=rlen)
        {
            //sequence
            bam_get_seq_string(srec, readseq);
            base = readseq->s[rpos];
            
            //qual
            bam_get_qual_string(srec, readqual, readseq->s);
            qual = readqual->s[rpos];
        }
        else
        {
            rpos = BAM_READ_INDEX_NA;
        }
    }
//    std::cout << "b: " << base << "\n";
//    std::cout << "q: " << s[bpos-1] << " " << q << "\n";
//    for (uint32_t i = 0; i < c->l_qseq; ++i) std::cerr << ((char)(s[i] + 33));
};
 
 
/**************
 *BCF HDR UTILS
 **************/

/**
 * Reads header of a VCF file and returns the bcf header object.
 * This wraps around vcf_hdr_read from the original htslib to
 * allow for an alternative header file to be read in.
 *
 * this searches for the alternative header saved as <filename>.hdr
 */
bcf_hdr_t *bcf_alt_hdr_read(htsFile *fp)
{
    bcf_hdr_t *h = NULL;
    
    //check for existence of alternative header
    kstring_t alt_hdr_fn = {0, 0, 0};
    kputs(fp->fn, &alt_hdr_fn);
    kputs(".hdr", &alt_hdr_fn); 
    FILE *file = fopen(alt_hdr_fn.s, "r");
    if (!file)
    {
        h = bcf_hdr_read(fp);
    }
    else
    {
        fprintf(stderr, "[I:%s:%d %s] read alternative header for %s\n", __FILE__, __LINE__, __FUNCTION__, fp->fn);
        fclose(file);
        htsFile *alt_hdr = hts_open(alt_hdr_fn.s, "r");
        h = bcf_hdr_read(alt_hdr);
        hts_close(alt_hdr);
    }
    
    if (alt_hdr_fn.m) free(alt_hdr_fn.s);
    return h;
}

/**
 * Get number of samples in bcf header
 */
int32_t bcf_hdr_get_n_sample(bcf_hdr_t *h)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    return kh_size(d);
}

/**
 * Gets sequence names and lengths
 */ 
void bcf_hdr_get_seqs_and_lens(const bcf_hdr_t *h, const char**& seqs, int32_t*& lens, int *n)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    seqs = (const char**) calloc(m,sizeof(const char*));
    lens = (int32_t*) calloc(m,sizeof(int32_t));
    khint_t k;
    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        assert( tid<m );
        seqs[tid] = kh_key(d,k);
        
        lens[tid] = 0;
        bcf_hrec_t *hrec = kh_val(d, k).hrec[0];
        for (int i=0; i<hrec->nkeys; ++i)
        {
            if (!strcmp(hrec->keys[i],"length"))
            {
                lens[tid] = atoi(hrec->vals[i]);
            }
        }
        assert(lens[tid]);
    }
    // sanity check: there should be no gaps
    for (tid=0; tid<m; tid++)
        assert(seqs[tid]);
    *n = m;
}

/**********
 *BCF UTILS
 **********/

/**
 * Gets a string representation of a variant.
 */
void bcf_variant2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;
    kputs(bcf_get_chrom(h, v), var);
    kputc(':', var);
    kputw(bcf_get_pos1(v), var);
    for (int32_t i=0; i<v->n_allele; ++i)
    {
        kputc(',', var);
        kputs(bcf_get_alt(v, i), var); 
    }
}

/**
 *Set chromosome name
 */
void bcf_set_chrom(bcf_hdr_t *h, bcf1_t *v, const char* chrom)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
    khint_t k = kh_get(vdict, d, chrom);
    if (k == kh_end(d)) 
    {
        fprintf(stderr, "[E:%s:%d %s] contig '%s' is not defined in the header\n", __FILE__, __LINE__, __FUNCTION__, chrom);
        kstring_t contig = {0,0,0};
        ksprintf(&contig, "##contig=<ID=%s,length=2147483647>", chrom);
        bcf_hdr_append(h, contig.s);
        if (contig.m) free(contig.s);
        k = kh_get(vdict, d, chrom);
    }
    v->rid = kh_val(d, k).id;       
};

/**
 *Check if variant is passed
 */
bool bcf_is_passed(bcf_hdr_t *h, bcf1_t *v)
{    
//    std::cerr << v->d.n_flt << ":" << v->d.flt[0] << "\n";
    return (v->d.n_flt==1 && !strcmp(h->id[BCF_DT_ID][v->d.flt[0]].key,"PASS"));    
}

