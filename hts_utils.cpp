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

/********
 *General
 ********/

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

/**********
 *FAI UTILS
 **********/

/**
 * An alternate sequence fetcher for upper case sequence.
 */
char *faidx_fetch_uc_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    int l;
    char c;
    khiter_t iter;
    faidx1_t val;
    char *seq=NULL;

    // Adjust position
    iter = kh_get(s, fai->hash, c_name);
    if(iter == kh_end(fai->hash)) return 0;
    val = kh_value(fai->hash, iter);
    if(p_end_i < p_beg_i) p_beg_i = p_end_i;
    if(p_beg_i < 0) p_beg_i = 0;
    else if(val.len <= p_beg_i) p_beg_i = val.len - 1;
    if(p_end_i < 0) p_end_i = 0;
    else if(val.len <= p_end_i) p_end_i = val.len - 1;

    // Now retrieve the sequence
    int ret = bgzf_useek(fai->bgzf, val.offset + p_beg_i / val.line_blen * val.line_len + p_beg_i % val.line_blen, SEEK_SET);
    if ( ret<0 )
    {
        *len = -1;
        fprintf(stderr, "[fai_fetch_seq] Error: fai_fetch failed. (Seeking in a compressed, .gzi unindexed, file?)\n");
        return NULL;
    }
    l = 0;
    seq = (char*)malloc(p_end_i - p_beg_i + 2);
    while ( (c=bgzf_getc(fai->bgzf))>=0 && l < p_end_i - p_beg_i + 1)
        if (isgraph(c)) seq[l++] = toupper(c);
    seq[l] = '\0';
    *len = l;
    return seq;
}

/**************
 *BAM HDR UTILS
 **************/

/**
 * Copies contigs found in bam header to bcf header.
 */
void bam_hdr_transfer_contigs_to_bcf_hdr(const bam_hdr_t *sh, bcf_hdr_t *vh)
{
    kstring_t s = {0,0,0};
    for (uint32_t i=0; i<bam_hdr_get_n_targets(sh); ++i)
    {
        s.l = 0;
        ksprintf(&s, "##contig=<ID=%s,length=%d>", bam_hdr_get_target_name(sh)[i], bam_hdr_get_target_len(sh)[i]);
        bcf_hdr_append(vh, s.s);
    }
    if (s.m) free(s.s);
}

/**********
 *BAM UTILS
 **********/

/**
 * Gets the read sequence from a bam record
 */
void bam_get_seq_string(bam1_t *s, kstring_t *seq)
{
    seq->l=0;
    uint8_t* sq = bam_get_seq(s);
    for (uint16_t i = 0; i < bam_get_l_qseq(s); ++i)
    {
        kputc("=ACMGRSVTWYHKDBN"[bam_seqi(sq, i)], seq);
    }
};

/**
 * Gets the base qualities from a bam record, when N is observed, a placeholder value of 0(!, 33 adjusted) is entered
 */
void bam_get_qual_string(bam1_t *s, kstring_t *qual)
{
    qual->l=0;
    uint32_t offset = 0;
    uint8_t* q = bam_get_qual(s);
    for (int32_t i = 0; i < bam_get_l_qseq(s); ++i)
    {
        kputc(q[i-offset] + 33, qual);
    }
};

/**
 * Gets the cigar from a BAM record
 */
void bam_get_cigar_string(bam1_t *s, kstring_t *cigar_string)
{
    cigar_string->l=0;
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(s);
        for (int32_t i = 0; i < n_cigar_op; ++i)
        {
            kputw(bam_cigar_oplen(cigar[i]), cigar_string);
            kputc(bam_cigar_opchr(cigar[i]), cigar_string);
        }
    }
}

/**
 * Gets the cigar string from a bam record
 */
void bam_get_cigar_expanded_string(bam1_t *s, kstring_t *cigar_expanded_string)
{
    kstring_t cigar_string = {0,0,0};
    bam_get_cigar_string(s, &cigar_string);

    cigar_expanded_string->l = 0;
    int32_t lastIndex = cigar_string.l;
    int32_t i = 0;
    kstring_t token = {0,0,0};

    if (lastIndex<0)
    {
        return;
    }
    char c;
    bool seenM = false;

    while (i<=lastIndex)
    {
        c = cigar_string.s[i];

        //captures the numeric count
        if (c<'A')
        {
            kputc(c, &token);
        }

        if (c>'A' || i==lastIndex)
        {
            //it is possible for I's to be observed before the first M's in the cigar string
            //in this case, we treat them as 'S'
            if (!seenM)
            {
                if (c=='I')
                {
                    c = 'S';
                }
                else if (c=='M')
                {
                    seenM = true;
                }
            }

            int32_t count = atoi(token.s);
            for (uint32_t j=0; j<count; ++j)
                kputc(c, cigar_expanded_string);
            token.l = 0;;
        }

        ++i;
    }

    if (cigar_string.m) free(cigar_string.s);
    if (token.m) free(token.s);
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
            bam_get_qual_string(srec, readqual);
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
 * Copies contigs found in bcf header to another bcf header.
 */
void bcf_hdr_transfer_contigs(const bcf_hdr_t *hsrc, bcf_hdr_t *hdest)
{
    vdict_t *d = (vdict_t*)hsrc->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    const char **names = (const char**) calloc(m,sizeof(const char*));
    int len[m];
    khint_t k;


    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        len[tid] = bcf_hrec_find_key(kh_val(d, k).hrec[0],"length");
        int j;
        if ( sscanf(kh_val(d, k).hrec[0]->vals[len[tid]],"%d",&j) )
            len[tid] = j;
        names[tid] = kh_key(d,k);
    }

    kstring_t s = {0,0,0};
    for (tid=0; tid<m; tid++)
    {
        s.l = 0;
        ksprintf(&s, "##contig=<ID=%s,length=%d>", names[tid], len[tid]);
        bcf_hdr_append(hdest, s.s);
    }
    if (s.m) free(s.s);
}

/**
 * Extracts sequence length by rid.
 */
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq)
{
    vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    int32_t *len = (int32_t*) malloc(m*sizeof(int32_t));
    khint_t k;

    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        len[tid] = bcf_hrec_find_key(kh_val(d, k).hrec[0],"length");
        int j;
        if ( sscanf(kh_val(d, k).hrec[0]->vals[len[tid]],"%d",&j) )
            len[tid] = j;
    }

    return len;
}

/**
 * Prints a VCF record to STDERR.
 */
void bcf_print(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    vcf_format(h, v, &s);
    std::cerr << s.s;
    if (s.m) free(s.s);
};

/**
 * Prints a VCF record in compact string representation to STDERR.
 */
void bcf_print_liten(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    bcf_variant2string(h, v, &s);
    std::cerr << s.s << "\n";
    if (s.m) free(s.s);
};

/**
 * Prints a VCF record in compact string representation to STDERR.
 */
void bcf_print_lite(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    bcf_variant2string(h, v, &s);
    std::cerr << s.s;
    if (s.m) free(s.s);
};

/**
 * Prints a VCF record in compact string representation to STDERR with alleles sorted.
 */
void bcf_print_lite_sorted(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    bcf_variant2string_sorted(h, v, &s);
    std::cerr << s.s;
    if (s.m) free(s.s);
};

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

        //helps move the pointer to the right place
        bcf_hdr_t *temp_h = bcf_hdr_read(fp);
        bcf_hdr_destroy(temp_h);
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
 * Checks if a info header exists.
 */
bool bcf_hdr_info_exists(bcf_hdr_t *h, const char* key)
{
    int id = bcf_hdr_id2int(h, BCF_DT_ID, key);
    if (bcf_hdr_idinfo_exists(h,BCF_HL_INFO,id))
    {
        return true;
    }

    return false;
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
    kputc(':', var);
    for (int32_t i=0; i<v->n_allele; ++i)
    {
        if (i) kputc('/', var);
        kputs(bcf_get_alt(v, i), var);
    }
}

/**
 * strcmp wrapper for qsort.
 */
int32_t cmpstr(void const *a, void const *b)
{
    char const *aa = *(char const **)a;
    char const *bb = *(char const **)b;

    return strcmp(aa, bb);
}

/**
 * Gets a sorted string representation of a variant.
 */
void bcf_variant2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_print_liten(h,v);

    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;
    kputs(bcf_get_chrom(h, v), var);
    kputc(':', var);
    kputw(bcf_get_pos1(v), var);
    kputc(':', var);

    if (v->n_allele==2)
    {
        kputs(bcf_get_alt(v, 0), var);
        kputc(',', var);
        kputs(bcf_get_alt(v, 1), var);
    }
    else
    {
        char** allele = bcf_get_allele(v);
        char** temp = (char**) malloc((bcf_get_n_allele(v)-1)*sizeof(char*));
        for (int32_t i=1; i<v->n_allele; ++i)
        {
            temp[i] = allele[i];
        }
        std::qsort(temp, bcf_get_n_allele(v), sizeof(char*), cmpstr);
        kputs(bcf_get_alt(v, 0), var);
        for (int32_t i=0; i<v->n_allele-1; ++i)
        {
            kputc(',', var);
            kputs(temp[i], var);
        }
        free(temp);
    }
}

/**
 * Gets a string representation of the alleles of a variant.
 */
void bcf_alleles2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;

    if (v->n_allele==2)
    {
        kputs(bcf_get_alt(v, 0), var);
        kputc(',', var);
        kputs(bcf_get_alt(v, 1), var);
    }
    else
    {
        char** allele = bcf_get_allele(v);
        for (int32_t i=0; i<v->n_allele; ++i)
        {
            if (i) kputc(',', var);
            kputs(allele[i], var);
        }
    }
}

/**
 * Gets a sorted string representation of the alleles of a variant.
 */
void bcf_alleles2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;

    if (v->n_allele==2)
    {
        kputs(bcf_get_alt(v, 0), var);
        kputc(',', var);
        kputs(bcf_get_alt(v, 1), var);
    }
    else
    {
        char** allele = bcf_get_allele(v);
        char** temp = (char**) malloc((bcf_get_n_allele(v)-1)*sizeof(char*));
        for (int32_t i=1; i<v->n_allele; ++i)
        {
            temp[i-1] = allele[i];
        }

        std::qsort(temp, bcf_get_n_allele(v)-1, sizeof(char*), cmpstr);

        kputs(bcf_get_alt(v, 0), var);
        for (int32_t i=0; i<v->n_allele-1; ++i)
        {
            kputc(',', var);
            kputs(temp[i], var);
        }

        free(temp);
    }
}

/**
 *Set chromosome name.
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
 * Set id.
 */
void bcf_set_id(bcf1_t *v, char* id)
{
    if (v->d.id)
    {
        free(v->d.id);
    }

    v->d.id = strdup(id);
};