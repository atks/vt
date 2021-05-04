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
#include "Rmath/Rmath.h"

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
    char* seq = faidx_fetch_seq(fai, c_name, p_beg_i, p_end_i, len);
   
    for (int32_t i=0; i<*len; ++i)
    {
        if (isgraph(seq[i])) seq[i] = toupper(seq[i]);
    }
    
    return seq;
}

/**********
 *HTS UTILS
 **********/

/**
 * Checks file extension for use in writing files.
 */
bool str_ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

/**
 * Checks file extension based on filename.
 */
int32_t hts_filename_type(std::string const& value)
{
    if (str_ends_with(value,".vcf") || str_ends_with(value,".vcf.gz"))
    {
        return vcf;
    }
    else if (str_ends_with(value,".bcf"))
    {
        return bcf;
    }
     else if (str_ends_with(value,".sam"))
    {
        return sam;
    }
    else if (str_ends_with(value,".bam"))
    {
        return bam;
    }
    else if (str_ends_with(value,".cram"))
    {
        return cram;
    }
    else
    {
        return unknown_format;
    }
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
    for (size_t i=0; i<bam_hdr_get_n_targets(sh); ++i)
    {
        s.l = 0;
        ksprintf(&s, "##contig=<ID=%s,length=%d>", bam_hdr_get_target_name(sh)[i], bam_hdr_get_target_len(sh)[i]);
        bcf_hdr_append(vh, s.s);
    }
    if (s.m) free(s.s);
}

/**
 * Gets sample names from bam header.
 */
std::string bam_hdr_get_sample_name(bam_hdr_t* hdr) {
  if ( !hdr )
    error("Failed to read the BAM header");

  const char *p = hdr->text;
  const char *q, *r;
  int n = 0;
  std::string sm;
  while( ( q = strstr(p, "@RG" ) ) != 0 ) {
    p = q + 3;
    r = q = 0;
    if ( ( q = strstr(p, "\tID:" ) ) != 0 ) q += 4;
    if ( ( r = strstr(p, "\tSM:" ) ) != 0 ) r += 4;
    if ( r && q ) {
      char *u, *v;
      for (u = (char*)q; *u && *u != '\t' && *u != '\n'; ++u);
      for (v = (char*)r; *v && *v != '\t' && *v != '\n'; ++v);
      *u = *v = '\0';
      if ( sm.empty() )
    sm = r;
      else if ( sm.compare(r) != 0 )
    error("Multiple sample IDs are included in one BAM file - %s, %s", sm.c_str(), r);
    }
    else break;
    p = q > r ? q : r;
    ++n;
  }
  if ( sm.empty() )
    error("Sample ID information cannot be found");
  return sm;
}


/**********
 *BAM UTILS
 **********/

/**
 * Gets the end position of the last mapped base in the read.
 */
int32_t bam_get_end_pos1(bam1_t *s)
{
    int32_t end_pos1 = bam_get_pos1(s);
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(s);
        for (int32_t i = 0; i < n_cigar_op; ++i)
        {
            int32_t opchr = bam_cigar_opchr(cigar[i]);

            if (opchr=='M' || opchr=='D' || opchr=='N' || opchr=='=' || opchr=='X')
            {
                end_pos1 += bam_cigar_oplen(cigar[i]);
            }
        }
    }

    return end_pos1-1;
}

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
 * Extracts the read sequence and qualities too.
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
            //assert(stop);

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

/**
 * Prints a bam.
 */
void bam_print(bam_hdr_t *h, bam1_t *s)
{
    const char* chrom = bam_get_chrom(h, s);
    uint32_t pos1 = bam_get_pos1(s);
    kstring_t seq = {0,0,0};
    bam_get_seq_string(s, &seq);
    uint32_t len = bam_get_l_qseq(s);
    kstring_t qual = {0,0,0};
    bam_get_qual_string(s, &qual);
    kstring_t cigar_string = {0,0,0};
    bam_get_cigar_string(s, &cigar_string);
    kstring_t cigar_expanded_string = {0,0,0};
    bam_get_cigar_expanded_string(s, &cigar_expanded_string);
    uint16_t flag = bam_get_flag(s);
    uint32_t mapq = bam_get_mapq(s);

    std::cerr << "##################" << "\n";
    std::cerr << "chrom-pos: " << chrom << "-" << pos1 << "\n";
    std::cerr << "read     : " << seq.s << "\n";
    std::cerr << "qual     : " << qual.s << "\n";
    std::cerr << "cigar_str: " << cigar_string.s << "\n";
    std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
    std::cerr << "len      : " << len << "\n";
    std::cerr << "mapq     : " << mapq << "\n";
    std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
    std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";

    if (seq.m) free(seq.s);
    if (qual.m) free(qual.s);
    if (cigar_string.m) free(cigar_string.s);
    if (cigar_expanded_string.m) free(cigar_expanded_string.s);
}

/**************
 *BCF HDR UTILS
 **************/

/**
 * Prints header.
 */
void bcf_hdr_print(bcf_hdr_t *h)
{
    if (h->dirty && bcf_hdr_sync(h)<0) 
    {
        fprintf(stderr, "[%s:%d %s] Cannot update header\n", __FILE__, __LINE__, __FUNCTION__);
        exit(1);
    }
    kstring_t htxt = {0,0,0};
    bcf_hdr_format(h, 0, &htxt);
    std::cerr << htxt.s;
    if (htxt.m) free(htxt.s);
}

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
 * Checks if a particular header type exists
 * @hdr  - header
 * @type - BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT, BCF_HL_CTG
 * @key  - the key name
 */
bool bcf_hdr_exists(bcf_hdr_t *hdr, int32_t type, const char *key)
{
    if (type>=BCF_HL_FLT && type<= BCF_HL_CTG) //FLT, INFO, FMT, CTG
    {
        for (uint32_t i=0; i<hdr->nhrec; ++i)
        {
            if (hdr->hrec[i]->type==type)
            {
                int j = bcf_hrec_find_key(hdr->hrec[i], "ID");
                if (j>=0)
                {
                    char* val = hdr->hrec[i]->vals[j];

                    if (strcmp(key, val)==0)
                    {
                        return true;
                    }
                }
            }
        }
    }

    return false;
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
 * Translates BCF_VL types to string.
 */
std::string bcf_hdr_vl2str(int32_t id)
{
    std::cerr << "vl id: " << id << "\n";

    if (id==BCF_VL_FIXED)
    {
        return "BCF_VL_FIXED";
    }
    else if (id==BCF_VL_VAR)
    {
        return "BCF_VL_VAR";
    }
    else if (id==BCF_VL_A)
    {
        return "BCF_VL_A";
    }
    else if (id==BCF_VL_G)
    {
        return "BCF_VL_G";
    }
    else if (id==BCF_VL_R)
    {
        return "BCF_VL_R";
    }
    else
    {
        return "UNRECOGNIZED BCF_VL_* TYPE";
    }
}

/**
 * Translates BCF_HT types to string.
 */
std::string bcf_hdr_ht2str(int32_t id)
{
    if (id==BCF_HT_FLAG)
    {
        return "BCF_HT_FLAG";
    }
    else if (id==BCF_HT_INT)
    {
        return "BCF_HT_INT";
    }
    else if (id==BCF_HT_REAL)
    {
        return "BCF_HT_REAL";
    }
    else if (id==BCF_HT_STR)
    {
        return "BCF_HT_STR";
    }
    else
    {
        return "UNRECOGNIZED BCF_HT_* TYPE";
    }
}


/**********
 *BCF UTILS
 **********/

/**
 * Copies an info fields from one record to another
 * @hsrc  - source header
 * @vsrc  - source bcf1_t
 * @hdest - destination header
 * @vdest - destination bcf1_t
 * @key   - the key name
 * @type  - BCF_HT_FLAG, BCF_HT_INT, BCF_HT_REAL, BCF_HT_STR
 */
void bcf_copy_info_field(bcf_hdr_t *hsrc, bcf1_t* vsrc, bcf_hdr_t *hdest, bcf1_t* vdest, const char* key, int32_t type)
{
    if (type==BCF_HT_FLAG)
    {
        if (bcf_get_info_flag(hsrc, vsrc, key, NULL, NULL)>0)
        {
            bcf_update_info_flag(hdest, vdest, key, NULL, 1);
        }
    }
    else if (type==BCF_HT_INT)
    {
        int32_t* str = NULL;
        int32_t n = 0;
        if (bcf_get_info_int32(hsrc, vsrc, key, &str, &n)>0)
        {
            bcf_update_info_int32(hdest, vdest, key, str, n);
            free(str);
        }
    }
    else if (type==BCF_HT_REAL)
    {
        float* str = NULL;
        int32_t n = 0;
        if (bcf_get_info_float(hsrc, vsrc, key, &str, &n)>0)
        {
            bcf_update_info_float(hdest, vdest, key, str, n);
            free(str);
        }
    }
    else if (type==BCF_HT_STR)
    {
        char** str = NULL;
        int32_t n = 0;
        if (bcf_get_info_string(hsrc, vsrc, key, &str, &n)>0)
        {
            bcf_update_info_string(hdest, vdest, key, str);
            free(str);
        }
    }
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
 * If the VCF files is BCF, any alternative header is ignored.
 */
bcf_hdr_t *bcf_alt_hdr_read(htsFile *fp)
{
    bcf_hdr_t *h = NULL;

    //check for existence of alternative header
    kstring_t alt_hdr_fn = {0, 0, 0};
    kputs(fp->fn, &alt_hdr_fn);
    kputs(".hdr", &alt_hdr_fn);

    FILE *file = NULL;
    struct stat st;
    if (stat(alt_hdr_fn.s, &st)==0 && st.st_size)
    {
        file = fopen(alt_hdr_fn.s, "r");
    }

    if (fp->format.format ==bcf || !file)
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
 * Help function for adding a header with a backup tag name.
 * If the <tag> is already present, a new tag is attempted
 * in the format <tag>_1 to <tag>_9.  If <tag>_9 failed,
 * the function will not add any new tag and will return
 * an empty string.
 *
 * Returns the tag that was inserted or updated.
 */
std::string bcf_hdr_append_info_with_backup_naming(bcf_hdr_t *h, std::string tag, std::string number, std::string type, std::string description, bool rename)
{
    if (bcf_hdr_id2int(h,  BCF_DT_ID, tag.c_str())==-1)
    {
        std::string meta_hdr = "##INFO=<ID=" + tag +
                                 ",Number=" + number +
                                 ",Type=" + type +
                                 ",Description=\"" + description + "\">";

        bcf_hdr_append(h, meta_hdr.c_str());
    }
    else
    {
        if (rename)
        {
            std::string new_tag = "";
            for (uint32_t i=0; i<=9; ++i)
            {
                char c = 49+i;
                new_tag = tag +  "_" + c;
                if (bcf_hdr_id2int(h,  BCF_DT_ID, new_tag.c_str())==-1)
                {
                    std::string meta_hdr = "##INFO=<ID=" + new_tag +
                                 ",Number=" + number +
                                 ",Type=" + type +
                                 ",Description=\"" + description + "\">";

                    bcf_hdr_append(h, meta_hdr.c_str());

                    break;
                }

                new_tag = "";
            }

            return new_tag;
        }
        else
        {
            return tag;
        }
    }

    return tag;
}

/**
 *  bcf_get_format_*() - same as bcf_get_info*() above
 *
 *  The function bcf_get_format_string() is a higher-level (slower) variant of bcf_get_format_char().
 *  see the description of bcf_update_format_string() and bcf_update_format_char() above.
 *  Unlike other bcf_get_format__*() functions, bcf_get_format_string() allocates two arrays:
 *  a single block of \0-terminated strings collapsed into a single array and an array of pointers
 *  to these strings. Both arrays must be cleaned by the user.
 *
 *  Returns negative value on error or the number of written values on success.
 *
 *  Example:
 *      int ndst = 0; char **dst = NULL;
 *      if ( bcf_get_format_string(hdr, line, "XX", &dst, &ndst) > 0 )
 *          for (i=0; i<bcf_hdr_nsamples(hdr); i++) printf("%s\n", dst[i]);
 *      free(dst[0]); free(dst);
 *
 *  Example:
 *      int ngt, *gt_arr = NULL, ngt_arr = 0;
 *      ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
 *
 *  todo: modify to allow direct reading, instead of copying to char*
 */
int32_t bcf_get_format_string_ro(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char ***dst, int *ndst)
{
    int i,tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id) ) return -1;    // no such FORMAT field in the header
    if ( bcf_hdr_id2type(hdr,BCF_HL_FMT,tag_id)!=BCF_HT_STR ) return -2;     // expected different type

    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==tag_id ) break;
    if ( i==line->n_fmt ) return -3;                               // the tag is not present in this record
    bcf_fmt_t *fmt = &line->d.fmt[i];

    int nsmpl = bcf_hdr_nsamples(hdr);
    if ( !*dst )
    {
        *dst = (char**) malloc(sizeof(char*)*nsmpl);
        if ( !*dst ) return -4;     // could not alloc
        (*dst)[0] = NULL;
    }
    int n = (fmt->n+1)*nsmpl;
    if ( *ndst < n )
    {
        (*dst)[0] = (char*) realloc((*dst)[0], n);
        if ( !(*dst)[0] ) return -4;    // could not alloc
        *ndst = n;
    }
    for (i=0; i<nsmpl; i++)
    {
        uint8_t *src = fmt->p + i*fmt->n;
        uint8_t *tmp = (uint8_t*)(*dst)[0] + i*(fmt->n+1);
        memcpy(tmp,src,fmt->n);
        tmp[fmt->n] = 0;
        (*dst)[i] = (char*) tmp;
    }
    return n;
}

/**********
 *BCF UTILS
 **********/

/**
 * Gets number of genotypes from number of alleles and ploidy.
 */
uint32_t bcf_ap2g(uint32_t no_allele, uint32_t no_ploidy)
{
    if (!no_ploidy || !no_allele)
    {
        return 0;
    }
    else if (no_ploidy==1)
    {
        return no_allele;
    }
    else if (no_ploidy==2)
    {
        return (((no_allele+1)*(no_allele))>>1); ;
    }
    else
    {
        return choose(no_ploidy+no_allele-1, no_allele-1);
    }
}

/**
 * Gets alleles from number of ploidy and genotype index.
 */
void bcf_pg2a(uint32_t no_ploidy, uint32_t genotype_index, std::vector<int32_t>& alleles)
{

}

/**
 * Gets number of ploidy from number of alleles and genotypes.
 *
 * Returns 0 if number of alleles and genotypes are not consistent.
 */
uint32_t bcf_ag2p(uint32_t no_alleles, uint32_t no_genotypes)
{
    //hard coded simple cases
    if (no_alleles==2 && no_genotypes==3)
    {
        return 2;
    }
    else if (no_alleles==3 && no_genotypes==6)
    {
        return 2;
    }
    else if (no_alleles==4 && no_genotypes==10)
    {
        return 2;
    }
    else if (no_alleles == no_genotypes)
    {
        return 1;
    }

    //this works in general for all number of alleles and number of genotypes
    uint32_t no_ploidy = 1;
    while (true)
    {
        uint32_t k = choose(no_ploidy+no_alleles-1, no_alleles-1);

        if (k==no_genotypes)
        {
            return no_ploidy;
        }
        //number of alleles and genotypes are not consistent for any ploidy
        else if (k>no_genotypes)
        {
            return 0;
        }

        ++no_ploidy;
    }
}

/**
 * Gets genotype from genotype index and ploidy.
 *
 * The genotype index is computed by a summation of a series which is
 * monotonically decreasing.  This allows you to compute the inverse function
 * from index to the ordered genotypes by using a "water rapids algorithm" with
 * decreasing height of each mini water fall.
 */
std::vector<int32_t> bcf_ip2g(int32_t genotype_index, uint32_t no_ploidy)
{
    std::vector<int32_t> genotype(no_ploidy, 0);

    int32_t pth = no_ploidy;
    int32_t max_allele_index = genotype_index;
    int32_t leftover_genotype_index = genotype_index;

    while (pth>0)
    {
        for (int32_t allele_index=0; allele_index <= max_allele_index; ++allele_index)
        {
            int32_t i = choose(pth+allele_index-1, pth);

            if (i>=leftover_genotype_index || allele_index==max_allele_index)
            {
                if (i>leftover_genotype_index) --allele_index;
                leftover_genotype_index -= choose(pth+allele_index-1, pth);
                --pth;
                max_allele_index = allele_index;
                genotype[pth] = allele_index;
                break;
            }
        }
    }

    return genotype;
}

/**
 * Gets index of a genotype of n ploidy.
 */
uint32_t bcf_g2i(int32_t* g, uint32_t n)
{
    if (n==1)
    {
        return g[0];
    }
    if (n==2)
    {
        return g[0] + (((g[1]+1)*(g[1]))>>1);
    }
    else
    {
        uint32_t index = 0;
        for (uint32_t i=0; i<n; ++i)
        {
            index += bcf_ap2g(g[i], i+1);
        }
        return index;
    }
}

/**
 * Gets index of a genotype of n ploidy.
 */
uint32_t bcf_g2i(std::vector<int32_t>& g)
{
    int32_t n = g.size();

    if (n==1)
    {
        return g[0];
    }
    if (n==2)
    {
        return g[0] + (((g[1]+1)*(g[1]))>>1);
    }
    else
    {
        uint32_t index = 0;
        for (uint32_t i=0; i<n; ++i)
        {
            index += bcf_ap2g(g[i], i+1);
        }
        return index;
    }
}

/**
 * Gets a string representation of a variant.
 */
std::string bcf_variant2string(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t var = {0,0,0};
    bcf_variant2string(h, v, &var);
    std::string var2(var.s);
    if (var.m) free(var.s);
    return var2;
}

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
    for (size_t i=0; i<bcf_get_n_allele(v); ++i)
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
 * Returns true if a is before b, false otherwise.
 */
bool bcf_is_in_order(bcf1_t *a, bcf1_t *b)
{
    if (bcf_get_rid(a)==bcf_get_rid(b))
    {
        return bcf_get_pos0(a)<=bcf_get_pos0(b);
    }

    return bcf_get_rid(a)<bcf_get_rid(b);
}

/**
 * Returns a copy v that only has the chr:pos1:ref:alt information.
 */
bcf1_t* bcf_copy_variant(bcf_hdr_t *h, bcf1_t *v)
{
    bcf1_t* nv = bcf_init1();
    bcf_clear(nv);
    bcf_set_n_sample(nv, bcf_get_n_sample(v));

    bcf_set_rid(nv, bcf_get_rid(v));
    bcf_set_pos1(nv, bcf_get_pos1(v));
    kstring_t s = {0,0,0};
    bcf_alleles2string(h, v, &s);
    bcf_update_alleles_str(h, nv, s.s);
    if (s.m) free(s.s);

    return nv;
}

/**
 * Gets a sorted string representation of a variant.
 */
void bcf_variant2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
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
        for (size_t i=1; i<v->n_allele; ++i)
        {
            temp[i] = allele[i];
        }
        std::qsort(temp, bcf_get_n_allele(v), sizeof(char*), cmpstr);
        kputs(bcf_get_alt(v, 0), var);
        for (size_t i=0; i<v->n_allele-1; ++i)
        {
            kputc('/', var);
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
 * Get chromosome name
 */
const char* bcf_get_chrom(bcf_hdr_t *h, bcf1_t *v)
{
    if (v->rid >= h->n[BCF_DT_CTG])
    {
        fprintf(stderr, "[E:%s:%d %s] rid '%d' does not have an associated contig defined in the header.  Try tabix workaround or just add the header.\n", __FILE__, __LINE__, __FUNCTION__, v->rid);
        exit(1);
    }
    else
    {
        return h->id[BCF_DT_CTG][v->rid].key;
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

/**
 * Gets an info flag.
 */
bool bcf_get_info_flg(bcf_hdr_t *h, bcf1_t *v, const char* tag)
{
    return (bcf_get_info_flag(h, v, tag, 0, 0) ? true : false);
}

/**
 * Sets an info flag.
 */
void bcf_set_info_flg(bcf_hdr_t *h, bcf1_t *v, const char* tag, bool value)
{
    bcf_update_info_flag(h, v, tag, "", (value ? 1 : 0));
}

/**
 * Gets an info integer.
 */
int32_t bcf_get_info_int(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t default_value)
{
    int32_t i = 0;
    int32_t* temp_i = NULL;
    int32_t n = 0;
    if (bcf_get_info_int32(h, v, tag, &temp_i, &n)>0)
    {
        i = *temp_i;
        free(temp_i);
    }
    else
    {
        return default_value;
    }

    return i;
}

/**
 * Sets an info integer.
 */
void bcf_set_info_int(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t value)
{
    bcf_update_info_int32(h, v, tag, &value, 1);
}

/**
 * Gets an info integer vector.
 */
std::vector<int32_t> bcf_get_info_int_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t default_size, int32_t default_value)
{
    std::vector<int32_t> i_vec;
    int32_t* temp_i = NULL;
    int32_t n = 0;
    if (bcf_get_info_int32(h, v, tag, &temp_i, &n)>0)
    {
        for (int32_t j=0; j<n; ++j) i_vec.push_back(temp_i[j]);
        free(temp_i);
    }
    else
    {
        i_vec.resize(default_size, default_value);
    }

    return i_vec;
}

/**
 * Sets an info integer vector.
 */
void bcf_set_info_int_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::vector<int32_t>& values)
{
    bcf_update_info_int32(h, v, tag, &values[0], values.size());
}

/**
 * Sets an info integer vector.
 */
void bcf_set_info_int(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::vector<int32_t>& values)
{
    bcf_update_info_int32(h, v, tag, &values[0], values.size());
}

/**
 * Gets an info float.
 */
float bcf_get_info_flt(bcf_hdr_t *h, bcf1_t *v, const char* tag, float default_value)
{
    float f = 0;
    float* temp_f = NULL;
    int32_t n = 0;
    if (bcf_get_info_float(h, v, tag, &temp_f, &n)>0)
    {
        f = *temp_f;
        free(temp_f);
    }
    else
    {
        return default_value;
    }

    return f;
}

/**
 * Sets an info float.
 */
void bcf_set_info_flt(bcf_hdr_t *h, bcf1_t *v, const char* tag, float value)
{
    bcf_update_info_float(h, v, tag, &value, 1);
}

/**
 * Gets an info integer vector.
 */
std::vector<float> bcf_get_info_flt_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, int32_t default_size, float default_value)
{
    std::vector<float> vec;

    return vec;
}

/**
 * Sets an info float vector.
 */
void bcf_set_info_flt_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::vector<float>& values)
{
}

/**
 * Gets an info string.
 */
std::string bcf_get_info_str(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::string default_value)
{
    std::string str = "";
    char* s = NULL;
    int32_t n = 0;
    if (bcf_get_info_string(h, v, tag, &s, &n)>0)
    {
        str.assign(s);
        free(s);
    }
    else
    {
        return default_value;
    }

    return str;
}

/**
 * Sets an info string.
 */
void bcf_set_info_str(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::string default_value)
{
}

/**
 * Gets an info string vector.
 */
std::vector<std::string> bcf_get_info_str_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::string default_value)
{
    std::vector<std::string> vec;
    std::string str = "";
    char* s = NULL;
    int32_t n = 0;
    int32_t ret = 0;
    if ((ret=bcf_get_info_string(h, v, tag, &s, &n))>0)
    {
        str.assign(s);
        free(s);
        vec.push_back(str);
    }
    else
    {
        vec.push_back(default_value);
    }

    return vec;
}
/**
 * Sets an info string vector.
 */
void bcf_set_info_str_vec(bcf_hdr_t *h, bcf1_t *v, const char* tag, std::vector<std::string> values)
{
}

/**
 * Creates a dummy header with hs37d5 contigs for testing purposes.
 */
bcf_hdr_t* bcf_create_dummy_hdr()
{
    bcf_hdr_t* h = bcf_hdr_init("r");
    bcf_hdr_append(h, "##fileformat=VCFv4.1");
    bcf_hdr_append(h, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    bcf_hdr_append(h, "##contig=<ID=1,assembly=b37,length=249250621>");
    bcf_hdr_append(h, "##contig=<ID=2,assembly=b37,length=243199373>");
    bcf_hdr_append(h, "##contig=<ID=3,assembly=b37,length=198022430>");
    bcf_hdr_append(h, "##contig=<ID=4,assembly=b37,length=191154276>");
    bcf_hdr_append(h, "##contig=<ID=5,assembly=b37,length=180915260>");
    bcf_hdr_append(h, "##contig=<ID=6,assembly=b37,length=171115067>");
    bcf_hdr_append(h, "##contig=<ID=7,assembly=b37,length=159138663>");
    bcf_hdr_append(h, "##contig=<ID=8,assembly=b37,length=146364022>");
    bcf_hdr_append(h, "##contig=<ID=9,assembly=b37,length=141213431>");
    bcf_hdr_append(h, "##contig=<ID=10,assembly=b37,length=135534747>");
    bcf_hdr_append(h, "##contig=<ID=11,assembly=b37,length=135006516>");
    bcf_hdr_append(h, "##contig=<ID=12,assembly=b37,length=133851895>");
    bcf_hdr_append(h, "##contig=<ID=13,assembly=b37,length=115169878>");
    bcf_hdr_append(h, "##contig=<ID=14,assembly=b37,length=107349540>");
    bcf_hdr_append(h, "##contig=<ID=15,assembly=b37,length=102531392>");
    bcf_hdr_append(h, "##contig=<ID=16,assembly=b37,length=90354753>");
    bcf_hdr_append(h, "##contig=<ID=17,assembly=b37,length=81195210>");
    bcf_hdr_append(h, "##contig=<ID=18,assembly=b37,length=78077248>");
    bcf_hdr_append(h, "##contig=<ID=19,assembly=b37,length=59128983>");
    bcf_hdr_append(h, "##contig=<ID=20,assembly=b37,length=63025520>");
    bcf_hdr_append(h, "##contig=<ID=21,assembly=b37,length=48129895>");
    bcf_hdr_append(h, "##contig=<ID=22,assembly=b37,length=51304566>");
    bcf_hdr_append(h, "##contig=<ID=X,assembly=b37,length=155270560>");
    bcf_hdr_append(h, "##contig=<ID=Y,assembly=b37,length=59373566>");
    bcf_hdr_append(h, "##contig=<ID=MT,assembly=b37,length=16569>");

    return h;
}

/**
 * Creates a dummy bcf record representing the variant for testing purposes.
 *
 * @variant - 1:123:ACT:AC/ACCCC
 */
bcf1_t* bcf_create_dummy_record(bcf_hdr_t* h, std::string& variant)
{
    bcf1_t* v = bcf_init1();

    std::vector<std::string> var;
    split(var, ":", variant);

    bcf_set_rid(v, bcf_hdr_name2id(h, var[0].c_str()));
    bcf_set_pos1(v, str2int32(var[1]));

    std::vector<std::string> alleles;
    split(alleles, ":", var[3]);

    std::string new_alleles = var[2];

    for (uint32_t i=0; i<alleles.size(); ++i)
    {
        new_alleles.append(1, ',');
        new_alleles.append(alleles[i]);
    }

    bcf_update_alleles_str(h, v, new_alleles.c_str());

    return v;
}

/**
 * Prints a message to STDERR and abort.
 */
void error(const char * msg, ...)
{
    va_list ap;
    va_start(ap, msg);

    fprintf(stderr, "\nFATAL ERROR - \n");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, "\n\n");

    va_end(ap);

    abort();
    //throw pexception;
    //exit(EXIT_FAILURE);
}

/**
 * Prints a message to STDERR.
 */
void notice(const char * msg, ...)
{
    va_list ap;
    va_start(ap, msg);

    time_t current_time;
    char buff[255];
    current_time = time(NULL);

    strftime(buff, 120, "%Y/%m/%d %H:%M:%S", localtime(&current_time));

    fprintf(stderr,"NOTICE [%s] - ", buff);
    vfprintf(stderr, msg, ap);
    fprintf(stderr,"\n");

    va_end(ap);
}
