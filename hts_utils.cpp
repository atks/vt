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
 
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t) typedef khash_t(vdict) vdict_t;
static bcf_idinfo_t bcf_idinfo_def = { { 15, 15, 15 }, { NULL, NULL, NULL}, -1 };

/**
 *Adds the contigs required for build version hs37b5
 */
void bcf_add_hs37d5_contig_headers(bcf_hdr_t *hdr)
{
    bcf_hdr_append(hdr, "##contig=<ID=1,length=249250621,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=2,length=243199373,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=3,length=198022430,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=4,length=191154276,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=5,length=180915260,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=6,length=171115067,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=7,length=159138663,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=8,length=146364022,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=9,length=141213431,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=10,length=135534747,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=11,length=135006516,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=12,length=133851895,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=13,length=115169878,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=14,length=107349540,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=15,length=102531392,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=16,length=90354753,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=17,length=81195210,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=18,length=78077248,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=19,length=59128983,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=20,length=63025520,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=21,length=48129895,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=22,length=51304566,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=X,length=155270560,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=Y,length=59373566,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=MT,length=16569,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000207.1,length=4262,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000226.1,length=15008,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000229.1,length=19913,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000231.1,length=27386,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000210.1,length=27682,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000239.1,length=33824,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000235.1,length=34474,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000201.1,length=36148,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000247.1,length=36422,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000245.1,length=36651,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000197.1,length=37175,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000203.1,length=37498,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000246.1,length=38154,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000249.1,length=38502,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000196.1,length=38914,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000248.1,length=39786,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000244.1,length=39929,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000238.1,length=39939,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000202.1,length=40103,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000234.1,length=40531,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000232.1,length=40652,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000206.1,length=41001,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000240.1,length=41933,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000236.1,length=41934,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000241.1,length=42152,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000243.1,length=43341,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000242.1,length=43523,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000230.1,length=43691,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000237.1,length=45867,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000233.1,length=45941,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000204.1,length=81310,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000198.1,length=90085,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000208.1,length=92689,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000191.1,length=106433,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000227.1,length=128374,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000228.1,length=129120,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000214.1,length=137718,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000221.1,length=155397,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000209.1,length=159169,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000218.1,length=161147,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000220.1,length=161802,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000213.1,length=164239,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000211.1,length=166566,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000199.1,length=169874,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000217.1,length=172149,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000216.1,length=172294,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000215.1,length=172545,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000205.1,length=174588,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000219.1,length=179198,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000224.1,length=179693,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000223.1,length=180455,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000195.1,length=182896,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000212.1,length=186858,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000222.1,length=186861,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000200.1,length=187035,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000193.1,length=189789,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000194.1,length=191469,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000225.1,length=211173,assembly=b37>");
	bcf_hdr_append(hdr, "##contig=<ID=GL000192.1,length=547496,assembly=b37>");
}

/**
 * Gets a string representation of a variant.
 */
void bcf_get_variant(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
	var->l = 0;
	kputs(bcf_get_chrom(h, v), var);
	kputc(':', var);
	kputw(bcf_get_pos1(v), var);
	for (int32_t i=0; i<v->n_allele; ++i)
	{
		kputc(':', var);
		kputs(bcf_get_alt(v, i), var); 
	}
	
	//std::cerr << var->s << "\n";
}

/**
 *Set chromosome name
 */
void bcf_set_chrom(bcf_hdr_t *h, bcf1_t *v, char* chrom)
{
	vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
	khint_t k = kh_get(vdict, d, chrom);
	if (k == kh_end(d)) 
    {
        // Simple error recovery for chromosomes not defined in the header. It will not help when VCF header has
        // been already printed, but will enable tools like vcfcheck to proceed.
        fprintf(stderr, "[W::%s] contig '%s' is not defined in the header\n", __func__, chrom);
        kstring_t tmp = {0,0,0};
        int l;
        ksprintf(&tmp, "##contig=<ID=%s,length=2147483647>", chrom);
        bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
        free(tmp.s);
        if ( bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) ) bcf_hdr_sync((bcf_hdr_t*)h);
        k = kh_get(vdict, d, chrom);
	}
    v->rid = kh_val(d, k).id;		
};


const char *hts_parse_reg1(const char *s, int *beg, int *end)
{
	int i, k, l, name_end;
	*beg = *end = -1;
	name_end = l = strlen(s);
	// determine the sequence name
	for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
	if (i >= 0) name_end = i;
	if (name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for (i = name_end + 1; i < l; ++i) {
			if (s[i] == '-') ++n_hyphen;
			else if (!isdigit(s[i]) && s[i] != ',') break;
		}
		if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
	}
	// parse the interval
	if (name_end < l) {
		char *tmp;
		tmp = (char*)alloca(l - name_end + 1);
		for (i = name_end + 1, k = 0; i < l; ++i)
			if (s[i] != ',') tmp[k++] = s[i];
		tmp[k] = 0;
		if ((*beg = strtol(tmp, &tmp, 10) - 1) < 0) *beg = 0;
		*end = *tmp? strtol(tmp + 1, &tmp, 10) : 1<<29;
		if (*beg > *end) name_end = l;
	}
	if (name_end == l) *beg = 0, *end = 1<<29;
	return s + name_end;
}

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
 *Gets the read sequence from a bam record
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
 *Gets the base qualities from a bam record, when N is observed, a placeholder value of 0(!, 33 adjusted) is entered
 */
void bam_get_qual_string(bam1_t *srec, kstring_t *qual, char* seq)
{
    qual->l=0;
    uint32_t offset = 0;
    uint8_t* s = bam_get_qual(srec);
    for (int32_t i = 0; i < bam_get_l_qseq(srec); ++i) 
    {
        if (seq[i]=='N')
        {
            kputc('!', qual);
            ++offset;
        }
        else
        {
            kputc(s[i-offset] + 33, qual);
        }
    }
};

/**
 *Gets the cigar string from a bam record
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
 *Gets the base in the read that is mapped to a genomic position.
 Extracts the read sequence and aualities too.
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


/**
 *Check if variant is passed
 */
bool bcf_is_passed(bcf_hdr_t *h, bcf1_t *v)
{    
//    std::cerr << v->d.n_flt << ":" << v->d.flt[0] << "\n";
    return (v->d.n_flt==1 && !strcmp(h->id[BCF_DT_ID][v->d.flt[0]].key,"PASS"));    
}

/**
 *Add single sample to BCF header
 */
void bcf_hdr_add_sample(bcf_hdr_t *h, char *s)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    int ret;
    int k = kh_put(vdict, d, s, &ret);
    if (ret) // absent 
    {   
        kh_val(d, k) = bcf_idinfo_def;
        kh_val(d, k).id = kh_size(d) - 1;
        int n = kh_size(d);
        h->samples = (char**) realloc(h->samples,sizeof(char*)*n);
        h->samples[n-1] = s;
        bcf_hdr_sync(h);
    } 
    else 
    {
        fprintf(stderr, "[W::%s] Duplicated sample name '%s'. Skipped.\n", __func__, s);
    }
}


/**
 *Get number of samples in bcf header
 */
int32_t bcf_hdr_get_n_sample(bcf_hdr_t *h)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    return kh_size(d);
}

/**
 *Gets sequence names and lengths
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


/**
 * Get genotype for individual i.
 * @i ith individual
 *
 * returns false when the value is missing.
 */
bool bcf_get_fmt_gt(kstring_t* s, int32_t i, bcf_fmt_t* f)
{    	
	//for an individual array
	// bcf_fmt_array(s, f->n, f->type, f->p + j * f->size);
	
	if (f->n==0) 
	{
		return false;
	}
//	else if (j>=fmt->n)
//	{
//		std::cerr << "[e] out of index for format value " << j << ">" << fmt->n << "\n";
//		return false;
//	}
    else
    {
		int8_t *x = (int8_t*)(f->p + i * f->size);
		
		s->l = 0;
		int32_t l;
        for (l = 0; l < f->n && x[l] != INT8_MIN; ++l) 
        {
            if (l) kputc("/|"[x[l]&1], s);
            if (x[l]>>1) kputw((x[l]>>1) - 1, s);
            else kputc('.', s);
        }
        if (l == 0) kputc('.', s);
		
		return true;
	}
};

/**
 * Returns c string of a single values info tag
 */
bcf_fmt_t* bcf_get_format(const bcf_hdr_t *h, bcf1_t *v, const char *tag)
{
    bcf_unpack(v, BCF_UN_FMT);
    kstring_t s;
    s.l = s.m = 0;
    s.s = 0;

	for (uint32_t i = 0; i < v->n_fmt; ++i)
	{
		bcf_fmt_t *z = &v->d.fmt[i];

        s.l=0;
		kputs(h->id[BCF_DT_ID][z->id].key, &s);
				
		if (!strcmp(tag,s.s))
		{
		    s.l=0;
            free(s.s);			
			return z;
	    }
	}
	
    if (s.s) free(s.s);
    return 0;
};

/**
 * Returns c string of a single values info tag
 */
char* bcf_get_info1(const bcf_hdr_t *h, bcf1_t *v, const char *tag)
{
    bcf_unpack(v, BCF_UN_INFO);
    kstring_t s;
    s.l = s.m = 0;
    s.s = 0;
    if (v->n_info)
    {
		for (uint32_t i = 0; i < v->n_info; ++i)
		{
			bcf_info_t *z = &v->d.info[i];

            s.l=0;
			kputs(h->id[BCF_DT_ID][z->key].key, &s);
			
			//std::cerr << "INSIDE KEY " << s.s << " " << tag << "\n";
			if (!strcmp(tag,s.s))
			{
			    s.l=0;
    			bcf_fmt_array(&s, z->len, z->type, z->vptr);
    			
    			return (char*) strdup(s.s);
		    }
		}
    }

    return (char*) s.s;
};

/**
 * Adds an INFO field.
 */
int32_t bcf_add_info(bcf_hdr_t *h, bcf1_t *v, int32_t type, char *key, uint8_t *value, int32_t len)
{
    if (type == BCF_BT_CHAR)
    {
        bcf_dec_t *d = &v->d;
    
        vdict_t *dict = (vdict_t*)h->dict[BCF_DT_ID];
        khint_t k = kh_get(vdict, dict, key);
		//info tag not in header, add a header record
        if (k == kh_end(dict) || kh_val(dict, k).info[BCF_HL_INFO] == 15)
        {
            fprintf(stderr, "[W::%s] INFO '%s' is not defined in the header, assuming Type=String\n", __func__, key);
            kstring_t tmp = {0,0,0};
            int l;
            ksprintf(&tmp, "##INFO=<ID=%s,Number=1,Type=String,Description=\"Dummy\">", key);
            bcf_hrec_t *hrec = bcf_hdr_parse_line(h,tmp.s,&l);
            free(tmp.s);
            if ( bcf_hdr_add_hrec((bcf_hdr_t*)h, hrec) ) bcf_hdr_sync((bcf_hdr_t*)h);
        }
        else
        {
            //check for duplicate info fields, 
            int32_t key = kh_val(dict, k).id;
            bcf_info_t *z = NULL;
            uint32_t i = 0;
            for (i=0; i<v->n_info; ++i)
            {
                z = &d->info[i];
                if (z->key == key)
                {
                  break;  
                }    
            }
            
            if (i==v->n_info)
            {    
                ++v->n_info;
      		    hts_expand(bcf_info_t, v->n_info, d->m_info, d->info);
        	}
        	z = &d->info[i];
            z->type = type;
            z->key =  kh_val(dict, k).id;
            z->len = len;
        	z->vptr = value;
        	
        }    
    }
    
    return 1;
}

/**
 *Format variant, returns in chr:pos:ref:alts(,) form
 */
void bcf_format_variant(bcf_hdr_t *h, bcf1_t *v, kstring_t *s) 
{
    s->l = 0;
    kputs(bcf_get_chrom(h, v), s);
    kputs(":", s);
    kputuw(bcf_get_pos1(v), s);
    kputs(":", s);
    kputs(bcf_get_alt(v, 0), s);
    kputs(":", s);
    
    for (uint32_t i=1; i<v->n_allele; ++i)
    {
        if (i>1) kputc(',', s);
        kputs(bcf_get_alt(v, i), s);
    }
}
