/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#include "decompose2.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string ref_fasta_file;
    bool smart;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    
    kstring_t s;
    kstring_t new_alleles;
    kstring_t old_alleles;

    /////////
    //stats//
    /////////
    uint32_t no_variants;
    uint32_t new_no_variants;
    uint32_t no_biallelic;
    uint32_t no_multiallelic;
    uint32_t no_additional_biallelic;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "decomposes multiallelic variants into biallelic in a VCF file.";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::SwitchArg arg_smart("s", "s", "smart decomposition [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            smart = arg_smart.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);

        odw = new BCFOrderedWriter(output_vcf_file);
        odw->link_hdr(odr->hdr);
        //bcf_hdr_append_info_with_backup_naming(odw->hdr, "OLD_MULTIALLELIC", "1", "String", "Original chr:pos:ref:alt encoding", false);
        bcf_hdr_append(odw->hdr, "##INFO=<ID=,Number=1,Type=String,Description=\"Original chr:pos:ref:alt encoding\">\n");
        odw->write_hdr();

        s = {0,0,0};
        old_alleles = {0,0,0};
        new_alleles = {0,0,0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;
        no_biallelic = 0;
        no_multiallelic = 0;

        no_additional_biallelic = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
    }

    /**
     * n choose r.
     */
    uint32_t choose(uint32_t n, uint32_t r)
    {
        if (r>n)
        {
            return 0;
        }
        else if (r==n)
        {
            return 1;
        }
        else if (r==0)
        {
            return 1;
        }
        else
        {
            if (r>(n>>1))
            {
                r = n-r;
            }

            uint32_t num = n;
            uint32_t denum = 1;

            for (uint32_t i=1; i<r; ++i)
            {
                num *= n-i;
                denum *= i+1;
            }

            return num/denum;
        }
    }

    void decompose2()
    {
        bcf1_t* v = bcf_init();
        Variant variant;

        while (odr->read(v))
        {
//            std::cerr << "=============================\n";
//            bcf_print(odr->hdr, v);
//            std::cerr << "=============================\n";

            int32_t n_allele = bcf_get_n_allele(v);

            if (n_allele > 2)
            {
                if (!smart)
                {
                    ++no_multiallelic;
                    no_additional_biallelic += n_allele-2;

                    old_alleles.l = 0;
                    bcf_variant2string(odw->hdr, v, &old_alleles);

                    int32_t rid = bcf_get_rid(v);
                    int32_t pos1 = bcf_get_pos1(v);
                    char** allele = bcf_get_allele(v);

                    int32_t *gt = NULL;
                    int32_t *pl = NULL;
                    float *gl = NULL;
                    int32_t *dp = NULL;
                    size_t no_samples = bcf_hdr_nsamples(odr->hdr);
                    bool has_GT = false;
                    bool has_PL = false;
                    bool has_GL = false;
                    bool has_DP = false;
                    int32_t ploidy = 0;
                    int32_t n_genotype;
                    int32_t n_genotype2;

                    int32_t *gts = NULL;
                    int32_t *pls = NULL;
                    float *gls = NULL;
                    int32_t *dps = NULL;

                    if (no_samples)
                    {
                        bcf_unpack(v, BCF_UN_FMT);

                        int32_t n = 0;
                        int32_t ret = bcf_get_genotypes(odr->hdr, v, &gt, &n);
                        if (ret>0) has_GT = true;
                        ploidy = n/bcf_hdr_nsamples(odr->hdr);
                        n_genotype = bcf_ap2g(n_allele, ploidy);
                        n_genotype2 = bcf_ap2g(2, ploidy);
                        if (ret>0) gts = (int32_t*) malloc(no_samples*n_genotype2*sizeof(int32_t));

                        n=0;
                        ret = bcf_get_format_int32(odr->hdr, v, "PL", &pl, &n);
                        if (ret>0)
                        {
                            has_PL = true;
                            pls = (int32_t*) malloc(no_samples*n_genotype2*sizeof(int32_t));
                        }

                        n=0;
                        ret = bcf_get_format_float(odr->hdr, v, "GL", &gl, &n);
                        if (ret>0)
                        {
                            has_GL = true;
                            gls = (float*) malloc(no_samples*n_genotype2*sizeof(float));
                        }

                        n=0;
                        ret = bcf_get_format_int32(odr->hdr, v, "DP", &dp, &n);
                        if (ret>0)
                        {
                            has_DP = true;
                            dps = (int32_t*) malloc(no_samples*sizeof(int32_t));
                        }
                    }

//                    std::cerr << "*****************************\n";
//                    bcf_print(odw->hdr, v);
//                    std::cerr << "*****************************\n";

                    for (size_t i=1; i<n_allele; ++i)
                    {
                        //bcf_clear(v);
                        bcf1_t* nv = bcf_dup(v);
                        new_alleles.l=0;
                        kputs(allele[0], &new_alleles);
                        kputc(',', &new_alleles);
                        kputs(allele[i], &new_alleles);

                        bcf_update_info_string(odw->hdr, nv, "OLD_MULTIALLELIC", old_alleles.s);
                        bcf_update_alleles_str(odw->hdr, nv, new_alleles.s);

                        if (no_samples)
                        {
                            //get array genotypes
                            for (size_t j=0; j<no_samples; ++j)
                            {
                                if (has_GT)
                                {
                                    for (size_t k=0; k<ploidy; ++k)
                                    {
                                        int32_t _a = gt[j*ploidy+k];
                                        if (_a<0)
                                        {
                                            gts[j*ploidy+k] = _a;
                                        }
                                        else
                                        {
                                            int32_t a = bcf_gt_allele(_a);

                                            if (a)
                                            {
                                                a = a==i ? 1 : -1;
                                            }

                                            gts[j*ploidy+k] =  ((a+1)<<1) | bcf_gt_is_phased(_a);
                                        }
                                    }
                                }

                                if (has_PL)
                                {
                                    pls[j*n_genotype2] = pl[j*n_genotype];
                                    uint32_t index = 0;
                                    for (uint32_t k = 1; k<n_genotype2; ++k)
                                    {
                                        index += choose(ploidy-(k-1)+i-1,i-1);
                                        pls[j*n_genotype2+k] = pl[j*n_genotype+index];
                                    }
                                }

                                if (has_GL)
                                {
                                    gls[j*n_genotype2] = gl[j*n_genotype];
                                    uint32_t index = 0;
                                    for (uint32_t k = 1; k<n_genotype2; ++k)
                                    {
                                        index += choose(ploidy-(k-1)+i-1,i-1);
                                        gls[j*n_genotype2+k] = gl[j*n_genotype+index];
                                    }
                                }

                                if (has_DP)
                                {
                                    dps[j] = dp[j];
                                }
                            }

                            //remove other format values except for GT, PL, GL and DP
                            if (i==1)
                            {
                                bcf_fmt_t *fmt = v->d.fmt;
                                for (size_t j = 0; j < v->n_fmt; ++j)
                                {
                                    const char* tag = odw->hdr->id[BCF_DT_ID][fmt[j].id].key;

                                    if (strcmp(tag,"GT")&&strcmp(tag,"PL")&&strcmp(tag,"GL")&&strcmp(tag,"DP"))
                                    {
                                        bcf_update_format_int32(odw->hdr, nv, tag, 0, 0);
                                    }
                                }
                            }

                            if (has_GT) bcf_update_genotypes(odw->hdr, nv, gts, no_samples*ploidy);
                            if (has_PL) bcf_update_format_int32(odw->hdr, nv, "PL", pls, no_samples*n_genotype2);
                            if (has_GL) bcf_update_format_float(odw->hdr, nv, "GL", gls, no_samples*n_genotype2);
                            if (has_DP) bcf_update_format_int32(odw->hdr, nv, "DP", dps, no_samples);
                        }

                        //remove other format values except for GT, PL, GL and DP
                        bcf_fmt_t *fmt = v->d.fmt;
                        for (size_t j = 0; j < v->n_fmt; ++j)
                        {
                            const char* tag = odw->hdr->id[BCF_DT_ID][fmt[j].id].key;

                            if (strcmp(tag,"GT")&&strcmp(tag,"PL")&&strcmp(tag,"GL")&&strcmp(tag,"DP"))
                            {
                                bcf_update_format_int32(odw->hdr, nv, tag, 0, 0);
                            }
                        }

                        odw->write(nv);
                        bcf_destroy(nv);
                    }

                    if (has_GT) {free(gt); free(gts);}
                    if (has_PL) {free(pl); free(pls);}
                    if (has_GL) {free(gl); free(gls);}
                    if (has_DP) {free(dp); free(dps);}
                }
                else //smart decomposition
                {
                    ++no_multiallelic;
                    no_additional_biallelic += n_allele-2;

                    old_alleles.l = 0;
                    bcf_variant2string(odw->hdr, v, &old_alleles);

                    int32_t rid = bcf_get_rid(v);
                    int32_t pos1 = bcf_get_pos1(v);
                    char** allele = bcf_get_allele(v);
                    uint32_t no_samples = bcf_hdr_nsamples(odr->hdr);

                    char** alleles = (char**) malloc(n_allele*sizeof(char*));
                    for (size_t i=0; i<n_allele; ++i)
                    {
                        alleles[i] = strdup(allele[i]);
                    }

//                    std::cerr << "=============================\n";
//                    bcf_print(odr->hdr, v);
//                    std::cerr << "=============================\n";
//
                    for (size_t i=1; i<n_allele; ++i)
                    {
                        bcf1_t *nv = bcf_dup(v);
                        bcf_set_rid(nv, rid);
                        bcf_set_pos1(nv, pos1);
                        new_alleles.l=0;
                        kputs(allele[0], &new_alleles);
                        kputc(',', &new_alleles);
                        kputs(allele[i], &new_alleles);

                        //////////////////////
                        //split up INFO fields
                        //////////////////////
                        if (nv->n_info)
                        {
                            bcf_unpack(v, BCF_UN_INFO);

                            bcf_hdr_t* hdr = odr->hdr;
                            bcf_info_t *info = v->d.info;

                            for (uint32_t j = 0; j < v->n_info; ++j)
                            {
                                int32_t key = info[j].key;

                                if (key<0)
                                {
                                    fprintf(stderr, "[E::%s] invalid BCF, the INFO key key=%d not present in the header.\n", __func__, key);
                                    abort();
                                }

                                const char* tag = hdr->id[BCF_DT_ID][key].key;
                                int32_t var_len = bcf_hdr_id2length(hdr,BCF_HL_INFO,key);
                                int32_t type = info[j].type;

                                if (var_len==BCF_VL_G)
                                {
                                    if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                                    {
                                        int32_t n = 0;
                                        int32_t* g = 0;
                                        int32_t ret = bcf_get_info_int32(odr->hdr, v, tag, &g, &n);

                                        if (ret>0)
                                        {
                                            int32_t n_genotype = n;
                                            int32_t ploidy = bcf_ag2p(n_allele, n_genotype);
                                            int32_t n_genotype2 = bcf_ap2g(2, ploidy);
                                            int32_t* gs = (int32_t*) malloc(n_genotype2*sizeof(int32_t));

                                            gs[0] = g[0];
                                            uint32_t index = 0;
                                            for (uint32_t k = 1; k<n_genotype2; ++k)
                                            {
                                                index += choose(ploidy-(k-1)+i-1,i-1);
                                                gs[k] = g[index];
                                            }

                                            bcf_update_info_int32(odw->hdr, nv, tag, gs, no_samples*n_genotype2);
                                            free(gs);
                                            free(g);
                                        }
                                    }
                                    else if (type==BCF_BT_FLOAT)
                                    {
                                        int32_t n = 0;
                                        float* g = 0;
                                        int32_t ret = bcf_get_info_float(odr->hdr, v, tag, &g, &n);

                                        if (ret>0)
                                        {
                                            int32_t n_genotype = n;
                                            int32_t ploidy = bcf_ag2p(n_allele, n_genotype);
                                            int32_t n_genotype2 = bcf_ap2g(2, ploidy);
                                            float* gs = (float*) malloc(n_genotype2*sizeof(float));

                                            gs[0] = g[0];
                                            uint32_t index = 0;
                                            for (uint32_t k = 1; k<n_genotype2; ++k)
                                            {
                                                index += choose(ploidy-(k-1)+i-1,i-1);
                                                gs[k] = g[index];
                                            }

                                            bcf_update_info_float(odw->hdr, nv, tag, gs, n_genotype2);
                                            free(gs);
                                            free(g);
                                        }
                                    }
                                    else if (type==BCF_BT_CHAR)
                                    {
                                        //not necessary
                                    }
                                }
                                else if (var_len == BCF_VL_A)
                                {
                                    if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                                    {
                                        int32_t n = 0;
                                        int32_t* a = 0;
                                        int32_t ret = bcf_get_info_int32(odr->hdr, v, tag, &a, &n);
                                        if (ret>0)
                                        {
                                            int32_t as = a[i-1];
                                            bcf_update_info_int32(odw->hdr, nv, tag, &as, 1);
                                        }
                                        free(a);
                                    }
                                    else if (type==BCF_BT_FLOAT)
                                    {
                                        int32_t n = 0;
                                        float* a = 0;
                                        int32_t ret = bcf_get_info_float(odr->hdr, v, tag, &a, &n);
                                        if (ret>0)
                                        {
                                            float as = a[i-1];
                                            bcf_update_info_float(odw->hdr, nv, tag, &as, 1);
                                        }
                                        free(a);
                                    }
                                    else if (type==BCF_BT_CHAR)
                                    {
                                        int32_t n = 0;
                                        char* a = 0;
                                        int32_t ret = bcf_get_info_string(odr->hdr, v, tag, &a, &n);

                                        if (ret>0)
                                        {
                                            char* b = a;
                                            char* c = a;
                                            int32_t e = 0;
                                            while (*b)
                                            {
                                                if (*b==',')
                                                {
                                                    *b = 0;
                                                    ++e;

                                                    if (e==i)
                                                    {
                                                        bcf_update_info_string(odw->hdr, nv, tag, c);
                                                        break;
                                                    }

                                                    c = b+1;
                                                }
                                                else if (*(b+1)==0)
                                                {
                                                    bcf_update_info_string(odw->hdr, nv, tag, c);
                                                    break;
                                                }

                                                ++b;
                                            }

                                            if (e>i)
                                            {
                                                bcf_update_info_string(odw->hdr, nv, tag, NULL);
                                            }

                                            free(a);
                                        }
                                    }
                                }
                                else if (var_len == BCF_VL_R)
                                {
                                    if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                                    {
                                        int32_t n = 0;
                                        int32_t* a = 0;
                                        int32_t ret = bcf_get_info_int32(odr->hdr, v, tag, &a, &n);
                                        if (ret)
                                        {
                                            int32_t* as = (int32_t*) malloc(sizeof(int32_t)*2);
                                            as[0] = a[0];
                                            as[1] = a[i];
                                            bcf_update_info_int32(odw->hdr, nv, tag, as, 2);
                                            free(as);
                                            free(a);
                                        }
                                    }
                                    else if (type==BCF_BT_FLOAT)
                                    {
                                        int32_t n = 0;
                                        float* a = 0;
                                        int32_t ret = bcf_get_info_float(odr->hdr, v, tag, &a, &n);
                                        if (ret>0)
                                        {
                                            float* as = (float*) malloc(sizeof(float)*2);
                                            as[0] = a[0];
                                            as[1] = a[i];
                                            bcf_update_info_float(odw->hdr, nv, tag, as, 2);
                                            free(as);
                                            free(a);
                                        }
                                    }
                                    else if (type==BCF_BT_CHAR)
                                    {
                                        int32_t n = 0;
                                        char* a = 0;
                                        int32_t ret = bcf_get_info_string(odr->hdr, v, tag, &a, &n);

                                        if (ret>0)
                                        {
                                            char* b = a;
                                            char* c = a;
                                            int32_t e = 0;
                                            while (*b)
                                            {
                                                if (*b==',')
                                                {
                                                    *b = 0;
                                                    ++e;

                                                    if (e==i)
                                                    {
                                                        bcf_update_info_string(odw->hdr, nv, tag, c);
                                                        break;
                                                    }

                                                    c = b+1;
                                                }
                                                //last occurrence
                                                else if (*(b+1)==0)
                                                {
                                                    bcf_update_info_string(odw->hdr, nv, tag, c);
                                                    break;
                                                }

                                                ++b;
                                            }

                                            if (e>i)
                                            {
                                                bcf_update_info_string(odw->hdr, nv, tag, NULL);
                                            }

                                            free(a);
                                        }
                                    }
                                }
                                else if (var_len == BCF_VL_FIXED)
                                {
                                    //leave it there
                                }
                                else if (var_len == BCF_VL_VAR)
                                {
                                    //leave it there
                                }
                            }
                        }

                        ////////////////////////
                        //split up FORMAT fields
                        ////////////////////////
                        if (no_samples)
                        {
                            bcf_unpack(v, BCF_UN_FMT);

                            if (nv->n_fmt)
                            {
                                bcf_hdr_t* hdr = odr->hdr;
                                bcf_fmt_t *fmt = v->d.fmt;

                                for (uint32_t j = 0; j < (int32_t)v->n_fmt; ++j)
                                {
                                    int32_t id = fmt[j].id;

                                    if (id<0)
                                    {
                                        fprintf(stderr, "[E::%s] invalid BCF, the FORMAT tag id=%d not present in the header.\n", __func__, id);
                                        abort();
                                    }

                                    const char* tag = hdr->id[BCF_DT_ID][id].key;
                                    int32_t var_len = bcf_hdr_id2length(hdr,BCF_HL_FMT,id);
                                    int32_t type = fmt[j].type;

                                    if (var_len==BCF_VL_G)
                                    {
                                        if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                                        {
                                            int32_t n = 0;
                                            int32_t* g = 0;
                                            int32_t ret = bcf_get_format_int32(odr->hdr, v, tag, &g, &n);
                                            int32_t n_genotype = n/no_samples;

                                            int32_t ploidy = bcf_ag2p(n_allele, n_genotype);
                                            int32_t n_genotype2 = bcf_ap2g(2, ploidy);
                                            int32_t* gs = (int32_t*) malloc(no_samples*n_genotype2*sizeof(int32_t));

                                            //get array genotypes
                                            for (uint32_t j=0; j<no_samples; ++j)
                                            {
                                                gs[j*n_genotype2] = g[j*n_genotype];
                                                uint32_t index = 0;
                                                for (uint32_t k = 1; k<n_genotype2; ++k)
                                                {
                                                    index += choose(ploidy-(k-1)+i-1,i-1);
                                                    gs[j*n_genotype2+k] = g[j*n_genotype+index];
                                                }
                                            }

                                            bcf_update_format_int32(odw->hdr, nv, tag, gs, no_samples*n_genotype2);
                                            free(gs);
                                            free(g);
                                        }
                                        else if (type==BCF_BT_FLOAT)
                                        {
                                            int32_t n = 0;
                                            float* g = 0;
                                            int32_t ret = bcf_get_format_float(odr->hdr, v, tag, &g, &n);
                                            int32_t n_genotype = n/no_samples;
                                            int32_t ploidy = bcf_ag2p(n_allele, n_genotype);
                                            int32_t n_genotype2 = bcf_ap2g(2, ploidy);
                                            float* gs = (float*) malloc(no_samples*n_genotype2*sizeof(float));

                                            for (uint32_t j=0; j<no_samples; ++j)
                                            {
                                                gs[j*n_genotype2] = g[j*n_genotype];
                                                uint32_t index = 0;
                                                for (uint32_t k = 1; k<n_genotype2; ++k)
                                                {
                                                    index += choose(ploidy-(k-1)+i-1,i-1);
                                                    gs[j*n_genotype2+k] = g[j*n_genotype+index];
                                                }
                                            }

                                            bcf_update_format_float(odw->hdr, nv, tag, gs, no_samples*n_genotype2);
                                            free(gs);
                                            free(g);
                                        }
                                        else if (type==BCF_BT_CHAR)
                                        {
                                            //to be implemented
                                        }
                                    }
                                    else if (var_len == BCF_VL_A)
                                    {
                                        if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                                        {
                                            int32_t n = 0;
                                            int32_t* a = 0;
                                            int32_t ret = bcf_get_format_int32(odr->hdr, v, tag, &a, &n);
                                            int32_t* as = (int32_t*) malloc(no_samples*sizeof(int32_t));

                                            if (ret!=(n_allele-1)*no_samples)
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    as[j] = bcf_int32_missing;
                                                }
                                            }
                                            else
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    if (a[j*(n_allele-1)+(i-1)]==bcf_int32_missing ||
                                                        a[j*(n_allele-1)+(i-1)]==bcf_int32_vector_end)
                                                    {
                                                        as[j] = bcf_int32_missing;
                                                    }
                                                    else
                                                    {
                                                        as[j] = a[j*(n_allele-1)+(i-1)];
                                                    }
                                                }
                                            }
                                            bcf_update_format_int32(odw->hdr, nv, tag, as, no_samples);
                                            free(as);
                                            free(a);
                                        }
                                        else if (type==BCF_BT_FLOAT)
                                        {
                                            int32_t n = 0;
                                            float* a = 0;
                                            int32_t ret = bcf_get_format_float(odr->hdr, v, tag, &a, &n);
                                            float* as = (float*) malloc(no_samples*sizeof(float));

                                            if (ret!=(n_allele-1)*no_samples)
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    bcf_float_set_missing(as[j]);;
                                                }
                                            }
                                            else
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    if (bcf_float_is_missing(a[j*(n_allele-1)+(i-1)]) ||
                                                        bcf_float_is_vector_end(a[j*(n_allele-1)+(i-1)]))
                                                    {
                                                        bcf_float_set_missing(as[j]);
                                                    }
                                                    else
                                                    {
                                                        as[j] = a[j*(n_allele-1)+(i-1)];
                                                    }
                                                }
                                            }

                                            bcf_update_format_float(odw->hdr, nv, tag, as, no_samples);
                                            free(as);
                                            free(a);
                                        }
                                        else if (type==BCF_BT_CHAR)
                                        {
                                            //to be implemented
                                        }
                                    }
                                    else if (var_len == BCF_VL_R)
                                    {
                                        if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32)
                                        {
                                            int32_t n = 0;
                                            int32_t* a = 0;
                                            int32_t ret = bcf_get_format_int32(odr->hdr, v, tag, &a, &n);
                                            int32_t* as = (int32_t*) malloc(no_samples*sizeof(int32_t)*2);

                                            if (ret!=n_allele*no_samples)
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    as[j*2] = bcf_int32_missing;
                                                    as[j*2+1] = bcf_int32_vector_end;
                                                }
                                            }
                                            else
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    as[j*2] = a[j*n_allele];
                                                    as[j*2+1] = a[j*n_allele+i];
                                                }
                                            }

                                            bcf_update_format_int32(odw->hdr, nv, tag, as, no_samples*2);
                                            free(as);
                                            free(a);
                                        }
                                        else if (type==BCF_BT_FLOAT)
                                        {
                                            int32_t n = 0;
                                            float* a = 0;
                                            int32_t ret = bcf_get_format_float(odr->hdr, v, tag, &a, &n);
                                            float* as = (float*) malloc(no_samples*sizeof(float));

                                            if (ret!=n_allele*no_samples)
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    bcf_float_set_missing(as[j*2]);
                                                    bcf_float_set_vector_end(as[j*2+1]);
                                                }
                                            }
                                            else
                                            {
                                                for (uint32_t j=0; j<no_samples; ++j)
                                                {
                                                    as[j*2] = a[j*n_allele];
                                                    as[j*2+1] = a[j*n_allele+i];
                                                }
                                            }

                                            bcf_update_format_float(odw->hdr, nv, tag, as, no_samples*2);
                                            free(as);
                                            free(a);
                                        }
                                        else if (type==BCF_BT_CHAR)
                                        {
                                            //to be implemented
                                        }
                                    }
                                    else if (var_len == BCF_VL_FIXED)
                                    {
                                        if (strcmp(tag,"GT")==0)
                                        {
                                            int32_t n = 0;
                                            int32_t* gt = 0;
                                            int32_t ret = bcf_get_genotypes(odr->hdr, v, &gt, &n);

                                            int32_t ploidy = n/bcf_hdr_nsamples(odr->hdr);
                                            int32_t n_genotype = bcf_ap2g(n_allele, ploidy);
                                            int32_t n_genotype2 = bcf_ap2g(2, ploidy);

                                            int32_t* gts = (int32_t*) malloc(no_samples*n_genotype2*sizeof(int32_t));

                                            for (uint32_t j=0; j<no_samples; ++j)
                                            {
                                                for (size_t k=0; k<ploidy; ++k)
                                                {
                                                    int32_t a = gt[j*ploidy+k];
                                                    if (a<0)
                                                    {
                                                        gts[j*ploidy+k] = a;
                                                    }
                                                    else
                                                    {
                                                        int32_t na = bcf_gt_allele(a);

                                                        if (na)
                                                        {
                                                            na = na==i ? 1 : -1;
                                                        }

                                                        gts[j*ploidy+k] =  ((na+1)<<1) | bcf_gt_is_phased(a);
                                                    }
                                                }
                                            }

                                            bcf_update_genotypes(odw->hdr, nv, gts, no_samples*ploidy);
                                            free(gts);
                                            free(gt);
                                        }
                                        else
                                        {
                                            //leave it there...
                                        }
                                    }
                                    else if (var_len == BCF_VL_VAR)
                                    {
                                        //leave it there
                                    }
                                }
                            }
                        }

                        bcf_update_info_string(odw->hdr, nv, "OLD_MULTIALLELIC", old_alleles.s);
                        bcf_update_alleles_str(odw->hdr, nv, new_alleles.s);

                        odw->write(nv);
                        bcf_destroy(nv);
                    }
                }
            }
            else
            {
                odw->write(v);
                ++no_biallelic;
            }

            ++no_variants;
        }

        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "decompose v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [s] smart decomposition   " << (smart ? "true" : "false") << " (experimental)\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        print_int_op("         [i] intervals             ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: no. variants                 : " << no_variants << "\n";
        std::clog << "       no. biallelic variants       : " << no_biallelic << "\n";
        std::clog << "       no. multiallelic variants    : " << no_multiallelic << "\n";
        std::clog << "\n";
        std::clog << "       no. additional biallelics    : " << no_additional_biallelic << "\n";
        std::clog << "       total no. of biallelics      : " << no_additional_biallelic + no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void decompose2(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.decompose2();
    igor.print_stats();
};
