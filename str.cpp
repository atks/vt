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

#include "str.h"

/**
 * Constructor.
 */
STRMotif::STRMotif(std::string& ref_fasta_file)
{
    vm = new VariantManip(ref_fasta_file.c_str());
    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL) 
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }
    rfhmm = new RFHMM();
    lfhmm = new LFHMM();

    motifs = kh_init(mdict);

    //update factors
    factors = new int32_t*[32];
    for (size_t i=1; i<=32; ++i)
    {
        factors[i] = new int32_t[32];
        int32_t count = 0;
        
        for (size_t j=1; j<=i; ++j)
        {
            if ((i%j)==0)
            {
                factors[i][count++] = j;
            }
        }
    }
};

/**
 * Destructor.
 */
STRMotif::~STRMotif()
{
    delete vm;
    fai_destroy(fai);
    delete rfhmm;
    delete lfhmm;

    for (size_t i=1; i<=32; ++i)
    {
        free(factors[i]);
    }
    free(factors);
}


/**
 * Annotates STR characteristics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
void STRMotif::annotate(bcf_hdr_t* h, bcf1_t* v)
{
    int32_t no_candidate_motifs;

    const char* chrom = bcf_get_chrom(h, v);
    char* ref = bcf_get_alt(v, 0);
    char* alt = bcf_get_alt(v, 1);
    int32_t start1 = bcf_get_pos1(v);

    if (bcf_get_n_allele(v)==2)
    {
        char** candidate_motifs = suggest_motifs(bcf_get_allele(v), bcf_get_n_allele(v), no_candidate_motifs);

        for (size_t i=0; i<no_candidate_motifs; ++i)
        {
            //should not overwrite model parameters here.

            //        std::cerr << "lflank           : " << lflank << "\n";
//        std::cerr << "RU               : " << ru << "\n";
//        std::cerr << "ref_genome       : " << ref_genome << "\n";
//        std::cerr << "CANDIDATE MOTIFS : ";
//        for (size_t i=0; i<no_candidate_motifs; ++i)
//        {
//            std::cerr << (i?",":"") << candidate_motifs[i];
//        }


          search_flanks(chrom, start1, candidate_motifs[i]);

            //check if fit is good, if good, exit loop;
            //at this point, do not store candidates.
            //idea is to choose only really good fits.
            //
            //as to what exactly a good fit.
            //motif discordance is perfect
            //motif concordance is perfect?
            //loosen definition for long sets.
        }
    }
}

/**
 * Suggests a set of repeat motif candidates in a set of alleles.
 */
char** STRMotif::suggest_motifs(char** alleles, int32_t n_allele, int32_t &no_candidate_motifs)
{
    char *motif;

    //grab all candidate alleles
    for (size_t i=1; i<n_allele; ++i)
    {
        char* ref = alleles[0];
        int32_t ref_len = strlen(ref);
        char *alt = alleles[i];
        int32_t alt_len = strlen(alleles[i]);

        int32_t dlen = alt_len-ref_len;

        //extract fragment
        if (dlen>0)
        {
            motif = &alt[alt_len-dlen];
        }
        else
        {
            motif = &ref[ref_len+dlen];
        }

        int32_t len = abs(dlen);

        std::cerr << dlen << " " << len << " "  << alleles[0] << " " << alleles[i] << "\n";

        char* m = get_shortest_repeat_motif(motif, len);

        int32_t ret;
        khiter_t k;
        if (kh_get(mdict, motifs, m)==kh_end(motifs))
        {
            k = kh_put(mdict, motifs, m, &ret);
            kh_value(motifs, k) = 1;
        }
        else
        {
            kh_value(motifs, k) += 1;
        }
    }

    no_candidate_motifs = kh_size(motifs);

    char** candidate_motifs = (char**) malloc(no_candidate_motifs*sizeof(char*));

    khiter_t k;
    int32_t i = 0;
    for (k=kh_begin(motifs); k!=kh_end(motifs); ++k)
    {
        if (kh_exist(motifs, k))
        {
            candidate_motifs[i] = (char*) kh_key(motifs, k);
        }
    }
    kh_clear(mdict, motifs);

    return candidate_motifs;
}

/**
 * Detect candidate flanks given a motif fit.
 * Update model atttributes.
 */
void STRMotif::search_flanks(const char* chrom, int32_t start1, char* motif)
{
    //given chromosome position and motif
    //attempt to

    //lfhmm
    int32_t lflank_len, ref_len;
    char* lflank, *ru, *rflank, *ref;
    int32_t ru_len;
    std::string qual(2048, 'K');

    //make left flank

    lflank = faidx_fetch_seq(fai, chrom, start1-10, start1-1, &lflank_len);
    
    ref = faidx_fetch_seq(fai, chrom, start1-10, start1+100, &ref_len);
    ++ru;

    lfhmm->set_model(lflank, motif);
    lfhmm->set_mismatch_penalty(5);
    lfhmm->align(ref, qual.c_str());
    lfhmm->print_alignment();

    //get genome position of right flank
    //lflank->

    //extract rflank
    //pick flank
    lflank = faidx_fetch_uc_seq(fai, chrom, start1-10, start1-1, &lflank_len);
    
    
    
    
    
    //rfhmm
    rfhmm->set_model(motif, rflank);
    rfhmm->set_mismatch_penalty(5);
    rfhmm->align(ref, qual.c_str());
    rfhmm->print_alignment();


//
//    //if insert, infuse a repeat unit
//    //choose most appropriate motif
//    if (!ru_len)
//    {
//        lfhmm->set_model(lflank, ru);
//        lfhmm->set_mismatch_penalty(5);
//        lfhmm->align(ref_genome, qual.c_str());
//        lfhmm->print_alignment();
//
//        if (strlen(ref)>strlen(alt))
//        {
//            lflank = faidx_fetch_uc_seq(fai, chrom, start1-10, start1-1, &lflank_len);
//            ////bcf_get_info_string(odr->hdr, v, "RU", &ru, &ru_len);
//            ref_genome = faidx_fetch_uc_seq(fai, chrom, start1-10, start1+100, &ref_genome_len);
//            ru = ref;
//            ++ru;
//        }
//        else //deletion
//        {
//            lflank = faidx_fetch_uc_seq(fai, chrom, start1-10, start1-1, &lflank_len);
//            ru = alt;
//            kstring_t str = {(size_t)lflank_len, (size_t)lflank_len, lflank};
//            ++ru;
//            kputs(ru, &str);
//            lflank_len = str.m;
//            lflank = str.s;
//
//            ref_genome = faidx_fetch_uc_seq(fai, chrom, start1, start1+100, &ref_genome_len);
//            str.l=0; str.s=0; str.m=0;
//            kputs(lflank, &str);
//            kputs(ru, &str);
//            kputs(ref_genome, &str);
//
//            free(ref_genome);
//            ref_genome = str.s;
//            ref_genome_len = str.l;
//        }
//    }

}

/**
 * Extracts the shortest repeat unit in a sequence.
 */
char* STRMotif::get_shortest_repeat_motif(char* allele, int32_t len)
{
    std::cerr << "get shortest repeatmotif " << allele << " : " << len << "\n";

    size_t i = 0;
    size_t sub_motif_len;
    while ((sub_motif_len=factors[len][i])!=len)
    {
        std::cerr << "sub motif len : " << sub_motif_len << " " <<  i << "\n";

        bool exact = true;

        size_t n_sub_motif = len/sub_motif_len;

        for (size_t i=0; i<sub_motif_len; ++i)
        {
            char b;
            for (size_t j=0; j<n_sub_motif; ++j)
            {
                if (j)
                {
                    if (b != allele[j*sub_motif_len+i])
                    {
                        exact = false;
                        break;
                    }
                }
                else
                {
                    b = allele[j*sub_motif_len+i];
                }
            }

            if (!exact) break;
        }

        if (exact) break;
        ++i;
    }

    char *motif = allele+len-sub_motif_len;

    return motif;
};