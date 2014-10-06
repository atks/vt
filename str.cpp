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
    factors = NULL;
    initialize_factors(32);
};

/**
 * Constructor.
 */
void STRMotif::initialize_factors(int32_t max_len)
{
    if (factors)
    {
        for (size_t i=1; i<=max_len; ++i)
        {
            free(factors[i]);
        }
        free(factors);
    }

    this->max_len = max_len;

    factors = new int32_t*[max_len+1];
    for (size_t i=1; i<=max_len; ++i)
    {
        factors[i] = new int32_t[max_len];
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

    if (factors)
    {
        for (size_t i=1; i<=max_len; ++i)
        {
            free(factors[i]);
        }
        free(factors);
    }
}

/**
 * Annotates STR characteristics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
void STRMotif::annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant)
{
    int32_t no_candidate_motifs;

    if (bcf_get_n_allele(v)==2)
    {
        const char* chrom = bcf_get_chrom(h, v);
        std::string ref(bcf_get_alt(v, 0));
        std::string alt(bcf_get_alt(v, 1));
        int32_t pos1 = bcf_get_pos1(v);

        /////////////
        //EXACT STRs
        ////////////
        trim(pos1, ref, alt);
//        std::cerr << "after trimming : " << pos1 << " " << ref << "\n";
//        std::cerr << "               : " << pos1 << " " << alt << "\n";

        int32_t end1 = pos1 + ref.size() - 1;
        right_align(chrom, end1, ref, alt);
//        std::cerr << "after right alignment : " << end1 << " " << ref << "\n";
//        std::cerr << "                      : " << end1 << " " << alt << "\n";

        //do this later as you want to pefrom the alignment shortly.
        int32_t beg1 = end1 - ref.size() + 1;
        left_align(chrom, beg1, ref, alt);
//        std::cerr << "after left alignment : " << beg1 << " " << ref << "\n";
//        std::cerr << "                     : " << beg1 << " " << alt << "\n";

        variant.eregion.beg1 = beg1;
        variant.eregion.end1 = end1;


        std::string emotif = pick_motif(ref, alt);

        std::cerr << "exact motif : " << emotif << "\n";

        variant.emotif = emotif;

        //use exact motif to perform left and right alignment

        ///////////////
        //INEXACT STRs
        //////////////

       bool debug = true;

        if (debug)
        {
            std::cerr << "//////////////////////\n";
            std::cerr << "//FUZZY LEFT ALIGNMENT\n";
            std::cerr << "//////////////////////\n";
        }

        int32_t lflank_len, rflank_len, lref_len, rref_len;
        char* lflank, *ru, *rflank, *lrefseq, *rrefseq;
        int32_t ru_len;
        std::string qual(2048, 'K');

        int32_t flank_len = 10;
        int32_t check_len = 100;
        
        //left flank
        lflank = faidx_fetch_seq(fai, chrom, pos1-flank_len, pos1-1, &lflank_len);
        lrefseq = faidx_fetch_seq(fai, chrom, pos1-flank_len, pos1+check_len, &lref_len);


        int32_t penalty = emotif.size()>6? 5 : emotif.size();

        std::cerr << "PENALTY : " << penalty << "\n";

        lfhmm->set_model(lflank, emotif.c_str());
        lfhmm->set_mismatch_penalty(penalty);
        lfhmm->align(lrefseq, qual.c_str());
        if (debug) lfhmm->print_alignment();

//    std::cerr << "model: " << "(" << lflank_start[MODEL] << "~" << lflank_end[MODEL] << ") "
//                          << "[" << motif_start[MODEL] << "~" << motif_end[MODEL] << "]\n";
//    std::cerr << "read : " << "(" << lflank_start[READ] << "~" << lflank_end[READ] << ") "
//                          << "[" << motif_start[READ] << "~" << motif_end[READ] << "]"
//                          << "[" << rflank_start[READ] << "~" << rflank_end[READ] << "]\n";

//        int32_t lflank_start[2], lflank_end[2], motif_start[2], motif_end[2], rflank_start[2], rflank_end[2];
//        int32_t motif_count, exact_motif_count, motif_m, motif_xid;


        if (debug) std::cerr << "\n#################################\n\n";

        if (debug)
        {
            std::cerr << "///////////////////////\n";
            std::cerr << "//FUZZY RIGHT ALIGNMENT\n";
            std::cerr << "///////////////////////\n";
        }

        
        int32_t motif_end1 = pos1+lfhmm->get_motif_read_epos1()-flank_len;

        //right flank
        rflank = faidx_fetch_seq(fai, chrom, motif_end1, motif_end1+flank_len-1, &rflank_len);
        rrefseq = faidx_fetch_seq(fai, chrom, motif_end1-check_len, motif_end1+flank_len-1, &rref_len);
        rfhmm->set_model(emotif.c_str(), rflank);
        rfhmm->set_mismatch_penalty(penalty);
        rfhmm->align(rrefseq, qual.c_str());
        if (debug) rfhmm->print_alignment();

        int32_t motif_beg1 = motif_end1 - (pos1+lfhmm->get_motif_read_epos1() - pos1+lfhmm->get_motif_read_spos1());

        
        variant.iregion.beg1 = motif_beg1;        
        variant.iregion.end1 = motif_end1;


        if (lflank_len) free(lflank);
        if (lref_len) free(lrefseq);
        if (rflank_len) free(rflank);
        if (rref_len) free(rrefseq);
            


//        char** candidate_motifs = suggest_motifs(bcf_get_allele(v), bcf_get_n_allele(v), no_candidate_motifs);
//
//        for (size_t i=0; i<no_candidate_motifs; ++i)
//        {
            //should not overwrite model parameters here.

            //        std::cerr << "lflank           : " << lflank << "\n";
//        std::cerr << "RU               : " << ru << "\n";
//        std::cerr << "ref_genome       : " << ref_genome << "\n";
//        std::cerr << "CANDIDATE MOTIFS : ";
//        for (size_t i=0; i<no_candidate_motifs; ++i)
//        {
//            std::cerr << (i?",":"") << candidate_motifs[i];
//        }


  //        search_flanks(chrom, start1, candidate_motifs[i]);

            //check if fit is good, if good, exit loop;
            //at this point, do not store candidates.
            //idea is to choose only really good fits.
            //
            //as to what exactly a good fit.
            //motif discordance is perfect
            //motif concordance is perfect?
            //loosen definition for long sets.
//        }
    }
}

/**
 * Pick shortest motif.
 */
std::string STRMotif::pick_motif(std::string& ref, std::string& alt)
{
    std::cerr << "setting motifs " << ref << " " << alt << "\n";
    
    std::string sequence;
    if (ref.size()==1 && alt.size()!=1)
    {
        if (ref.at(0)==alt.at(0))
        {
            sequence = alt.substr(1, alt.size()-1);
        }
        else if (ref.at(ref.size()-1)==alt.at(alt.size()-1))
        {
            sequence = alt.substr(0, alt.size()-1);
        }
    }
    else if (ref.size()!=1 && alt.size()==1)
    {
        if (ref.at(0)==alt.at(0))
        {
            std::cerr << "set sequence here\n";
            sequence = ref.substr(1, ref.size()-1);
        }
        else if (ref.at(ref.size()-1)==alt.at(alt.size()-1))
        {
            sequence = ref.substr(0, ref.size()-1);
        }
    }

    std::cerr << "checking " << sequence << "\n";

    if (sequence.size()>max_len)
    {
        initialize_factors(sequence.size()+1);
    }

    size_t i = 0;
    size_t sub_motif_len;
    size_t len = sequence.size();
    const char* seq = sequence.c_str();
    while ((sub_motif_len=factors[sequence.size()][i])!=len)
    {
        size_t n_sub_motif = len/sub_motif_len;

        bool exact = true;

        for (size_t j=0; j<n_sub_motif; ++j)
        {
            if (strncmp(&seq[0], &seq[j*sub_motif_len], sub_motif_len))
            {
                exact = false;
                break;
            }
        }

        if (exact) break;
        ++i;
    }

    return sequence.substr(0, sub_motif_len);
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
 * Trim alleles.
 */
void STRMotif::trim(int32_t& pos1, std::string& ref, std::string& alt)
{
    while (true)
    {
        if (ref.size()==1 || alt.size()==1)
        {
            break;
        }
        else if (ref.at(0)!=alt.at(0) && ref.at(ref.size()-1)!=alt.at(alt.size()-1))
        {
            break;
        }
        else
        {
            if (ref.at(ref.size()-1)==alt.at(alt.size()-1))
            {
                ref.erase(ref.size()-1,1);
                alt.erase(alt.size()-1,1);
            }
            else if (ref.at(0)==alt.at(0))
            {
                ref.erase(0,1);
                alt.erase(0,1);
                ++pos1;
            }
        }
    }
}

/**
 * Left align alleles.
 */
void STRMotif::left_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt)
{
    int32_t seq_len;
    char* seq;
    while (ref.at(ref.size()-1)==alt.at(alt.size()-1) && pos1>1)
    {
        seq = faidx_fetch_seq(fai, chrom, pos1-2, pos1-2, &seq_len);
        if (seq_len)
        {
            ref.erase(ref.size()-1,1);
            alt.erase(alt.size()-1,1);
            ref.insert(0, 1, seq[0]);
            alt.insert(0, 1, seq[0]);
            free(seq);
            --pos1;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }
}

/**
 * Right align alleles.
 */
void STRMotif::right_align(const char* chrom, int32_t& pos1, std::string& ref, std::string& alt)
{
    int32_t seq_len;
    char* seq;
    while (ref.at(0)==alt.at(0))
    {
        seq = faidx_fetch_seq(fai, chrom, pos1, pos1, &seq_len);
        if (seq_len)
        {
            ref.erase(0,1);
            alt.erase(0,1);
            ref.push_back(seq[0]);
            alt.push_back(seq[0]);
            free(seq);
            ++pos1;
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] Cannot read from sequence file\n", __FILE__,__LINE__,__FUNCTION__);
            exit(1);
        }
    }
}

/**
 * Detect allele lower bound extent.
 */
void STRMotif::detect_lower_bound_allele_extent(const char* chrom, int32_t& pos1, std::vector<std::string>& alleles, int32_t& start1, int32_t& end1)
{
    if (alleles.size()==2)
    {
        std::string ref = alleles[0];
        std::string alt = alleles[1];

        trim(pos1, ref, alt);

    }
    else
    {
    }
}

/**
 * Detect candidate flanks given a motif fit.
 * Update model atttributes.
 */
void STRMotif::search_flanks(const char* chrom, int32_t start1, char* motif)
{
    //given chromosome position and motif
    //attempt to


    ////////////////////
    //detect right flank
    ////////////////////

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

    //extract rflank
    //pick flank
    lflank = faidx_fetch_uc_seq(fai, chrom, start1-10, start1-1, &lflank_len);

    //rfhmm
    rfhmm->set_model(motif, rflank);
    rfhmm->set_mismatch_penalty(5);
    rfhmm->align(ref, qual.c_str());
    rfhmm->print_alignment();

    ///////////////////
    //detect left flank
    ///////////////////
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

    //////////////////////////////
    //annotate STR characteristics
    //////////////////////////////

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