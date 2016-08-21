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

#include "flank_detector.h"

/**
 * Constructor.
 */
FlankDetector::FlankDetector(std::string& ref_fasta_file, bool debug)
{
    //////////////////////
    //initialize variables
    //////////////////////
    this->debug = debug;

    ///////////////////
    //initialize raHMMs
    ///////////////////
//    float delta = 0.0001;
//    float epsilon = 0.0005;
//    float tau = 0.01;
//    float eta = 0.01;
//    float mismatch_penalty = 3;
    float delta = 0.0000001;
    float epsilon = 0.0000001;
    float tau = 0.01;
    float eta = 0.01;
    float mismatch_penalty = 5;

    ahmm = new AHMM(false);
    ahmm->set_delta(delta);
    ahmm->set_epsilon(epsilon);
    ahmm->set_tau(tau);
    ahmm->set_eta(eta);
    ahmm->set_mismatch_penalty(mismatch_penalty);
    ahmm->initialize_T();

    lfhmm = new LFHMM(false);
    lfhmm->set_delta(delta);
    lfhmm->set_epsilon(epsilon);
    lfhmm->set_tau(tau);
    lfhmm->set_eta(eta);
    lfhmm->set_mismatch_penalty(mismatch_penalty);
    lfhmm->initialize_T();

    rfhmm = new RFHMM(false);
    rfhmm->set_delta(delta);
    rfhmm->set_epsilon(epsilon);
    rfhmm->set_tau(tau);
    rfhmm->set_eta(eta);
    rfhmm->set_mismatch_penalty(mismatch_penalty);
    rfhmm->initialize_T();

    qual.assign(1024, 'K');

    //////////////////
    //initialize tools
    //////////////////
    rs = new ReferenceSequence(ref_fasta_file);
};

/**
 * Destructor.
 */
FlankDetector::~FlankDetector()
{
    delete rs;
}

/**
 * Detect repeat region.
 *
 * updates
 * beg1
 * end1
 * repeat_tract
 */
void FlankDetector::detect_flanks(Variant& variant, uint32_t mode)
{
    VNTR& vntr = variant.vntr;

    //simple single base pair clipping of ends
    if (mode==EXACT)
    {
        if (debug)
        {
            std::cerr << "********************************************\n";
            std::cerr << "EXACT\n";
            std:: cerr << "\n";
        }

        if (debug)
        {
            std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
            std::cerr << "Detect and remove anchor base\n";
        }

        if (vntr.exact_repeat_tract.size()>2)
        {
            //removing the anchor bases
            if (vntr.mlen==1)
            {
                int32_t offset = 0;
                int32_t length = vntr.exact_repeat_tract.size();
                if (vntr.exact_repeat_tract.at(0)!=vntr.motif.at(0))
                {
                    offset = 1;
                    ++vntr.exact_beg1;
                }

                if (vntr.exact_repeat_tract.at(vntr.exact_repeat_tract.size()-1)!=vntr.motif.at(0))
                {
                    length -= offset+1;
                    --vntr.exact_end1;
                }

                vntr.exact_repeat_tract = vntr.exact_repeat_tract.substr(offset, length);
            }
            else
            {
                //this works only for simple indels
                if (vntr.exact_repeat_tract.size()>=3)
                {
                    vntr.exact_repeat_tract = vntr.exact_repeat_tract.substr(1, vntr.exact_repeat_tract.size()-2);
                    ++vntr.exact_beg1;
                    --vntr.exact_end1;
                }
            }
        }
        //this is for nonexistent repeat units
        // RU : T
        // repeat_tract : G[T]C where T is an insert
        else if (vntr.exact_repeat_tract.size()==2)
        {

        }

        if (debug)
        {
            std::cerr << "\n";
            std::cerr << "repeat_tract              : " << vntr.exact_repeat_tract << "\n";
            std::cerr << "position                  : [" << vntr.exact_beg1 << "," << vntr.exact_end1 << "]\n";
            std::cerr << "score                     : " << vntr.exact_score << "\n";
            std::cerr << "trf score                 : " << vntr.exact_trf_score << "\n";
            std::cerr << "repeat units              : " << vntr.exact_rl << "\n";
            std::cerr << "perfect repeat units      : " << vntr.exact_no_perfect_ru << "\n";
            std::cerr << "no. of repeat units       : " << vntr.exact_no_ru << "\n";
            std::cerr << "\n";
        }
    }
    else if (mode==FUZZY)
    {
        if (debug)
        {
            std::cerr << "********************************************\n";
            std::cerr << "DETECTING REPEAT TRACT FUZZILY\n";
        }

        ///////////////////////
        //fuzzy right alignment
        ///////////////////////
        if (debug)
        {
            std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
            std::cerr << "Fuzzy right alignment\n";
        }

        int32_t slen = 100;

        std::string rflank;
        std::string lflank;

        int32_t lflank_end1;
        int32_t rflank_beg1;

        std::string seq;
        int32_t seq_len;

        while (true)
        {
            //fetch sequences for modeling
            rs->fetch_seq(variant.chrom, vntr.exact_end1+1, vntr.exact_end1+5, rflank);
            rs->fetch_seq(variant.chrom, vntr.exact_end1-slen, vntr.exact_end1, seq);
            vntr.fuzzy_ru = choose_fuzzy_3prime_repeat_unit(seq, vntr.fuzzy_motif);
            seq.append(rflank);

            rfhmm->set_model(vntr.fuzzy_ru.c_str(), rflank.c_str());
            rfhmm->align(seq.c_str(), qual.c_str());
            if (debug) rfhmm->print_alignment();

            //////////////////////
            //fuzzy left alignment
            //////////////////////
            if (debug)
            {
                std::cerr << "\n";
                std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
                std::cerr << "Fuzzy left alignment\n";
            }

            //this is a hack around rfhmm rigidity in modeling the RUs
            //todo: we should change this to a reverse version of LFHMM!!!!
            if (rfhmm->get_lflank_read_epos1()>std::min((int32_t)(10*vntr.fuzzy_ru.size()), 50))
            {
                lflank_end1 = vntr.exact_end1-slen-1+1 + rfhmm->get_lflank_read_epos1() - 1;
                break;
            }
            else if (slen==1000)
            {
                lflank_end1 = vntr.exact_end1 - 1000 - 1;
                vntr.is_large_repeat_tract = true;
                break;
            }
            else
            {
                slen +=100;

                if (debug)
                    std::cerr << "extending the reference sequence for RFHMM : " << slen << "\n";
            }
        }

//        slen = 100;

        //pick 5 bases to right
        while(true)
        {
            //fetch sequences for modeling
            rs->fetch_seq(variant.chrom, lflank_end1-5+1, lflank_end1, lflank);
            rs->fetch_seq(variant.chrom, lflank_end1+1, lflank_end1+slen, seq);
            vntr.fuzzy_ru = choose_fuzzy_5prime_repeat_unit(seq, vntr.fuzzy_motif);
            seq =  lflank + seq;

            lfhmm->set_model(lflank.c_str(), vntr.fuzzy_ru.c_str());
            lfhmm->align(seq.c_str(), qual.c_str());
            if (debug) lfhmm->print_alignment();

//            if (lfhmm->get_rflank_read_epos1()!=INT32_MAX ||
            if (lfhmm->get_rflank_read_spos1()<slen-std::min((int32_t)(10*vntr.fuzzy_ru.size()), 50))
            {
                rflank_beg1 = lflank_end1 - 5 + lfhmm->get_rflank_read_spos1();
                break;
            }
            else if (slen==1000)
            {
                rflank_beg1 = lflank_end1 + 1000;
                vntr.is_large_repeat_tract = true;
                break;
            }
            else
            {
                slen +=100;
                if (debug)
                    std::cerr << "extending the reference sequence for LFHMM : " << slen << "\n";
            }
        }

        vntr.fuzzy_beg1 = lflank_end1+1;
        vntr.fuzzy_end1 = rflank_beg1-1;
        rs->fetch_seq(variant.chrom, vntr.fuzzy_beg1, vntr.fuzzy_end1, vntr.fuzzy_repeat_tract);
    }
};

/**
 * Shifts a string.
 */
std::string FlankDetector::shift_str(std::string& seq, uint32_t i)
{
    std::string sseq = seq;
    if (i)
    {
        sseq = seq.substr(i) + seq.substr(0,i);
    }

    return sseq;
}

/**
 * Score string.
 */
int32_t FlankDetector::compute_score(int32_t start, int32_t len, std::string& a, std::string& b)
{
    if (len>b.size())
    {
        len = b.size();
    }

    int32_t score = 0;
    int32_t i = start;
    int32_t j = 0;
    while (j<len)
    {
        if (a.at(i)==b.at(j))
        {
            score += 2;
        }
        else
        {
            score -=7;
        }

        ++i;
        ++j;
    }

    return score;
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment from the 5 prime end.
 * This differs from choose_exact_repeat_unit() where the motif is returned
 * if not suitable repeat unit is found.
 */
std::string FlankDetector::choose_5prime_repeat_unit(std::string& seq, std::string& motif)
{
//    std::cerr << "seq   : " << seq << "\n";
//    std::cerr << "motif : " << motif << "\n";

    int32_t best_score = -10000;
    std::string best_motif = motif;
    int32_t max_score = 2*motif.size();

    int32_t mlen = motif.size();
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = shift_str(motif, i);

        if (seq.compare(0, mlen, smotif)==0)
        {
            return smotif;
        }

        if (seq.compare(0, mlen, smotif)==0)
        {
            return smotif;
        }

        std::string rc_smotif = VNTR::reverse_complement(smotif);
        if (seq.compare(0, mlen, rc_smotif)==0)
        {
            return rc_smotif;
        }
    }

    //if no exact match, perform best fit.
    //use a priority queue here

    return motif;
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment from the 5 prime end.
 * If no exact match is available, the best possible match is returned.
 */
std::string FlankDetector::choose_fuzzy_5prime_repeat_unit(std::string& seq, std::string& motif)
{
    int32_t best_score = -10000;
    std::string best_motif = motif;
    int32_t max_score = 2*motif.size();

    int32_t mlen = motif.size();
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = shift_str(motif, i);
        int32_t sscore = compute_score(0, mlen, seq, smotif);
        if (sscore==max_score)
        {
            return smotif;
        }
        else if (sscore>best_score)
        {
            best_score = sscore;
            best_motif = smotif;
        }

        std::string rc_smotif = VNTR::reverse_complement(smotif);
        int32_t rc_score = compute_score(0, mlen, seq, rc_smotif);
        if (rc_score==max_score)
        {
            return rc_smotif;
        }
        else if (rc_score>best_score)
        {
            best_score = rc_score;
            best_motif = rc_smotif;
        }
    }

    return best_motif;
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment from the 3 prime end.
 * This differs from choose_exact_repeat_unit() where the motif is returned
 * if not suitable repeat unit is found.
 */
std::string FlankDetector::choose_3prime_repeat_unit(std::string& seq, std::string& motif)
{
    int32_t mlen = motif.size();
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = shift_str(motif, i);
        if (seq.compare(seq.size()-mlen, mlen, smotif)==0)
        {
            return smotif;
        }

        std::string rc_smotif = VNTR::reverse_complement(smotif);
        if (seq.compare(seq.size()-mlen, mlen, rc_smotif)==0)
        {
            return rc_smotif;
        }
    }

    return motif;
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment from the 3 prime end.
 * If no exact match is available, the best possible match is returned.
 */
std::string FlankDetector::choose_fuzzy_3prime_repeat_unit(std::string& seq, std::string& motif)
{
    int32_t best_score = -10000;
    std::string best_motif = motif;
    int32_t max_score = 2*motif.size();

    int32_t mlen = motif.size();
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = shift_str(motif, i);
        int32_t sscore = compute_score(seq.size()-mlen+1, mlen, seq, smotif);
        if (sscore==max_score)
        {
            return smotif;
        }
        else if (sscore>best_score)
        {
            best_score = sscore;
            best_motif = smotif;
        }

        std::string rc_smotif = VNTR::reverse_complement(smotif);
        int32_t rc_score = compute_score(seq.size()-mlen+1, mlen, seq, rc_smotif);
        if (rc_score==max_score)
        {
            return rc_smotif;
        }
        else if (rc_score>best_score)
        {
            best_score = rc_score;
            best_motif = rc_smotif;
        }
    }

    return best_motif;
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment.
 * This returns the empty string if the motif does not have an exact
 * match in all its phases.
 */
std::string FlankDetector::choose_exact_repeat_unit(std::string& seq, std::string& motif)
{
//    if (debug)
//    {
//        std::cerr << "choose_repeat_unit\n";
//        std::cerr << "\tseq    " << seq << "\n";
//        std::cerr << "\tmotif  " << motif << "\n";
//    }

    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = shift_str(motif, i);

        if (seq.compare(0, smotif.size(), smotif)==0)
        {
            return smotif;
        }
    }

    return "";
}

/**
 * Polish repeat tract ends.
 */
void FlankDetector::polish_repeat_tract_ends(Variant& variant)
{
}

/**
 * Polish repeat tract ends.
 *
 *
 */
void FlankDetector::polish_repeat_tract_ends(std::string& repeat_tract, std::string& motif, bool debug)
{
    if (debug)
    {
        std::cerr << "===================\n";
        std::cerr << "Polish repeat tract\n";
        std::cerr << "===================\n";
        std::cerr << "repeat tract : " << repeat_tract << " (" << repeat_tract.size() << ")\n";
        std::cerr << "motif        : " << motif  << "\n";
    }

    min_beg0 = repeat_tract.size();
    max_end0 = -1;
    int32_t mlen = motif.size();
    int32_t rlen = repeat_tract.size();

    //todo:  we can use a FSA for substring matching.
    //is there a way check all phases simultaneously?
    for (uint32_t i = 0; i<mlen; ++i)
    {
        std::string smotif = shift_str(motif, i);

        int32_t temp_min_beg0 = 0;
        int32_t temp_max_end0 = 0;

        for (int32_t j=0; j<rlen; ++j)
        {
            if (repeat_tract.compare(j, mlen, smotif)==0)
            {
                temp_min_beg0 = j;
                min_beg0 = std::min(j, min_beg0);
                break;
            }
        }

        for (int32_t j=rlen-mlen; j>=0; --j)
        {
            if (repeat_tract.compare(j, mlen, smotif)==0)
            {
                temp_max_end0 = j+mlen-1;
                max_end0 = std::max(j+mlen-1, max_end0);
                break;
            }
        }

        if (debug) std::cerr << "\t" << smotif  << " " << temp_min_beg0 << " " << temp_max_end0 << "\n";
    }

    polished_repeat_tract = repeat_tract.substr(min_beg0, max_end0-min_beg0+1);

    if (debug)
    {
        std::cerr << "min beg      : " << min_beg0 << "\n";
        std::cerr << "max end      : " << max_end0 << "\n";
        std::cerr << "polished     : " << polished_repeat_tract << "\n";
    }
}

/**
 * Computes purity score of a sequence with respect to a motif.
 *
 * updates
 * 1. ru
 * 2. score
 * 3. trf_score
 * 4. no_exact_ru
 * 5. total_no_ru
 * 6. ref
 * 7. rl
 * 8. ll
 */
void FlankDetector::compute_purity_score(Variant& variant, int32_t amode)
{
    VNTR& vntr = variant.vntr;

    if (amode&FINAL)
    {
        compute_purity_score(vntr.repeat_tract, vntr.ru);

        vntr.ru = ru;
        vntr.score = score;
        vntr.trf_score = trf_score;
        vntr.no_perfect_ru = no_perfect_ru;
        vntr.no_ru = no_ru;
        vntr.ref = ref;
        vntr.rl = rl;
        vntr.ll = rl + variant.max_dlen;
    }

    if (amode&EXACT)
    {
        compute_purity_score(vntr.exact_repeat_tract, vntr.exact_ru);

        vntr.exact_ru = ru;
        vntr.exact_score = score;
        vntr.exact_trf_score = trf_score;
        vntr.exact_no_perfect_ru = no_perfect_ru;
        vntr.exact_no_ru = no_ru;
        vntr.exact_ref = ref;
        vntr.exact_rl = rl;
        vntr.exact_ll = rl + variant.max_dlen;
    }

    if (amode&FUZZY)
    {
        compute_purity_score(vntr.fuzzy_repeat_tract, vntr.fuzzy_ru);

        vntr.fuzzy_ru = ru;
        vntr.fuzzy_score = score;
        vntr.fuzzy_trf_score = trf_score;
        vntr.fuzzy_no_perfect_ru = no_perfect_ru;
        vntr.fuzzy_no_ru = no_ru;
        vntr.fuzzy_ref = ref;
        vntr.fuzzy_rl = rl;
        vntr.fuzzy_ll = rl + variant.max_dlen;
    }
}

/**
 * Computes purity score of a sequence with respect to a motif.
 */
void FlankDetector::compute_purity_score(std::string& repeat_tract, std::string& motif)
{
    ru = choose_exact_repeat_unit(repeat_tract, motif);

    if (debug)
    {
        std::cerr << "ru           " << ru << "\n";
        std::cerr << "repeat tract " << repeat_tract << "\n";

    }

    ///////////////////
    //exact calculation
    ///////////////////
    if (ru!="")
    {
        uint32_t mlen = ru.size();
        uint32_t j=0;
        bool exact = true;
        for (uint32_t i=0; i<repeat_tract.size(); ++i)
        {
            if (ru.at(j)!=repeat_tract.at(i))
            {
                exact = false;
                break;
            }

            j = (j==mlen-1) ? 0 : j+1;
        }

        if (exact)
        {
            score = 1;
            no_perfect_ru = repeat_tract.size()/motif.size();
            no_ru = no_perfect_ru;
            ref = (float) repeat_tract.size()/motif.size();
            rl = repeat_tract.size();
            trf_score = repeat_tract.size() << 1;
            return; //done!
        }
    }

    ///////////////////
    //fuzzy calculation
    ///////////////////
    ru = motif;

    if (ru.size()>ahmm->max_len)
    {
        //compute by chunks
        //todo:: not the best way. update with a localized aligner

    }
    else
    {
        ahmm->set_model(ru.c_str());
        ahmm->align(repeat_tract.c_str(), qual.c_str());

        score = ahmm->motif_concordance;
        score = std::round(100*score)/100;
        no_perfect_ru = ahmm->exact_motif_count;
        no_ru = ahmm->motif_count;
        ref = ahmm->frac_no_repeats;
        rl = repeat_tract.size();
        trf_score = ahmm->trf_score;
    }
}

/**
 * Computes composition and entropy ofrepeat tract.
 *
 * updates
 * 1. comp
 * 2. entropy
 * 3. entropy2
 * 4. kl_divergence
 * 5. kl_divergence2
 */
void FlankDetector::compute_composition_and_entropy(Variant& variant, int32_t amode)
{
    VNTR& vntr = variant.vntr;

    if (amode&FINAL)
    {
        compute_composition_and_entropy(vntr.repeat_tract);
        vntr.comp[0] = comp[0];
        vntr.comp[1] = comp[1];
        vntr.comp[2] = comp[2];
        vntr.comp[3] = comp[3];

        vntr.entropy = entropy;
        vntr.entropy2 = entropy2;
        vntr.kl_divergence = kl_divergence;
        vntr.kl_divergence2 = kl_divergence2;
    }
    else if (amode&EXACT)
    {
        compute_composition_and_entropy(vntr.exact_repeat_tract);
        vntr.exact_comp[0] = comp[0];
        vntr.exact_comp[1] = comp[1];
        vntr.exact_comp[2] = comp[2];
        vntr.exact_comp[3] = comp[3];

        vntr.exact_entropy = entropy;
        vntr.exact_entropy2 = entropy2;
        vntr.exact_kl_divergence = kl_divergence;
        vntr.exact_kl_divergence2 = kl_divergence2;
    }
    else if (amode&FUZZY)
    {
        compute_composition_and_entropy(vntr.fuzzy_repeat_tract);
        vntr.fuzzy_comp[0] = comp[0];
        vntr.fuzzy_comp[1] = comp[1];
        vntr.fuzzy_comp[2] = comp[2];
        vntr.fuzzy_comp[3] = comp[3];
        vntr.fuzzy_entropy = entropy;
        vntr.fuzzy_entropy2 = entropy2;
        vntr.fuzzy_kl_divergence = kl_divergence;
        vntr.fuzzy_kl_divergence2 = kl_divergence2;
    }
}

/**
 * Computes composition and entropy ofrepeat tract.
 */
void FlankDetector::compute_composition_and_entropy(std::string& repeat_tract)
{
    int32_t aux_comp[4] = {0,0,0,0};
    int32_t aux_comp2[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int32_t b2i[10] = {0,1,0,2,0,0,0,0,0,3};

//    ACGT x ACGT
//    AA - 0    = 0*0 + 0
//    AC - 1
//    AG - 2
//    AT - 3
//    CA - 4    = 1*4 + 0 = 1<<2
//    CC - 5
//    CG - 6
//    CT - 7
//    GA - 8    = 2*4 + 0  = 2<<1 = 10=>100
//    GC - 9
//    GG - 10
//    GT - 11
//    TA - 12   = 3*4 + 0  = 3<<1  11=>110=>1100
//    TC - 13
//    TG - 14
//    TT - 15    = 3*4+3 = 15

    int32_t n = repeat_tract.size();

    for (uint32_t i=0; i<n; ++i)
    {
        uint32_t b1 = b2i[(repeat_tract.at(i)-65)>>1];

        ++aux_comp[b1];
        if (i<n-1)
        {
            uint32_t b2 = b2i[(repeat_tract.at(i+1)-65)>>1];
            uint32_t bb = (b1<<2) + b2;
            ++aux_comp2[bb];
        }
    }

    float p[4] = {(float) aux_comp[0]/n, (float) aux_comp[1]/n, (float) aux_comp[2]/n, (float) aux_comp[3]/n};

    entropy = 0;
    if (p[0]) entropy += p[0] * std::log2(p[0]);
    if (p[1]) entropy += p[1] * std::log2(p[1]);
    if (p[2]) entropy += p[2] * std::log2(p[2]);
    if (p[3]) entropy += p[3] * std::log2(p[3]);
    kl_divergence = p[0] + p[1] + p[2] + p[3];
    kl_divergence = entropy + 2*kl_divergence;
    kl_divergence = std::round(100*kl_divergence)/100;
    kl_divergence = fix_neg_zero(kl_divergence);
    entropy = -entropy;
    entropy =  std::round(100*entropy)/100;
    entropy =  fix_neg_zero(entropy);

    comp[0] = std::round(p[0] * 100);
    comp[1] = std::round(p[1] * 100);
    comp[2] = std::round(p[2] * 100);
    comp[3] = std::round(p[3] * 100);

//    std::cerr << "tract: " << repeat_tract << "\n";
//    std::cerr << "A: " << comp[0] << " " << aux_comp[0] << "\n";
//    std::cerr << "C: " << comp[1] << " " << aux_comp[1] << "\n";
//    std::cerr << "G: " << comp[2] << " " << aux_comp[2] << "\n";
//    std::cerr << "T: " << comp[3] << " " << aux_comp[3] << "\n";
//    std::cerr << "\n";
//    std::cerr << "entropy       : " << entropy << "\n";
//    std::cerr << "kl_divergence : " << kl_divergence << "\n";

    entropy2 = 0;
    kl_divergence2 = 0;
    float log2q = -4;
    if (n!=1)
    {
        float p2[16];
        for (uint32_t i=0; i<16; ++i)
        {
            p2[i] = (float)aux_comp2[i]/(n-1);
        }

        for (uint32_t i=0; i<16; ++i)
        {
            if (p2[i])
            {
                entropy2 += p2[i]* std::log2(p2[i]);
                kl_divergence2 += p2[i];
            }
        }
        kl_divergence2 = entropy2 + 4*kl_divergence2;
        kl_divergence2 = std::round(100*kl_divergence2)/100;
        kl_divergence2 = fix_neg_zero(kl_divergence2);
        entropy2 = -entropy2;
        entropy2 = std::round(100*entropy2)/100;
        entropy2 = fix_neg_zero(entropy2);

//        std::cerr << "tract: " << repeat_tract << "\n";
//        std::cerr << "AA: " << aux_comp2[0] << " " << p2[0] << "\n";
//        std::cerr << "AC: " << aux_comp2[1] << " " << p2[1] << "\n";
//        std::cerr << "AG: " << aux_comp2[2] << " " << p2[2] << "\n";
//        std::cerr << "AT: " << aux_comp2[3] << " " << p2[3]  << "\n";
//        std::cerr << "CA: " << aux_comp2[4] << " " << p2[4] << "\n";
//        std::cerr << "CC: " << aux_comp2[5] << " " << p2[5] << "\n";
//        std::cerr << "CG: " << aux_comp2[6] << " " << p2[6] << "\n";
//        std::cerr << "CT: " << aux_comp2[7] << " " << p2[7] << "\n";
//        std::cerr << "GA: " << aux_comp2[8] << " " << p2[8] << "\n";
//        std::cerr << "GC: " << aux_comp2[9] << " " << p2[9] << "\n";
//        std::cerr << "GG: " << aux_comp2[10] << " " << p2[10] << "\n";
//        std::cerr << "GT: " << aux_comp2[11] << " " << p2[11] << "\n";
//        std::cerr << "TA: " << aux_comp2[12] << " " << p2[12] << "\n";
//        std::cerr << "TC: " << aux_comp2[13] << " " << p2[13] << "\n";
//        std::cerr << "TG: " << aux_comp2[14] << " " << p2[14] << "\n";
//        std::cerr << "TT: " << aux_comp2[15] << " " << p2[15] << "\n";
//        std::cerr << "\n";
//        std::cerr << "entropy2       : " << entropy2 << "\n";
//        std::cerr << "kl_divergence2 : " << kl_divergence2 << "\n";
    }

}