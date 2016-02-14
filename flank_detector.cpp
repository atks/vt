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
    float delta = 0.0001;
    float epsilon = 0.0005;
    float tau = 0.01;
    float eta = 0.01;
    float mismatch_penalty = 3;

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
 */
void FlankDetector::detect_flanks(bcf_hdr_t* h, bcf1_t *v, Variant& variant, uint32_t mode)
{
    VNTR& vntr = variant.vntr;

    //simple single base pair clipping of ends
    if (mode==CLIP_ENDS)
    {
        if (debug)
        {
            std::cerr << "********************************************\n";
            std::cerr << "CLIP ENDS\n";
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
        //
        // RU : T
        // repeat_tract : G[T]C where T is an insert
        else if (vntr.exact_repeat_tract.size()==2)
        {

        }

        vntr.ru = choose_repeat_unit(vntr.exact_repeat_tract, vntr.motif);
        ahmm->set_model(vntr.ru.c_str());
        ahmm->align(vntr.exact_repeat_tract.c_str(), qual.c_str());

        vntr.exact_score = ahmm->motif_concordance;
        vntr.exact_trf_score = ahmm->trf_score;
        vntr.exact_no_exact_ru = ahmm->exact_motif_count;
        vntr.exact_total_no_ru = ahmm->motif_count;
        vntr.exact_rl = vntr.exact_repeat_tract.size();
        vntr.exact_ll = vntr.exact_rl + variant.max_dlen;
        
        compute_composition_and_entropy(vntr.exact_repeat_tract);
        vntr.exact_comp[0] = comp[0];
        vntr.exact_comp[1] = comp[1];
        vntr.exact_comp[2] = comp[2];
        vntr.exact_comp[3] = comp[3];
        vntr.exact_entropy = entropy;
            
        if (debug)
        {
            std::cerr << "\n";
            std::cerr << "repeat_tract              : " << vntr.exact_repeat_tract << "\n";
            std::cerr << "position                  : [" << vntr.exact_beg1 << "," << vntr.exact_end1 << "]\n";
            std::cerr << "score                     : " << vntr.exact_score << "\n";
            std::cerr << "trf score                 : " << vntr.exact_trf_score << "\n";
            std::cerr << "repeat units              : " << vntr.exact_rl << "\n";
            std::cerr << "exact repeat units        : " << vntr.exact_no_exact_ru << "\n";
            std::cerr << "total no. of repeat units : " << vntr.exact_total_no_ru << "\n";
            std::cerr << "\n";
        }
    }
    else if (mode==FRAHMM)
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

        char* rflank = NULL;
        char* lflank = NULL;

        int32_t lflank_end1;
        int32_t rflank_beg1;

        char* seq;
        int32_t seq_len;
//        bool encountered_N = false;

        while (true)
        {
            //pick 5 bases to the right
            rflank = rs->fetch_seq(variant.chrom.c_str(), vntr.exact_end1+1, vntr.exact_end1+5);

            //pick 105 bases for aligning

            seq = rs->fetch_seq(variant.chrom.c_str(), vntr.exact_end1-slen, vntr.exact_end1+5);


            rfhmm->set_model(vntr.ru.c_str(), rflank);
            rfhmm->align(seq, qual.c_str());
            if (debug) rfhmm->print_alignment();

            if (rflank) free(rflank);
            if (seq) free(seq);

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
            if (rfhmm->get_lflank_read_epos1()>2*vntr.ru.size())
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

        slen = 100;

        //pick 5 bases to right
        while(true)
        {
            lflank = rs->fetch_seq(variant.chrom.c_str(), lflank_end1-5, lflank_end1);

            //pick 105 bases for aligning
            char* seq = rs->fetch_seq(variant.chrom.c_str(), lflank_end1-5, lflank_end1+slen-1);

            lfhmm->set_model(lflank, vntr.ru.c_str());
            lfhmm->align(seq, qual.c_str());
            if (debug) lfhmm->print_alignment();

            if (seq_len) free(seq);

            if (lfhmm->get_rflank_read_epos1()!=INT32_MAX)
            {
                rflank_beg1 = lflank_end1 - 5 + lfhmm->get_rflank_read_spos1() - 1;
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
        if (lflank) free(lflank);

        lflank = rs->fetch_seq(variant.chrom.c_str(), lflank_end1-10, lflank_end1);
        rflank = rs->fetch_seq(variant.chrom.c_str(), rflank_beg1, rflank_beg1 +10-1);

        vntr.fuzzy_beg1 = lflank_end1+1;
        vntr.fuzzy_end1 = rflank_beg1-1;
        int32_t repeat_tract_len;
        char* repeat_tract = rs->fetch_seq(variant.chrom.c_str(), lflank_end1, rflank_beg1-1-1);
        vntr.fuzzy_repeat_tract.assign(repeat_tract);
        if (repeat_tract) free(repeat_tract);
        vntr.fuzzy_score = lfhmm->motif_concordance;
        vntr.fuzzy_trf_score = lfhmm->trf_score;
        vntr.fuzzy_no_exact_ru = lfhmm->exact_motif_count;
        vntr.fuzzy_total_no_ru = lfhmm->motif_count;
        vntr.fuzzy_rl = rflank_beg1-lflank_end1-1;
        vntr.fuzzy_ll = vntr.fuzzy_rl + variant.max_dlen;
        
        if (lflank) free(lflank);
        if (rflank) free(rflank);

        compute_composition_and_entropy(vntr.fuzzy_repeat_tract);
        vntr.fuzzy_comp[0] = comp[0];
        vntr.fuzzy_comp[1] = comp[1];
        vntr.fuzzy_comp[2] = comp[2];
        vntr.fuzzy_comp[3] = comp[3];
        vntr.fuzzy_entropy = entropy;

        if (debug)
        {
            std::cerr << "\n";
            vntr.print();
            std::cerr << "\n";
        }
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
 * Chooses a phase of the motif that is appropriate for the alignment.
 * This differs from choose_exact_repeat_unit() where the motif is returned
 * if not suitable repeat unit is found.
 */
std::string FlankDetector::choose_repeat_unit(std::string& seq, std::string& motif)
{
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = shift_str(motif, i);
        if (seq.compare(0, smotif.size(), smotif)==0)
        {
            return smotif;
        }
    }

    //should return empty string
    return motif;
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment.
 * This returns the empty string if the motif does not have an exact
 * match in all its phases.
 */
std::string FlankDetector::choose_exact_repeat_unit(std::string& seq, std::string& motif)
{
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
 */
void FlankDetector::compute_purity_score(Variant& variant, int32_t amode)
{
    VNTR& vntr = variant.vntr;
     
    if (amode&FINAL)
    {
        compute_purity_score(vntr.repeat_tract, vntr.motif);
        
        vntr.ru = ru;
        vntr.score = score;
        vntr.trf_score = trf_score;
        vntr.no_exact_ru = no_exact_ru;
        vntr.total_no_ru = total_no_ru;
        vntr.rl = rl;
        vntr.ll = rl + variant.max_dlen;
    }
    
    if (amode&EXACT)
    {    
        compute_purity_score(vntr.exact_repeat_tract, vntr.motif);
        
        vntr.ru = ru;
        vntr.exact_score = score;
        vntr.exact_trf_score = trf_score;
        vntr.exact_no_exact_ru = no_exact_ru;
        vntr.exact_total_no_ru = total_no_ru;
        vntr.exact_rl = rl;
        vntr.exact_ll = rl + variant.max_dlen;
    }
    
    if (amode&FUZZY)
    {
        compute_purity_score(vntr.fuzzy_repeat_tract, vntr.motif);
        
        vntr.ru = ru;
        vntr.fuzzy_score = score;
        vntr.fuzzy_trf_score = trf_score;
        vntr.fuzzy_no_exact_ru = no_exact_ru;
        vntr.fuzzy_total_no_ru = total_no_ru;
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
    
    float exact_motif_count = 0;
    uint32_t mlen = ru.size();

    if (ru!="")
    {
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
            no_exact_ru = repeat_tract.size()/motif.size();
            total_no_ru = no_exact_ru;
            rl = (float) repeat_tract.size()/(float) motif.size();            
            trf_score = repeat_tract.size() << 2;
            return;
        }
    }
    else
    {
        ru = motif;
    }

    //fall through to computing inexact sequence
    ahmm->set_model(ru.c_str());
    ahmm->align(repeat_tract.c_str(), qual.c_str());

    score = ahmm->motif_concordance;
    no_exact_ru = ahmm->exact_motif_count;
    total_no_ru = ahmm->motif_count;
    rl = repeat_tract.size();
    trf_score = ahmm->trf_score;
}

/**
 * Computes composition and entropy ofrepeat tract.
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
    entropy = -entropy;
    entropy = std::round(100*entropy)/100; 
        
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
//    std::cerr << "entropy: " << entropy << "\n";    

    entropy2 = 0;
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
            }
        }
        entropy2 = -entropy2;
        entropy2 = std::round(100*entropy2)/100; 
   
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
//        std::cerr << "entropy2: " << entropy2 << "\n";   
    }
        
}