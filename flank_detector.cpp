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

        /////////////////
        //exact alignment
        /////////////////
        if (debug)
        {
            std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
            std::cerr << "Exact left/right alignment\n";
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

        vntr.exact_motif_concordance = ahmm->motif_concordance;
        vntr.exact_no_exact_ru = ahmm->exact_motif_count;
        vntr.exact_total_no_ru = ahmm->motif_count;
        vntr.exact_rl = ahmm->repeat_tract_len;

        if (debug)
        {
            std::cerr << "\n";
//            std::cerr << "repeat_tract              : " << vntr.exact_repeat_tract << "\n";
            std::cerr << "position                  : [" << vntr.exact_beg1 << "," << vntr.exact_end1 << "]\n";
            std::cerr << "motif_concordance         : " << vntr.exact_motif_concordance << "\n";
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

        /////////////////
        //exact alignment
        /////////////////
        if (debug)
        {
            std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
            std::cerr << "Exact left/right alignment\n";
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

        vntr.exact_motif_concordance = ahmm->motif_concordance;
        vntr.exact_no_exact_ru = ahmm->exact_motif_count;
        vntr.exact_total_no_ru = ahmm->motif_count;
        vntr.exact_rl = ahmm->repeat_tract_len;
        vntr.exact_trf_score = ahmm->trf_score;

        if (debug)
        {
            std::cerr << "\n";
            std::cerr << "repeat_tract              : " << vntr.exact_repeat_tract << "\n";
            std::cerr << "position                  : [" << vntr.exact_beg1 << "," << vntr.exact_end1 << "]\n";
            std::cerr << "motif_concordance         : " << vntr.exact_motif_concordance << "\n";
            std::cerr << "repeat units              : " << vntr.exact_rl << "\n";
            std::cerr << "exact repeat units        : " << vntr.exact_no_exact_ru << "\n";
            std::cerr << "total no. of repeat units : " << vntr.exact_total_no_ru << "\n";
            std::cerr << "\n";
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
        vntr.fuzzy_motif_concordance = lfhmm->motif_concordance;
        vntr.fuzzy_no_exact_ru = lfhmm->exact_motif_count;
        vntr.fuzzy_total_no_ru = lfhmm->motif_count;
        vntr.fuzzy_rl = rflank_beg1-lflank_end1-1;
        
        if (lflank) free(lflank);
        if (rflank) free(rflank);

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
void FlankDetector::compute_purity_score(Variant& variant, std::string mode)
{
    std::string& repeat_tract = (mode=="e") ? variant.vntr.exact_repeat_tract : variant.vntr.fuzzy_repeat_tract;
    compute_purity_score(repeat_tract, variant.vntr.motif);
    
    variant.vntr.ru = ru;
    variant.vntr.exact_motif_concordance = motif_concordance;
    variant.vntr.exact_no_exact_ru = no_exact_ru;
    variant.vntr.exact_total_no_ru = total_no_ru;
    variant.vntr.exact_rl = rl;
}

/**
 * Computes purity score of a sequence with respect to a motif.
 */
void FlankDetector::compute_purity_score(std::string& repeat_tract, std::string motif)
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
            motif_concordance = 1;
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

    motif_concordance = ahmm->motif_concordance;
    no_exact_ru = ahmm->exact_motif_count;
    total_no_ru = ahmm->motif_count;
    rl = ahmm->repeat_tract_len;
    trf_score = ahmm->trf_score;

    return;
}