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
    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }

    ///////////////////
    //initialize raHMMs
    ///////////////////
    float delta = 0.0001;
    float epsilon = 0.05;
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

    this->debug = debug;
    qual.assign(256, 'K');
};

/**
 * Destructor.
 */
FlankDetector::~FlankDetector()
{
    fai_destroy(fai);
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

        if (vntr.repeat_tract.size()>2)
        {
            vntr.repeat_tract = vntr.repeat_tract.substr(1, vntr.repeat_tract.size()-2);
            ++vntr.rbeg1;
        }

        vntr.ru = choose_repeat_unit(vntr.repeat_tract, vntr.motif);
        vntr.rend1 = vntr.rbeg1+vntr.rl-1;

        ahmm->set_model(vntr.ru.c_str());
        ahmm->align(vntr.repeat_tract.c_str(), qual.c_str());

        vntr.motif_concordance = ahmm->motif_concordance;
        vntr.no_exact_ru = ahmm->exact_motif_count;
        vntr.total_no_ru = ahmm->motif_count;
        vntr.rl = ahmm->repeat_tract_len;

//statistics for repeat unit
//    float motif_score;          //motif score from motif tree
//    float motif_concordance;    //motif concordance from hmm
//    float rl;                   //number of repeat units on repeat tract
//    float no_exact_ru;          //number exact repeat units from hmm
//    float total_no_ru;          //total no of repeat units from hmm

        if (debug)
        {
            vntr.print();
        }
    }
    //simple single base pair clipping of ends
    else if (mode==CLIP_1L2R)
    {
        if (debug)
        {
            std::cerr << "********************************************\n";
            std::cerr << "CLIP ENDS\n";
            std:: cerr << "\n";
        }

        if (vntr.repeat_tract.size()>3)
        {
            vntr.repeat_tract = vntr.repeat_tract.substr(1, vntr.repeat_tract.size()-3);
            ++vntr.rbeg1;
        }

        vntr.ru = choose_repeat_unit(vntr.repeat_tract, vntr.motif);
        vntr.rl = (float)vntr.repeat_tract.size()/vntr.ru.size();
        vntr.rend1 = vntr.rbeg1 +  vntr.rl -1;

        if (debug)
        {
            vntr.print();
        }
    }
    else if (mode==FRAHMM)
    {
        if (debug)
        {
            std::cerr << "********************************************\n";
            std::cerr << "4. DETECTING REPEAT TRACT FUZZILY\n";
        }

        /////////////////
        //exact alignment
        /////////////////
        if (debug)
        {
            std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
            std::cerr << "4a. Exact left/right alignment\n";
        }
        
        if (vntr.repeat_tract.size()>2)
        {
            if (vntr.mlen==1)
            {
                int32_t offset = 0;
                int32_t length = vntr.repeat_tract.size();
                if (vntr.repeat_tract.at(0)!=vntr.motif.at(0))
                {
                    offset = 1;  
                    ++vntr.rbeg1;  
                }
                
                if (vntr.repeat_tract.at(vntr.repeat_tract.size()-1)!=vntr.motif.at(0))
                {
                    length -= offset+1;    
                    --vntr.rend1;
                }
                
                vntr.repeat_tract = vntr.repeat_tract.substr(offset, length);
            }
            else
            {
                
            }
        }
        
        vntr.ru = choose_repeat_unit(vntr.repeat_tract, vntr.motif);
        ahmm->set_model(vntr.ru.c_str());
        ahmm->align(vntr.repeat_tract.c_str(), qual.c_str());

        vntr.motif_concordance = ahmm->motif_concordance;
        vntr.no_exact_ru = ahmm->exact_motif_count;
        vntr.total_no_ru = ahmm->motif_count;
        vntr.rl = ahmm->repeat_tract_len;
        
        if (debug)
        {
            std::cerr << "\n";
            std::cerr << "repeat_tract              : " << vntr.repeat_tract << "\n";
            std::cerr << "position                  : [" << vntr.rbeg1 << "," << vntr.rend1 << "]\n";
            std::cerr << "motif_concordance         : " << vntr.motif_concordance << "\n";
            std::cerr << "repeat units              : " << vntr.rl << "\n";
            std::cerr << "exact repeat units        : " << vntr.no_exact_ru << "\n";
            std::cerr << "total no. of repeat units : " << vntr.total_no_ru << "\n";
            std::cerr << "\n";
        }
        
        ///////////////////////
        //fuzzy right alignment
        ///////////////////////
        if (debug)
        {
            std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
            std::cerr << "4a. Fuzzy right alignment\n";
        }
        int32_t slen = 100;
        
        //pick 5 bases to the right
        int32_t rflank_len;
        char* rflank = faidx_fetch_seq(fai, variant.chrom.c_str(), vntr.rbeg1-1, vntr.rbeg1+5-1-1, &rflank_len);

        //pick 105 bases for aligning
        int32_t seq_len;
        char* seq = faidx_fetch_seq(fai, variant.chrom.c_str(), vntr.rbeg1-slen-1, vntr.rbeg1+5-1-1, &seq_len);

        rfhmm->set_model(vntr.ru.c_str(), rflank);
        rfhmm->align(seq, qual.c_str());
        if (debug) rfhmm->print_alignment();

        if (rflank_len) free(rflank);
        if (seq_len) free(seq);

        //////////////////////
        //fuzzy left alignment
        //////////////////////
        if (debug)
        {
            std::cerr << "\n";
            std::cerr << "++++++++++++++++++++++++++++++++++++++++++++\n";
            std::cerr << "4b. Fuzzy left alignment\n";
        }
        int32_t lflank_end1;
        int32_t rflank_beg1;

        lflank_end1 = vntr.rbeg1-slen-1+1 + rfhmm->get_lflank_read_epos1() - 1;

        slen = 100;

        //pick 5 bases to right
        int32_t lflank_len;
        char* lflank;
        while(true)
        {
            lflank = faidx_fetch_seq(fai, variant.chrom.c_str(), lflank_end1-5-1, lflank_end1-1, &lflank_len);

            //pick 105 bases for aligning
            seq = faidx_fetch_seq(fai, variant.chrom.c_str(), lflank_end1-5-1, lflank_end1+slen-1-1, &seq_len);

            lfhmm->set_model(lflank, vntr.ru.c_str());
            lfhmm->align(seq, qual.c_str());
            if (debug) lfhmm->print_alignment();

            if (seq_len) free(seq);

            if (lfhmm->get_rflank_read_epos1()!=INT32_MAX)
            {
                rflank_beg1 = lflank_end1 - 5 + lfhmm->get_rflank_read_spos1() -1;
                break;
            }
            else if (slen==1000)
            {
                rflank_beg1 = lflank_end1 - 5 + 1000 -1;
                vntr.is_large_repeat_tract = true;
                break;
            }
            else
            {
                slen +=100;
            }
        }

        if (lflank_len) free(lflank);

        lflank = faidx_fetch_seq(fai, variant.chrom.c_str(), lflank_end1-10-1, lflank_end1-1, &lflank_len);
        rflank = faidx_fetch_seq(fai, variant.chrom.c_str(), rflank_beg1-1, rflank_beg1 +10 -1 -1, &rflank_len);

        vntr.fuzzy_rbeg1 = lflank_end1+1;
        vntr.fuzzy_rend1 = rflank_beg1-1;
        int32_t repeat_tract_len;
        char* repeat_tract = faidx_fetch_seq(fai, variant.chrom.c_str(), lflank_end1, rflank_beg1-1-1, &repeat_tract_len);
        vntr.fuzzy_repeat_tract.assign(repeat_tract);
        if (repeat_tract_len) free(repeat_tract);
        vntr.fuzzy_motif_concordance = lfhmm->motif_concordance;
        vntr.fuzzy_no_exact_ru = lfhmm->exact_motif_count;
        vntr.fuzzy_total_no_ru = lfhmm->motif_count;
        vntr.fuzzy_rl = rflank_beg1-lflank_end1-1;


//        std::cerr << "lflank_end1 : lflank "<< lflank_end1 << ":" << lflank << "\n";
//        std::cerr << "rflank_beg1 : rflank " << rflank_beg1 << ":" << rflank << "\n";


        if (lflank_len) free(lflank);
        if (rflank_len) free(rflank);

        if (debug)
        {
            std::cerr << "\n";
            vntr.print();
            std::cerr << "\n";
        }
    }

    //fill in flanks
    const char* chrom = variant.chrom.c_str();
    uint32_t pos1 = vntr.rbeg1;
    int32_t len = 0;
    faidx_fetch_seq(fai, chrom, pos1-10, pos1-1, &len);

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
 * Chooses a phase of the motif that is appropriate for the alignment
 */
std::string FlankDetector::choose_repeat_unit(std::string& ref, std::string& motif)
{
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = shift_str(motif, i);
        if (ref.compare(0, smotif.size(), smotif)==0)
        {
            return smotif;
        }
    }

    return motif;
}

/**
 * Checks if a vntr is a homopolymer.
 */
bool FlankDetector::is_homopolymer(bcf_hdr_t* h, bcf1_t* v)
{
    bool is_homopolymer = false;
    uint32_t ref_len = strlen(bcf_get_ref(v));
    for (size_t i=1; i<bcf_get_n_allele(v); ++i)
    {
        std::string ref(bcf_get_alt(v, 0));
        std::string alt(bcf_get_alt(v, i));
        int32_t pos1 = bcf_get_pos1(v);
    }

    return is_homopolymer;
}

/**
 * Trim alleles.
 */
void FlankDetector::trim(int32_t& pos1, std::string& ref, std::string& alt)
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
 * Gets motif of a repeat unit.
 */
std::string FlankDetector::get_motif(std::string& ru)
{
    std::string motif = "";
    for (size_t i=0; i<ru.size(); ++i)
    {
        std::string phase = shift_phase(ru, i);
        std::string rc = reverse_complement(phase);
        motif = phase < rc ? phase : rc;
    }

    return motif;
}

/**
 * Reverse complement a sequence.
 */
std::string FlankDetector::reverse_complement(std::string& seq)
{
    std::string rc = "";

    for (size_t i=seq.size()-1; i>0; --i)
    {
        char b = seq.at(i);

        switch (b)
        {
            case 'A':
                rc.append(1, 'A');
                break;
            case 'C':
                rc.append(1, 'C');
                break;
            case 'G':
                rc.append(1, 'G');
                break;
            case 'T':
                rc.append(1, 'T');
                break;
        }
    }

    return rc;
}

/**
 * Shifts a sequence to the right by i bases.
 */
std::string FlankDetector::shift_phase(std::string& seq, size_t i)
{
    i = i<seq.size() ? i : i%seq.size();
    std::string shifted = seq.substr(i, seq.size()-i);
    shifted.append(seq, 0, i);

    return shifted;
}