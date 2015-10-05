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

#include "candidate_motif_picker.h"

/**
 * Constructor.
 */
CandidateMotifPicker::CandidateMotifPicker(std::string& ref_fasta_file, bool debug)
{
    vm = new VariantManip(ref_fasta_file.c_str());

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

    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }

    cre = new CandidateRegionExtractor(ref_fasta_file, debug);
    fd = new FlankDetector(ref_fasta_file, debug);

    max_mlen = 10;
    mt = new MotifTree(max_mlen, debug);

    this->debug = debug;
    qual.assign(256, 'K');
};

/**
 * Destructor.
 */
CandidateMotifPicker::~CandidateMotifPicker()
{
    delete vm;
    fai_destroy(fai);

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
 * Pick candidate motifs in different modes.
 * Invokes motif tree and the candidate motifs are stored in a
 * heap within the motif tree.
 */
void CandidateMotifPicker::pick_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, Variant& variant)
{
    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "PICK CANDIDATE MOTIFS\n\n";
    }

    if (variant.ins && variant.alleles.size()==1)
    {   
        char** alleles = bcf_get_allele(v);
        
        if (debug)
        {     
            const char* repeat_tract = variant.vntr.repeat_tract.c_str();
            std::cerr << "Longest Allele : "   << alleles[0][0] << "[" <<  &alleles[1][1]  << "]" << &repeat_tract[1] << "\n"; 
        }
        
        //spike in inserted allele
        std::string spiked_seq(alleles[1]);
        std::string insertion = variant.vntr.repeat_tract.substr(strlen(alleles[0]), variant.vntr.repeat_tract.size()-strlen(alleles[0]));
        spiked_seq.append(insertion); 
           
        
        mt->detect_candidate_motifs(spiked_seq);   
    }
    else
    {
        mt->detect_candidate_motifs(variant.vntr.repeat_tract);
    }
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
void CandidateMotifPicker::choose_best_motif(bcf_hdr_t* h, bcf1_t* v, MotifTree* mt, VNTR& vntr, uint32_t mode)
{
    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "PICK BEST MOTIF\n\n";
    }

    if (mode==PICK_BEST_MOTIF)
    {
        //extract deleted or inserted fragments
        //this should work for multiallelics
//        std::vector<std::string> indels;
//        char* ref = bcf_get_ref_allele(v, 0);
//                
//        
//        char* ref = bcf_get_ref_allele(v, 0);
//        
        
        
        
        if (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            vntr.motif = cm.motif;
            vntr.mlen = cm.motif.size();
            vntr.motif_score = cm.score;

            //check for consistency between chosen motif and inserted/deleted fragment

        }

        if (debug)
        {
            printf("selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", mt->pcm.top().motif.c_str(),
                                                            mt->pcm.top().score,
                                                            mt->pcm.top().fit,
                                                            ahmm->get_motif_concordance(),
                                                            (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                            ahmm->get_exact_motif_count(),
                                                            ahmm->get_motif_count());
        }
    }
    else if (mode==11)
    {
        if (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            vntr.motif = cm.motif;
            vntr.mlen = cm.motif.size();
            vntr.motif_score = cm.score;

//            if (cm.motif.size()==1)
//            {    
//                if (cm.score>0.8)
//                {
//                    vntr.motif = cm.motif;
//                    vntr.mlen = cm.motif.size();
//                    vntr.motif_score = cm.score;
//                }
//                else
//                {
//                    while ()
//                    {
//                        mt->pcm.pop();
//                
//                    cm = mt->pcm.top();
//                    vntr.motif = cm.motif;
//                    vntr.mlen = cm.motif.size();
//                    vntr.motif_score = cm.score;
//                }
//            }
//            else if (cm.motif.size()>1)
//            {
//                if (cm.score>0.7)
//                {    
//                    vntr.motif = cm.motif;
//                    vntr.mlen = cm.motif.size();
//                    vntr.motif_score = cm.score;
//                }
//                else
//                {
//                    mt->pcm.pop();
//                
//                    cm = mt->pcm.top();
//                    vntr.motif = cm.motif;
//                    vntr.mlen = cm.motif.size();
//                    vntr.motif_score = cm.score; 
//                }
//            }
//            else
//            {
//                vntr.motif = cm.motif;
//                vntr.mlen = cm.motif.size();
//                vntr.motif_score = cm.score;
//            }
        }

        if (debug)
        {
            printf("selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", mt->pcm.top().motif.c_str(),
                                                            mt->pcm.top().score,
                                                            mt->pcm.top().fit,
                                                            ahmm->get_motif_concordance(),
                                                            (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                            ahmm->get_exact_motif_count(),
                                                            ahmm->get_motif_count());
        }
    }
    else if (mode==10) //backup plan to make sense of.
    {
        //choose candidate motif
        bool first = true;
        float cp = 0;
        float mscore = 0;
        uint32_t clen = 0;
        std::string ref(bcf_get_ref(v));
        if (ref.size()>256) ref = ref.substr(0,256);

        uint32_t no = 0;

        while (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            if (first)
            {
                std::string ru = choose_repeat_unit(ref, cm.motif);
                ahmm->set_model(ru.c_str());
                ahmm->align(ref.c_str(), qual.c_str());

                vntr.motif = cm.motif;
                vntr.mlen = cm.motif.size();

                vntr.motif_score = ahmm->get_motif_concordance();;
                vntr.ru = ru;

                cp = cm.score;
                clen = cm.len;
                first = false;
                mscore = ahmm->get_motif_concordance();

                if (debug)
                {
                    ahmm->print_alignment();

                    printf("selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n",
                            cm.motif.c_str(),
                            cm.score,
                            cm.fit,
                            ahmm->get_motif_concordance(),
                            (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                            ahmm->get_exact_motif_count(),
                            ahmm->get_motif_count());
                }
            }
            else
            {
                if (no<3 && (cp-cm.score<((float)1/cm.len)*cp || cm.score>0.4))
                {
                    //if score are not perfect, use AHMM
                    std::string ru = choose_repeat_unit(ref, cm.motif);
                    ahmm->set_model(ru.c_str());
                    ahmm->align(ref.c_str(), qual.c_str());

                    cp = cm.score;
                    clen = cm.len;

                    if (ahmm->get_motif_concordance()>mscore)
                    {
                        vntr.motif = cm.motif;
                        vntr.motif_score = ahmm->get_motif_concordance();;
                        vntr.ru = ru;
                    }

                    if (debug)
                    {
                        //ahmm->print_alignment();

                        printf("selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score,
                                                                        mt->pcm.top().fit,
                                                                        ahmm->get_motif_concordance(),
                                                                        (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                                        ahmm->get_exact_motif_count(),
                                                                        ahmm->get_motif_count());
                    }
                }
                else
                {
                    if (debug)
                    {
                        printf("not selected: %10s %.2f \n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score);
                    }
                }
            }

            if (!mt->pcm.empty()) mt->pcm.pop();

            ++no;
        }
    }
    else if (mode==ALLELE_FUZZY)
    {
        //choose candidate motif
        bool first = true;
        float cp = 0;
        float mscore = 0;
        uint32_t clen = 0;
        std::string ref(bcf_get_ref(v));
        if (ref.size()>256) ref = ref.substr(0,256);

        uint32_t no = 0;

        while (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            if (first)
            {
                std::string ru = choose_repeat_unit(ref, cm.motif);
                ahmm->set_model(ru.c_str());
                ahmm->align(ref.c_str(), qual.c_str());

                vntr.motif = cm.motif;
                vntr.motif_score = ahmm->get_motif_concordance();;
                vntr.ru = ru;

                cp = cm.score;
                clen = cm.len;
                first = false;
                mscore = ahmm->get_motif_concordance();

                if (debug)
                {
                    ahmm->print_alignment();

                    printf("    selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", cm.motif.c_str(),
                                                                    cm.score,
                                                                    cm.fit,
                                                                    ahmm->get_motif_concordance(),
                                                                    (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                                    ahmm->get_exact_motif_count(),
                                                                    ahmm->get_motif_count());
                }
            }
            else
            {
                if (no<3 && (cp-cm.score<((float)1/cm.len)*cp || cm.score>0.4))
                {
                    //if score are not perfect, use AHMM
                    std::string ru = choose_repeat_unit(ref, cm.motif);
                    ahmm->set_model(ru.c_str());
                    ahmm->align(ref.c_str(), qual.c_str());

                    cp = cm.score;
                    clen = cm.len;

                    if (ahmm->get_motif_concordance()>mscore)
                    {
                        vntr.motif = cm.motif;
                        vntr.motif_score = ahmm->get_motif_concordance();;
                        vntr.ru = ru;
                    }

                    if (debug)
                    {
                        ahmm->print_alignment();

                        printf("    selected: %10s %.2f %.2f %.2f %.2f (%d/%d)\n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score,
                                                                        mt->pcm.top().fit,
                                                                        ahmm->get_motif_concordance(),
                                                                        (float)ahmm->get_exact_motif_count()/ahmm->get_motif_count(),
                                                                        ahmm->get_exact_motif_count(),
                                                                        ahmm->get_motif_count());
                    }
                }
                else
                {
                    if (debug)
                    {
                        printf("not selected: %10s %.2f \n", mt->pcm.top().motif.c_str(),
                                                                        mt->pcm.top().score);
                    }
                }
            }

            if (!mt->pcm.empty()) mt->pcm.pop();

            ++no;
        }
    }
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
std::string CandidateMotifPicker::choose_repeat_unit(std::string& ref, std::string& motif)
{
    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string smotif = mt->shift_str(motif, i);
        if (ref.compare(0, smotif.size(), smotif)==0)
        {
            return smotif;
        }
    }

    return motif;
}

