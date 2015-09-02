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

#include "vntr_annotator.h"

/**
 * Constructor.
 */
VNTRAnnotator::VNTRAnnotator(std::string& ref_fasta_file, bool debug)
{
    vm = new VariantManip(ref_fasta_file.c_str());



    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }

    cre = new CandidateRegionExtractor(ref_fasta_file, debug);

    max_mlen = 10;
    mt = new MotifTree(max_mlen, debug);

    this->debug = debug;
    qual.assign(256, 'K');
};

/**
 * Destructor.
 */
VNTRAnnotator::~VNTRAnnotator()
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
 * Annotates VNTR characteristics.
 * @mode -
 */
void VNTRAnnotator::annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant, std::string mode)
{
    VNTR& vntr = variant.vntr;

    //update chromosome and position
    variant.rid = bcf_get_rid(v);
    variant.pos1 = bcf_get_pos1(v);

    if (variant.type==VT_VNTR)
    {
        if (debug) std::cerr << "ANNOTATING VNTR/STR \n";

        //1. pick candidate region
        vntr.repeat_tract.assign(bcf_get_ref(v));
        vntr.rbeg1 = bcf_get_pos1(v);
 //       pick_candidate_region(h, v, vntr, REFERENCE);

        //2. detect candidate motifs from a reference seqeuence
        pick_candidate_motifs(h, v, vntr);

        //3. choose the best candidate motif
        choose_best_motif(h, v, mt, vntr, REFERENCE);

    }
    else if (variant.type&VT_INDEL)
    {
        //EXACT MODE
        if (mode=="e")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL EXACTLY\n";

            //1. pick candidate region using exact left and right alignment
            cre->extract_regions_by_exact_alignment(h, v, vntr);

            //2. detect candidate motifs from a reference sequence
            pick_candidate_motifs(h, v, vntr);

            //3. choose the best candidate motif
            choose_best_motif(h, v, mt, vntr, PICK_BEST_MOTIF);

            //4. evaluate reference length
            detect_repeat_region(h, v, variant, CLIP_ENDS);

            if (debug) std::cerr << "============================================\n";
            return;
        }
        //FUZZY DETECTION
        else if (mode=="f")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL FUZZILY\n";

            //1. selects candidate region by fuzzy left and right alignment
            cre->extract_regions_by_fuzzy_alignment(h, v, vntr);

            //2. detect candidate motifs from candidate region
            pick_candidate_motifs(h, v, vntr);

            //3. choose the best candidate motif
            choose_best_motif(h, v, mt, vntr, PICK_BEST_MOTIF);

            //4. evaluate reference length
            detect_repeat_region(h, v, variant, CLIP_ENDS);

            if (debug) std::cerr << "============================================\n";
            return;
        }
        else if (mode=="p")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL FUZZILY WITH PENALTY\n";

            //1. selects candidate region by fuzzy left and right alignment
            cre->extract_regions_by_fuzzy_alignment_with_penalty(h, v, vntr);

            //2. detect candidate motifs from candidate region
            pick_candidate_motifs(h, v, vntr);

            //3. choose the best candidate motif
            choose_best_motif(h, v, mt, vntr, PICK_BEST_MOTIF);

            //4. evaluate reference length
            detect_repeat_region(h, v, variant, CLIP_1L2R);

            if (debug) std::cerr << "============================================\n";
            return;
        }
        else if (mode=="h")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL USING raHMMs\n";

            //1. selects candidate region by fuzzy left and right alignment
            cre->extract_regions_by_fuzzy_alignment_with_penalty(h, v, vntr);

            //2. detect candidate motifs from candidate region
            pick_candidate_motifs(h, v, vntr);

            //3. choose the best candidate motif
            choose_best_motif(h, v, mt, vntr, PICK_BEST_MOTIF);

            //4. evaluate reference length
            detect_repeat_region(h, v, variant, CLIP_1L2R);

            if (debug) std::cerr << "============================================\n";
            return;

        }
        else if (mode=="x")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "Integrated Methods\n";

            //1. selects candidate region by fuzzy left and right alignment
            cre->extract_regions_by_fuzzy_alignment_with_penalty(h, v, vntr);

            //2. detect candidate motifs from candidate region
            pick_candidate_motifs(h, v, vntr);

            //3. choose the best candidate motif
            choose_best_motif(h, v, mt, vntr, PICK_BEST_MOTIF);

            //4. evaluate reference length
            detect_repeat_region(h, v, variant, CLIP_1L2R);

            if (debug) std::cerr << "============================================\n";
            return;
        }
    }
}

/**
 * Pick candidate motifs in different modes.
 * Invokes motif tree and the candidate motifs are stored in a
 * heap within the motif tree.
 */
void VNTRAnnotator::pick_candidate_motifs(bcf_hdr_t* h, bcf1_t* v, VNTR& vntr)
{
    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "PICK CANDIDATE MOTIFS\n\n";
    }

    mt->detect_candidate_motifs(vntr.repeat_tract);
}

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
void VNTRAnnotator::choose_best_motif(bcf_hdr_t* h, bcf1_t* v, MotifTree* mt, VNTR& vntr, uint32_t mode)
{
    if (debug)
    {
        std::cerr << "********************************************\n";
        std::cerr << "PICK BEST MOTIF\n\n";
    }

    if (mode==PICK_BEST_MOTIF)
    {
        if (!mt->pcm.empty())
        {
            CandidateMotif cm = mt->pcm.top();

            vntr.motif = cm.motif;
            vntr.motif_score = cm.score;
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
 * Detect repeat region.
 */
void VNTRAnnotator::detect_repeat_region(bcf_hdr_t* h, bcf1_t *v, Variant& variant, uint32_t mode)
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
//        vntr.rl = (float)vntr.repeat_tract.size()/vntr.ru.size();
        vntr.rl = (float)vntr.repeat_tract.size();

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

        if (debug)
        {
            vntr.print();
        }
    }

    //fill in flanks
    const char* chrom = variant.chrom.c_str();
    uint32_t pos1 = vntr.rbeg1;
    int32_t len = 0;
    faidx_fetch_seq(fai, chrom, pos1-10, pos1-1, &len);

};

/**
 * Chooses a phase of the motif that is appropriate for the alignment
 */
std::string VNTRAnnotator::choose_repeat_unit(std::string& ref, std::string& motif)
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

/**
 * Checks if a vntr is a homopolymer.
 */
bool VNTRAnnotator::is_homopolymer(bcf_hdr_t* h, bcf1_t* v)
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
 * This is a quick scan for a motif that is exactly repeated.
 */
std::string VNTRAnnotator::scan_exact_motif(std::string& sequence)
{
    size_t i = 0;
    size_t len = sequence.size();
    size_t sub_motif_len;
    const char* seq = sequence.c_str();
    size_t d = 0;
    size_t max_sub_motif_len = len/2;
    for (sub_motif_len = 1; sub_motif_len<=max_sub_motif_len; ++sub_motif_len)
    {
        size_t n_sub_motif = len/sub_motif_len;

        bool exact = true;
        size_t concordant = 0;
        size_t c1 = 0;

        for (size_t j=0; j<n_sub_motif; ++j)
        {
            if ((strncmp(&seq[0], &seq[j*sub_motif_len], sub_motif_len)))
            {
                exact = false;

                for (size_t k=0; k<sub_motif_len; ++k)
                {
                    if (seq[j*sub_motif_len+k]==seq[k])
                    {
                       ++concordant;
                    }
                }
            }
            else
            {
                ++c1;
                concordant += sub_motif_len;
            }
        }

        if (n_sub_motif*sub_motif_len<len)
        {
            if (strncmp(&seq[0], &seq[n_sub_motif*sub_motif_len], len-n_sub_motif*sub_motif_len))
            {
                exact = false;
            }
        }

        if (exact)
        {
            return sequence.substr(0, sub_motif_len);
        }
    }

    return "";
}

/**
 * Trim alleles.
 */
void VNTRAnnotator::trim(int32_t& pos1, std::string& ref, std::string& alt)
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
 * Detect allele lower bound extent.
 */
void VNTRAnnotator::detect_lower_bound_allele_extent(const char* chrom, int32_t& pos1, std::vector<std::string>& alleles, int32_t& start1, int32_t& end1)
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
 * Extracts the shortest repeat unit in a sequence.
 */
char* VNTRAnnotator::get_shortest_repeat_motif(char* allele, int32_t len)
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

/**
 * Gets motif of a repeat unit.
 */
std::string VNTRAnnotator::get_motif(std::string& ru)
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
std::string VNTRAnnotator::reverse_complement(std::string& seq)
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
std::string VNTRAnnotator::shift_phase(std::string& seq, size_t i)
{
    i = i<seq.size() ? i : i%seq.size();
    std::string shifted = seq.substr(i, seq.size()-i);
    shifted.append(seq, 0, i);

    return shifted;
}

/**
 * Returns true if is to be classified as an STR
 */
bool VNTRAnnotator::is_vntr(Variant& variant, int32_t mode)
{
    uint32_t mlen = variant.vntr.mlen;
    uint32_t rlen = variant.vntr.rend1 - variant.vntr.rbeg1 + 1;

    if (mode==TAN_KANG_2015_VNTR)
    {
        return (rlen - mlen >= 6);

//        equivalent to    
//        (mlen==1 && rlen>=5)  ||
//        (mlen==2 && rlen>=11) ||
//        (mlen==3 && rlen>=14) ||
//        (mlen==4 && rlen>=14) ||
//        (mlen==5 && rlen>=16) ||
//        (mlen==6 && rlen>=17) ||
//        (mlen>=7 && rlen>=mlen*2)
    
    }
    else if (mode==WILLEMS_2014_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=11) ||
                (mlen==3 && rlen>=14) ||
                (mlen==4 && rlen>=14) ||
                (mlen==5 && rlen>=16) ||
                (mlen==6 && rlen>=17) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else if (mode==ANANDA_2013_STR)
    {
        return ((mlen==1 && rlen>=2) ||
                (mlen==2 && rlen>=4) ||
                (mlen==3 && rlen>=6) ||
                (mlen==4 && rlen>=8) ||
                (mlen==5 && rlen>=10) ||
                (mlen==6 && rlen>=12) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else if (mode==FONDON_2012_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=13) ||
                (mlen==3 && rlen>=20) ||
                (mlen==4 && rlen>=23) ||
                (mlen==5 && rlen>=27) ||
                (mlen==6 && rlen>=27));
    }
    else if (mode==KELKAR_2008_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=10) ||
                (mlen==3 && rlen>=6) ||
                (mlen==4 && rlen>=8) ||
                (mlen==5 && rlen>=10) ||
                (mlen==6 && rlen>=12) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else if (mode==LAI_2003_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=8) ||
                (mlen==3 && rlen>=12) ||
                (mlen==4 && rlen>=16) ||
                (mlen==5 && rlen>=20) ||
                (mlen==6 && rlen>=24) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else
    {
        fprintf(stderr, "[%s:%d %s] STR definition mode does not exist: %d\n", __FILE__,__LINE__,__FUNCTION__, mode);
        exit(1);
    }
}
