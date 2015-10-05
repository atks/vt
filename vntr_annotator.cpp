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
    cmp = new CandidateMotifPicker(debug);
    fd = new FlankDetector(ref_fasta_file, debug);
    
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

    //this is for reannotating an VNTR record
    //this is more for the purpose of evaluation to
    //check if vt's algorithm is concordant with
    //VNTRs from other sources.
    if (variant.type==VT_VNTR)
    {
        if (debug) std::cerr << "ANNOTATING VNTR/STR \n";

        //1. pick candidate region
        vntr.repeat_tract.assign(bcf_get_ref(v));
        vntr.rbeg1 = bcf_get_pos1(v);
        char** alleles = bcf_get_allele(v);
        vntr.rend1 = strlen(alleles[0]);
        vntr.fuzzy_rbeg1 = vntr.rbeg1;
        vntr.fuzzy_rend1 = vntr.rend1;
        
        //2. detect candidate motifs from a reference seqeuence
        cmp->generate_candidate_motifs(h, v, variant);
        cmp->next_motif(h, v, variant);
    }
    //main purpose - annotation of Indels.
    else if (variant.type&VT_INDEL)
    {
        //the basic steps in annotating a TR
        //
        //1. extract a region that has a chance of containing the repeat units
        //2. choose a set of candidate motifs
        //3. choose the motif
        //4. detect repeat region and evaluate

        //EXACT MODE
        if (mode=="e")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL EXACTLY\n";

            //1. pick candidate region using exact left and right alignment
            cre->extract_regions_by_exact_alignment(h, v, vntr);

            //2. evaluate reference length
            fd->detect_flanks(h, v, variant, CLIP_ENDS);

            if (debug) std::cerr << "============================================\n";
            return;
        }
        //FUZZY DETECTION
        else if (mode=="f")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL FUZZILY\n";

            //1. selects candidate region by fuzzy left and right alignment
            cre->extract_regions_by_exact_alignment(h, v, vntr);

            //2. detect candidate motifs from a reference sequence
            cmp->generate_candidate_motifs(h, v, variant);
            cmp->next_motif(h, v, variant);
           
            //3. evaluate reference length
            fd->detect_flanks(h, v, variant, FRAHMM);

            //introduce reiteration!!!!

            if (debug) std::cerr << "============================================\n";
            return;
        }
    }
}

/**
 * Returns true if is to be classified as an STR
 */
bool VNTRAnnotator::is_vntr(Variant& variant, int32_t mode, std::string& method)
{
    uint32_t mlen = 0;
    uint32_t rlen = 0;
    float motif_concordance = 0;
    uint32_t no_exact_ru = 0;

    if (method == "e")
    {
        mlen = variant.vntr.mlen;
        rlen = variant.vntr.rl;
        motif_concordance = variant.vntr.motif_concordance;
        no_exact_ru = variant.vntr.no_exact_ru;
    }
    else if (method == "f")
    {
        mlen = variant.vntr.mlen;
        rlen = variant.vntr.fuzzy_rl;
        motif_concordance = variant.vntr.fuzzy_motif_concordance;
        no_exact_ru = variant.vntr.fuzzy_no_exact_ru;
    }

    if (mode==TAN_KANG_2015_VNTR)
    {
//        if (rlen - mlen < 6)
//        {
//            variant.vntr.print();
//        }
//        std::cerr << "rlen " << rlen << "\n";
//        std::cerr << "mlen " << mlen << "\n";
//        std::cerr << "no_exact_ru " << variant.vntr.no_exact_ru << "\n";
        if ((rlen - mlen) >= 6 && no_exact_ru>=2)
        {
            if (mlen==1 && motif_concordance>0.9)
            {
                return true;
            }
            else if (mlen>1 || motif_concordance>0.75)
            {
                return true;
            }
        }

        return false;
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
