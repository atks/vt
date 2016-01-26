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
    cmp = new CandidateMotifPicker(debug);
    fd = new FlankDetector(ref_fasta_file, debug);

    this->debug = debug;
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
 * @mode - e for exact alignment annotation
 *       - f for fuzzy alignment annotation
 *       -
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
        if (mode=="r")
        {
            if (debug) std::cerr << "ANNOTATING VNTR/STR \n";
            //takes accepted MOTIF and RU

            //computes purity
            //1. pick candidate region
            cre->pick_candidate_region(h, v, variant, REFERENCE);

            //just use the motif annotator
            //2. detect candidate motifs from a reference seqeuence
            cmp->generate_candidate_motifs(h, v, variant);
            cmp->next_motif(h, v, variant);

            //3. compute purity
            fd->detect_flanks(h, v, variant, CLIP_ENDS);
        }
        else if (mode=="c")
        {
            if (debug) std::cerr << "ANNOTATING VNTR's purity score\n";

            //1. pick candidate region
            cre->pick_candidate_region(h, v, variant, REFERENCE);

            //2. set motifs from info field
            cmp->set_motif_from_info_field(variant);

            //3. compute purity scores
            fd->compute_purity_score(variant, "e");
        }
    }
    //main purpose - annotation of Indels.
    else if (variant.type&VT_INDEL)
    {
        //the basic steps in annotating a TR
        //
        //1. extract a region that has a chance of containing the repeat units
        //2. choose a set of candidate motifs and pick motif
        //3. detect repeat region and evaluate
        //4. iterate 2 and 3

        //EXACT MODE
        if (mode=="e")
        {
            if (debug) std::cerr << "============================================\n";
            if (debug) std::cerr << "ANNOTATING INDEL EXACTLY\n";

            //1. pick candidate region using exact left and right alignment
            cre->pick_candidate_region(h, v, variant, EXACT_LEFT_RIGHT_ALIGNMENT);

            //2. detect candidate motifs from a reference sequence
            cmp->generate_candidate_motifs(h, v, variant);

            if (!cmp->next_motif(h, v, variant))
            {
                std::cerr << "oops, no candidate motif for next step\n";
            }

            //3. detect flanks and evaluate reference tract
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
            cre->pick_candidate_region(h, v, variant, EXACT_LEFT_RIGHT_ALIGNMENT);

            //2. detect candidate motifs from a reference sequence
            cmp->generate_candidate_motifs(h, v, variant);

            if (!cmp->next_motif(h, v, variant))
            {
                std::cerr << "oops, no candidate motif for next step\n";
            }

            //3. evaluate reference length
            fd->detect_flanks(h, v, variant, FRAHMM);

            //introduce reiteration based on concordance and exact concordance.

            if (debug) std::cerr << "============================================\n";
            return;
        }
    }
}

/**
 * Return the cannonicalized representation of a motif.
 */
std::string VNTRAnnotator::cannonicalize(std::string& motif)
{
    std::string cmotif = motif;

    for (uint32_t i=0; i<motif.size(); ++i)
    {
        std::string shifted_motif = motif.substr(i) + motif.substr(0,i);

        if (shifted_motif > cmotif)
        {
            cmotif = shifted_motif;
        }
    }

    return cmotif;
}

/**
 * Return the string of unique bases in a motif.
 */
std::string VNTRAnnotator::basis(std::string& motif)
{
    bool bases[4];

    for (uint32_t i=0; i<motif.size(); ++i)
    {
        char base = motif.at(i);

        if (base<=67)
        {
            if (base==65)
            {
                if (!bases[0])
                {
                    bases[0] = true;
                }
            }
            else
            {
                if (!bases[1])
                {
                    bases[1] = true;
                }
            }
        }
        else
        {
            if (base==71)
            {
                if (!bases[2])
                {
                    bases[2] = true;
                }
            }
            else
            {
                if (!bases[3])
                {
                    bases[3] = true;
                }
            }
        }
    }

    std::string basis;
    if (bases[0]) basis.append(1, 'A');
    if (bases[1]) basis.append(1, 'C');
    if (bases[2]) basis.append(1, 'G');
    if (bases[3]) basis.append(1, 'T');

    return basis;
}

/**
 * Returns true if is to be classified as an STR
 */
bool VNTRAnnotator::is_vntr(Variant& variant, int32_t mode, std::string& method)
{
    uint32_t mlen = variant.vntr.mlen;
    uint32_t rlen = variant.vntr.exact_rl;
    float motif_concordance = 0;
    uint32_t no_exact_ru = 0;


    VNTR& vntr = variant.vntr;

    if (method == "e")
    {
        motif_concordance = variant.vntr.exact_motif_concordance;
        no_exact_ru = variant.vntr.exact_no_exact_ru;
    }
    else if (method == "f")
    {
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
        if (method == "f")
        {
            if ((vntr.fuzzy_rl - vntr.mlen)>=6 && vntr.fuzzy_no_exact_ru>=2)
            {
                //should we set separate cutoffs for motifs that contain only 2 types of nucleotides?
                //
                if (vntr.mlen==1 && vntr.fuzzy_motif_concordance>0.9)
                {
                    vntr.definition_support = "f";
                    return true;
                }
                else if (vntr.mlen>1 && vntr.fuzzy_motif_concordance>0.75)
                {
                    vntr.definition_support = "f";
                    return true;
                }
            }
        }

        //if the fuzzy definition is not caught above, this will fall through to the exact definition.

        if ((vntr.exact_rl - vntr.mlen)>=6 && vntr.exact_no_exact_ru>=2)
        {
            //should we set separate cutoffs for motifs that contain only 2 types of nucleotides?
            //
            if (vntr.mlen==1 && vntr.exact_motif_concordance>0.9)
            {
                vntr.definition_support = "e";
                return true;
            }
            else if (vntr.mlen>1 && vntr.exact_motif_concordance>0.75)
            {
                vntr.definition_support = "e";
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
