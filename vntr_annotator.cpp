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
void VNTRAnnotator::annotate(Variant& variant, int32_t amode)
{
    VNTR& vntr = variant.vntr;

    bcf_hdr_t* h = variant.h;
    bcf1_t* v = variant.v;

    //update chromosome and position
    variant.rid = bcf_get_rid(v);
    variant.pos1 = bcf_get_pos1(v);

    //this is for reannotating an VNTR record
    //this is more for the purpose of evaluation to
    //check if vt's algorithm is concordant with
    //VNTRs from other sources.
    if (variant.type==VT_VNTR)
    {
        if (debug) std::cerr << "============================================\n";
        if (debug) std::cerr << "ANNOTATING VNTR/STR \n";
        
        //1. pick candidate region
        cre->pick_candidate_region(variant, REFERENCE, FINAL);

        //2. set motifs from info field
        cmp->set_motif_from_info_field(variant);

        //3. compute purity scores
        fd->compute_purity_score(variant, FINAL);

        //4. compute composition and sequence statistics
        fd->compute_composition_and_entropy(variant, FINAL);
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
        
        if (debug) std::cerr << "============================================\n";
        if (debug) std::cerr << "ANNOTATING INDEL\n";

        //1. selects candidate region by fuzzy left and right alignment
        cre->pick_candidate_region(variant, EXACT_LEFT_RIGHT_ALIGNMENT, EXACT);
        cmp->update_exact_repeat_unit(variant);
        
        //2. detect candidate motifs from a reference sequence
        cmp->generate_candidate_motifs(variant);

        //this cannot possibly fail as next_motif() guarantees it
        if (!cmp->next_motif(variant, CHECK_MOTIF_PRESENCE_IN_ALLELE))
        {
            //fall back on exact motif chosen.
            VNTR& vntr = variant.vntr;
            vntr.fuzzy_motif = vntr.exact_motif;
            vntr.fuzzy_ru = vntr.exact_ru;
            vntr.fuzzy_basis = vntr.exact_basis;
            vntr.fuzzy_mlen = vntr.exact_mlen;
            vntr.fuzzy_blen = vntr.exact_blen;;
            
            if (debug) std::cerr << "updating fuzzy motif with exact motifs\n";
        }

        fd->detect_flanks(variant, EXACT);
        fd->compute_purity_score(variant, EXACT);
        fd->compute_composition_and_entropy(variant, EXACT);

        fd->detect_flanks(variant, FUZZY);
        fd->compute_purity_score(variant, FUZZY);
        fd->compute_composition_and_entropy(variant, FUZZY);

        //introduce reiteration based on concordance and exact concordance.

        if (debug)
        {
            std::cerr << "\n";
            vntr.print();
            std::cerr << "\n";
        }

        if (debug) std::cerr << "============================================\n";
        return;
        
    }
}