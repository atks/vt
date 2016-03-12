/* The MIT License

   Copyright (c) 2016 Adrian Tan <atks@umich.edu>

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

#include "multiallelics_consolidator.h"

MultiallelicsConsolidator::MultiallelicsConsolidator(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& fexp, std::string& ref_fasta_file)
{
    //////////////////////
    //i/o initialization//
    //////////////////////
    buffer_window_allowance = 5000;

    this->input_vcf_file = input_vcf_file;
    this->output_vcf_file = output_vcf_file;
    odr = new BCFOrderedReader(input_vcf_file, intervals);
    odw = new BCFOrderedWriter(output_vcf_file, 2*buffer_window_allowance);
    odw->link_hdr(odr->hdr);

    //for adding empty genotype fields for a VCF file with individual information
    no_samples = bcf_hdr_nsamples(odw->hdr);
    gts = NULL;
    if (no_samples)
    {
        gts = (int32_t*) malloc(no_samples*sizeof(int32_t));
        for (uint32_t i=0; i<no_samples; ++i)
        {
            gts[i] = 0;
        }
    }

    ////////////////////////////////////////////
    //add relevant field for adding VNTR records
    ////////////////////////////////////////////
    bool rename = false;
    ASSOCIATED_BIALLEIC_VARIANTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "ASSOCIATED_BIALLEIC_VARIANTS", ".", "String", "Indels that induced this VNTR.", rename);
    odw->write_hdr();

    /////////////////////////
    //filter initialization//
    /////////////////////////
    filter.parse(fexp.c_str(), false);
    filter_exists = fexp=="" ? false : true;

    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_variants = 0;
    no_added_multiallelics = 0;

    ////////////////////////
    //tools initialization//
    ////////////////////////
    rs = new ReferenceSequence(ref_fasta_file);
}

/**
 * Inserts a biallelic record that is already present.
 * Also inserts a multiallelic if it induces it.
 *
 * assumption: variants are ordered.
 *
 * Note that this function does not guaranteed order of variants,
 * the positions of the multiallelic variants can change but its position
 * in the buffer is not updated, the variants when flushed, will be
 * ordered via OrderedBCFReader.
 */
void MultiallelicsConsolidator::insert(Variant* var)
{
    flush(var);

    Variant& nvar = *var;
    std::list<Variant*>::iterator i =vbuffer.begin();
    bool inserted = false;
    Variant* mvar = NULL;

    while(i!=vbuffer.end())
    {
        Variant& cvar = **i;

        if (nvar.rid > cvar.rid)
        {
            break;
        }
        else if (nvar.rid == cvar.rid)
        {
            if (nvar.beg1 <= cvar.beg1 + buffer_window_allowance)
            {
                if (nvar.beg1 <= cvar.end1 && nvar.end1>=cvar.beg1)
                {
                    mvar = create_or_update_multiallelic(nvar, cvar);
                    break;
                }
                else
                {
                    ++i;
                }
            }
            else
            {
                break;
            }
        }
    }

    vbuffer.push_front(var);
    if (mvar) vbuffer.push_front(mvar);
}

/**
 * Flush variant buffer.
 */
void MultiallelicsConsolidator::flush(Variant* var)
{
//    std::cerr << "enter flush\n";

    if (vbuffer.empty())
    {
        return;
    }

    if (var)
    {
        Variant& nvar = *var;

        //search for point to start deleting from.
        std::list<Variant*>::iterator i = vbuffer.begin();
        while(i!=vbuffer.end())
        {
            Variant& cvar = **i;

            if (nvar.rid > cvar.rid)
            {
                break;
            }
            else if (nvar.rid == cvar.rid)
            {
                if ((cvar.end1+buffer_window_allowance) < nvar.beg1)
                {
                    break;
                }
            }

            ++i;
        }

        //delete all records beyond this point
        while (i!=vbuffer.end())
        {
            process_exit(*i);
            i = vbuffer.erase(i);
        }
    }
    else
    {
        std::list<Variant*>::iterator i = vbuffer.begin();
        while (i!=vbuffer.end())
        {
            var = *i;
            process_exit(*i);
            i = vbuffer.erase(i);
        }
    }
}

/**
 * Flush variant buffer.
 */
void MultiallelicsConsolidator::process()
{
    bcf_hdr_t *h = odr->hdr;
    bcf1_t *v = odw->get_bcf1_from_pool();

    while (odr->read(v))
    {
        Variant* var = new Variant(h, v);

        if (filter_exists)
        {
            if (!filter.apply(h, v, var, false))
            {
                delete(var);
                continue;
            }
        }

//        bcf_print(h, v);

        insert(var);

        ++no_variants;

        v = odw->get_bcf1_from_pool();
    }

    flush();
    close();
};

/**
 * Process exiting variant.run
 */
void MultiallelicsConsolidator::process_exit(Variant* var)
{
    if (var->is_new_multiallelic)
    {
        update_multiallelic_for_printing(*var);
        ++no_added_multiallelics;
    }

//    bcf_print(var->h, var->v);
    odw->write(var->v);
}

/**.
 * Close files.
 */
void MultiallelicsConsolidator::close()
{
    odw->close();
    odr->close();
}

/**
 * Creates or updates a multiallelic.
 * @nvar - always a biallelic
 * @cvar - biallelic or multiallelic
 */
Variant* MultiallelicsConsolidator::create_or_update_multiallelic(Variant& nvar, Variant& cvar)
{
    //create new multiallelic
    if (!cvar.is_new_multiallelic && !cvar.is_involved_in_a_multiallelic)
    {
        Variant *mvar = new Variant();

        mvar->h = nvar.h;
        mvar->v = bcf_init();
        bcf_clear(mvar->v);
        bcf_set_rid(mvar->v, nvar.rid);
        bcf_set_pos1(mvar->v, std::min(nvar.beg1, cvar.beg1));

        std::string s = "NNN,<MULTI>";
        bcf_update_alleles_str(mvar->h, mvar->v, s.c_str());

        mvar->chrom = nvar.chrom;
        mvar->rid = nvar.rid;
        mvar->beg1 = std::min(nvar.beg1, cvar.beg1);
        bcf_set_pos1(mvar->v, mvar->beg1);
        mvar->end1 = std::max(nvar.end1, cvar.end1);
        mvar->is_new_multiallelic = true;
        mvar->is_involved_in_a_multiallelic = false;

        nvar.is_new_multiallelic = false;
        nvar.is_involved_in_a_multiallelic = true;
        nvar.associated_new_multiallelic = mvar;
        cvar.is_new_multiallelic = false;
        cvar.is_involved_in_a_multiallelic = true;
        cvar.associated_new_multiallelic = mvar;

        bcf1_t *nv = bcf_dup(nvar.v);
        bcf1_t *cv = bcf_dup(cvar.v);

        mvar->vs.push_back(nv);
        mvar->vs.push_back(cv);

//        std::cerr << "=========================\n";
//        std::cerr << "NEW MULTIALLELIC\n";
//        nvar.print();
//        cvar.print();
//        std::cerr << "mvar.beg1 : " << mvar->beg1 << "\n";
//        std::cerr << "mvar.end1 : " << mvar->end1 << "\n";
//        std::cerr << "nvar.beg1 : " << nvar.beg1 << "\n";
//        std::cerr << "nvar.end1 : " << nvar.end1 << "\n";
//        std::cerr << "cvar.beg1 : " << cvar.beg1 << "\n";
//        std::cerr << "cvar.end1 : " << cvar.end1 << "\n";
//        bcf_print(mvar->h, mvar->v);
//        bcf_print(mvar->h, nvar.v);
//        bcf_print(mvar->h, cvar.v);
//        std::cerr << "no. overlapping variants : " << mvar->vs.size() << "\n";
//        std::cerr << "=========================\n";

        return mvar;
    }

    //fall through - update existing multiallelic

    Variant *mvar = NULL;
    if (cvar.is_new_multiallelic)
    {
        mvar = &cvar;
    }
    else
    {
        mvar = cvar.associated_new_multiallelic;
    }

//    bcf_print(mvar->h, mvar->v);

    mvar->beg1 = std::min(nvar.beg1, mvar->beg1);
    bcf_set_pos1(mvar->v, mvar->beg1);
    mvar->end1 = std::max(nvar.end1, mvar->end1);

    bcf1_t *nv = bcf_dup(nvar.v);
    mvar->vs.push_back(nv);
    nvar.is_new_multiallelic = false;
    nvar.associated_new_multiallelic = mvar;

    return NULL;
}

/**
 * Creates a new multiallelic.
 */
void MultiallelicsConsolidator::update_multiallelic_for_printing(Variant& mvar)
{
//        std::cerr << "=========================\n";
//        std::cerr << "create new multi allelic\n";
//        std::cerr << "number of overlapping variants: " << mvar.vs.size() << "\n";
//        nvar.print();
//        cvar.print();
//        std::cerr << "NREW MULTI\n";
//        std::cerr <<  std::min(nvar.beg1, cvar.beg1) << "\n";
//        std::cerr <<  nvar.beg1 << "\n";
//        std::cerr <<  cvar.beg1 << "\n";
//        bcf_print(mvar->h, mvar->v);
//        std::cerr << "=========================\n";

    std::string ref;
    rs->fetch_seq(mvar.chrom, mvar.beg1, mvar.end1, ref);

    std::string new_alleles = ref;

    for (uint32_t i=0; i<mvar.vs.size(); ++i)
    {
        bcf_unpack(mvar.vs[i], BCF_UN_STR);
        std::string alt(bcf_get_alt(mvar.vs[i], 1));
        int32_t beg1 = bcf_get_pos1(mvar.vs[i]);
        int32_t end1 = bcf_get_end1(mvar.vs[i]);
        
        //extend alt
        std::string left = ref.substr(0, beg1-mvar.beg1);
        std::string right = ref.substr(end1-mvar.beg1+1, mvar.end1-end1);

//        std::cerr << "\talt  : " << alt << "\n";
//        std::cerr << "\txalt : " <<left << " " << alt  << " "  << right << "\n";

        new_alleles += "," + left + alt + right;

//        bcf_print_liten(mvar.h, mvar.vs[i]);
    }

//    std::cerr << "NEW ALLELES : " << new_alleles << "\n";
    bcf_set_pos1(mvar.v, mvar.beg1);
    bcf_update_alleles_str(mvar.h, mvar.v, new_alleles.c_str());
    
//     bcf_print(mvar.h, mvar.v);
     
//    std::cerr << "=========================\n";
        
        
        
//        VNTR& vntr = nvar.vntr;
//
//        //create a new copy of bcf1_t
//        bcf_hdr_t* h = nvar.h;
//        nvar.update_vntr_from_info_fields();
//        bcf1_t* nv = bcf_init1();
//        bcf_clear(nv);
//
//        if (vntr.exact_repeat_tract == "")
//        {
//            rs->fetch_seq(nvar.chrom, vntr.exact_beg1, vntr.exact_end1, vntr.exact_repeat_tract);
//        }
//
//        bcf_set_rid(nv, nvar.rid);
//        bcf_set_pos1(nv, vntr.exact_beg1);
//        kstring_t s = {0,0,0};
//        kputs(vntr.exact_repeat_tract.c_str(), &s);
//        kputc(',', &s);
//        kputs("<VNTR>", &s);
//        bcf_update_alleles_str(h, nv, s.s);
//        if (s.m) free(s.s);
//
//        if (no_samples) bcf_update_genotypes(h, nv, gts, no_samples);
//
//
//
//        Variant *nvntr = new Variant(h, nv);
//
//        std::string indel = bcf_variant2string(nvar.h, nvar.v);
//        nvntr->vntr.add_associated_indel(indel);
//
//        insert(nvntr);
//
////        bcf_print(h, nvar.v);
////        bcf_print(h, nv);
//
//        ++no_added_vntrs;
//    }
}