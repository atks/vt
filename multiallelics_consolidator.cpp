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

MultiallelicsConsolidator::MultiallelicsConsolidator(std::string& input_vcf_file, std::vector<GenomeInterval>& intervals, std::string& output_vcf_file, std::string& fexp, std::string& ref_fasta_file, int32_t window_overlap)
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
    END = bcf_hdr_append_info_with_backup_naming(odw->hdr, "END", "1", "Integer", "End position of the variant.", rename);
    MERGED_NEW_MULTIALLELIC_VARIANT = bcf_hdr_append_info_with_backup_naming(odw->hdr, "MERGED_NEW_MULTIALLELIC_VARIANT", ".", "String", "New multiallelic variant that this variant is merged into.", rename);
    CONSTITUENT_VARIANTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "CONSTITUENT_VARIANTS", ".", "String", "Variants that are merged into this multiallelic variant.", rename);
    odw->write_hdr();

    //
    //
    //
    this->window_overlap = window_overlap;

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
    if (mvar)
    {
        if (mvar->vs.size()>1) vbuffer.remove(mvar);
        vbuffer.push_front(mvar);
    }
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

        //search for point to start deleting from
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
    delete var;
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
        //create multiallelic
        Variant *mvar = new Variant();
        mvar->h = nvar.h;
        mvar->v = bcf_init();
        bcf_clear(mvar->v);
        std::string s = "NNN,<MULTI>";
        bcf_update_alleles_str(mvar->h, mvar->v, s.c_str());

        mvar->chrom = nvar.chrom;
        mvar->rid = cvar.rid;
        mvar->beg1 = std::min(nvar.beg1, cvar.beg1);
        bcf_set_rid(mvar->v, mvar->rid);
        bcf_set_pos1(mvar->v, mvar->beg1);
        mvar->end1 = std::max(nvar.end1, cvar.end1);
        mvar->vs.push_back(cvar.v);
        mvar->vs.push_back(nvar.v);
        mvar->is_new_multiallelic = true;

        //update biallelics
        nvar.is_involved_in_a_multiallelic = true;
        nvar.associated_new_multiallelic = mvar;
        cvar.is_involved_in_a_multiallelic = true;
        cvar.associated_new_multiallelic = mvar;

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

//    std::cerr << "=============================\n";
//    std::cerr << "UPDATE MULTIALLELIC\n";
//    nvar.print();
//    cvar.print();

    //update existing multiallelic
    mvar->beg1 = std::min(nvar.beg1, mvar->beg1);
    mvar->end1 = std::max(nvar.end1, mvar->end1);
    mvar->vs.push_back(nvar.v);

//    std::cerr << "mvar.beg1 : " << mvar->beg1 << "\n";
//    std::cerr << "mvar.end1 : " << mvar->end1 << "\n";
//    std::cerr << "nvar.beg1 : " << nvar.beg1 << "\n";
//    std::cerr << "nvar.end1 : " << nvar.end1 << "\n";
//    std::cerr << "cvar.beg1 : " << cvar.beg1 << "\n";
//    std::cerr << "cvar.end1 : " << cvar.end1 << "\n";
//    bcf_print(mvar->h, mvar->v);
//    bcf_print(mvar->h, nvar.v);
//    bcf_print(mvar->h, cvar.v);
//    std::cerr << "no. overlapping variants : " << mvar->vs.size() << "\n";
//    std::cerr << "=============================\n";

    //update new record to point to multiallelic
    nvar.associated_new_multiallelic = mvar;

    return mvar;
}

/**
 * Creates a new multiallelic.
 *
 * Each alternate allele only contains an instance of one of the alternative alleles observed.
 *
 * todo: to add combinations of the alternate alleles to be merged
 *
 * The current implementation is works fine for pure SNPs and Tandem repeats
 *    
 *    SNP  :   A/G  + A/T    => A/G/T
 *    VNTR :   TG/T + T/TGGG => TG/T/TGGG
 *
 * But in the case of complex variants
 *
 *   SNP+INDEL : A/G  + A/AT  => A/G/AT
 *
 *  This is not sufficient as both alternate alleles may coexist
 *   
 *   SNP+INDEL : A/G  + A/AT  => A/G/AT/GT
 */
void MultiallelicsConsolidator::update_multiallelic_for_printing(Variant& mvar)
{
//    std::cerr << "================================\n";
//    std::cerr << "UPDATE MULTIALLELICF OR PRINTING\n";
//    std::cerr << "mvar.beg1 : " << mvar.beg1 << "\n";
//    std::cerr << "mvar.end1 : " << mvar.end1 << "\n";
//    bcf_print(mvar.h, mvar.v);
//    std::cerr << "no. overlapping variants : " << mvar.vs.size() << "\n";
//    std::cerr << "================================\n";

    std::string ref;
    rs->fetch_seq(mvar.chrom, mvar.beg1, mvar.end1, ref);

    std::vector<std::string> alleles;
    alleles.push_back(ref);
    std::string biallelics = "";

    for (uint32_t i=0; i<mvar.vs.size(); ++i)
    {
        bcf_unpack(mvar.vs[i], BCF_UN_STR);
        std::string alt(bcf_get_alt(mvar.vs[i], 1));
        int32_t beg1 = bcf_get_pos1(mvar.vs[i]);
        int32_t end1 = bcf_get_end1(mvar.vs[i]);

        //extend alt
        std::string left = ref.substr(0, beg1-mvar.beg1);
        std::string right = ref.substr(end1-mvar.beg1+1, mvar.end1-end1);

        alleles.push_back(left + alt + right);

        if (biallelics.size()!=0) biallelics += ",";
        biallelics +=  bcf_variant2string(odw->hdr, mvar.vs[i]);
    }

    //trim
    int32_t left_trimmed = 0;
    VariantManip::left_trim(alleles, mvar.beg1, left_trimmed);

    bcf_set_rid(mvar.v, mvar.rid);
    bcf_set_pos1(mvar.v, mvar.beg1);
    bcf_update_info_int32(odw->hdr, mvar.v, "END", &mvar.end1, 1);
    std::string new_alleles = join(alleles, ",");
    bcf_update_alleles_str(mvar.h, mvar.v, new_alleles.c_str());
    bcf_update_info_string(odw->hdr, mvar.v, CONSTITUENT_VARIANTS.c_str(), biallelics.c_str());

    mvar.updated_multiallelic = true;

    std::string multi = bcf_variant2string(odw->hdr, mvar.v);
    for  (uint32_t i=0; i<mvar.vs.size(); ++i)
    {
        bcf_update_info_string(odw->hdr, mvar.vs[i], MERGED_NEW_MULTIALLELIC_VARIANT.c_str(), multi.c_str());
    }

//    std::cerr << "================================\n";
//    std::cerr << "UPDATE MULTIALLELIC FOR PRINTING END\n";
//    std::cerr << "mvar.beg1 : " << mvar.beg1 << "\n";
//    std::cerr << "mvar.end1 : " << mvar.end1 << "\n";
//    bcf_print(mvar.h, mvar.v);
//    std::cerr << "no. overlapping variants : " << mvar.vs.size() << "\n";
//    std::cerr << "==================================\n";
}