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
    refseq = new ReferenceSequence(ref_fasta_file);
}

/**
 * Inserts a VNTR record.
 * Returns true if successful.
 */
void MultiallelicsConsolidator::insert(Variant* var)
{
//    std::cerr << "buffer size : " << vbuffer.size() << "\n";

    flush(var);

    Variant& nvar = *var;

//    std::cerr << "inside insert\n";

//    std::cerr << "\t ";
//    bcf_print_liten(odr->hdr, var->v);

//    bcf_print(nvar.h, nvar.v);

    std::list<Variant*>::iterator i =vbuffer.begin();
    while(i != vbuffer.end())
    {
        Variant& cvar = **i;
//        std::cerr << "\t vs ";
//        bcf_print_liten(odr->hdr, cvar.v);
//vs 20:1000347:T/TATCCATCCATTCAACCATCCACCCACCCTCCC

        if (nvar.rid > cvar.rid)
        {
            vbuffer.insert(i, &nvar);
            return;
        }
        else if (nvar.rid == cvar.rid)
        {
            if (nvar.beg1 > cvar.beg1)
            {
               vbuffer.insert(i, &nvar);
               return;
            }
            else if (nvar.beg1 > cvar.beg1)
            {
                if (nvar.end1 > cvar.end1)
                {
                    vbuffer.insert(i, &nvar);
                    return;
                }
                else if (cvar.end1 == nvar.end1)
                {
                    vbuffer.insert(i, &nvar);
                    return;
                }
                else //nvar.rend1 < cvar.end1
                {
                    ++i;
                }
            }
            else //nvar.beg1 < cvar.beg1
            {
                ++i;
            }
        }
        else
        {
            fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
            exit(1);
        }
    }

    vbuffer.push_back(var);
}

/**
 * Flush variant buffer.
 */
void MultiallelicsConsolidator::flush(Variant* var)
{
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
            else
            {
                fprintf(stderr, "[%s:%d %s] Buffer is unordered\n", __FILE__, __LINE__, __FUNCTION__);
                exit(1);
            }

            ++i;
        }

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
        create_new_multiallelic(*var);
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
 * Creates a new multiallelic.
 */
void MultiallelicsConsolidator::create_new_multiallelic(Variant& nvar)
{
//    //convert all indels into their corresponding VNTR representation using exact parameters
//    if (nvar->new_multiallelic)
//    {
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
//            refseq->fetch_seq(nvar.chrom, vntr.exact_beg1, vntr.exact_end1, vntr.exact_repeat_tract);
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