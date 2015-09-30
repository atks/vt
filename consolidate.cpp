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

#include "consolidate.h"

namespace
{

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    bool merge_by_pos;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    ////////////////
    //variant buffer
    ////////////////
    std::list<Variant *> variant_buffer; //front is most recent

    ////////////
    //filter ids
    ////////////
    char* overlap_snp;
    char* overlap_indel;
    char* overlap_vntr;
    int32_t overlap_snp_id;
    int32_t overlap_indel_id;
    int32_t overlap_vntr_id;
    
    /////////
    //stats//
    /////////
    int32_t no_total_variants;
    int32_t no_nonoverlap_variants;
    int32_t no_overlap_variants;
    int32_t no_new_multiallelic_snps;
    int32_t no_new_multiallelic_indels;
    
    /////////
    //tools//
    /////////
    VariantManip *vm;

    Igor(int argc, char **argv)
    {
        version = "0.57";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Consolidates variants in the following ways\n"
                 "              1.Updates each variant's filter field if it overlaps with another variant.\n"
                 "                    SNP   - overlap at position \n"
                 "                    Indel - overlap with region bounded by exact flanks\n"
                 "                    VNTR  - overlap with region bounded by fuzzy flanks\n"
                 "              2.Adds candidate multiallelic SNPs that do not overlap with Indels or VNTRs\n"
                 "              3.Adds candidate multiallelic Indels if and only if\n"
                 "                    a. they do not overlap with SNPs and VNTRs\n"
                 "                    b. the exact flanks and fuzzy flanks are equal to one another\n"
                 "              4.Adds INFO fields indicating the number of variants that overlap with this variant";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        odw->link_hdr(odr->hdr);
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_snp,Description=\"Overlaps with snp\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_vntr,Description=\"Overlaps with VNTR\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_indel,Description=\"Overlaps with indel\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=shorter_vntr,Description=\"A longer VNTR overlaps with this VNTR\">");
        odw->write_hdr();

        overlap_snp = const_cast<char*>("overlap_snp");
        overlap_indel = const_cast<char*>("overlap_indel");
        overlap_vntr = const_cast<char*>("overlap_vntr");

        overlap_snp_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_snp");
        overlap_indel_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_indel");
        overlap_vntr_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr");

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_total_variants = 0;
        no_nonoverlap_variants = 0;
        no_overlap_variants = 0;
        no_new_multiallelic_snps = 0;
        no_new_multiallelic_indels = 0;
    
        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip();
    }

    /**
     * Inserts a Variant record.
     */
    void insert_variant_record_into_buffer(Variant* variant)
    {
        char* tr = NULL;
        int32_t n = 0;
        
        //this is for indels that are annotated with TR - you cannot
        //detect by overlap of the reference sequence because in the normalized
        //representation for insertions, only the anchor base is present.
        //This is not an issue for deletions.
        //reminder: consolidate is for a VCF file that is produced by annotate_indels.
        if (bcf_get_info_string(odw->hdr, variant->v, "TR", &tr, &n)>0)
        {
            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
            free(tr);
        }
        
        //update filters
        std::list<Variant *>::iterator i = variant_buffer.begin();
        while(i!=variant_buffer.end())
        {
            Variant *cvariant = *i;

            if (variant->rid > cvariant->rid)
            {
                break;
            }
            else if (variant->rid == cvariant->rid)
            {
                if (variant->end1 < cvariant->pos1) //not possible
                {
                    fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
                    exit(1);
                }
                else if (variant->pos1 > cvariant->end1) //does not overlap
                {
                    break;
                }
                else //overlaps
                {
                    ///////
                    //SNP//
                    ///////
                    if (variant->type==VT_SNP)
                    {
                        if (cvariant->type==VT_SNP)
                        {    
                            //compile into multiallelic
                            //
                            if (bcf_get_n_filter(cvariant->v)==0)
                            {
                                //make a new a multiallelic
                                Variant* multiallelic_variant = new Variant(cvariant, variant);
                                multiallelic_variant->no_overlapping_snps = 2;
                                variant_buffer.push_front(multiallelic_variant);
                            }
                                                        
                            bcf_add_filter(odw->hdr, variant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_snp_id);
                            ++cvariant->no_overlapping_snps;
                        }   
                        else if (cvariant->type==VT_INDEL)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_indel_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_snp_id);
                            ++cvariant->no_overlapping_snps;
                        }
                        else if (cvariant->type==VT_VNTR)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
                            ++variant->no_overlapping_vntrs;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                        }
                        else if (cvariant->type==VT_UNDEFINED)
                        {        
                            variant->vs.push_back(bcf_dup(variant->v));
                            ++variant->no_overlapping_snps;
                        }                        
                    }
                    /////////
                    //INDEL//
                    /////////
                    else if (variant->type==VT_INDEL)
                    {
                        if (cvariant->type==VT_SNP)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_indel_id);
                            ++cvariant->no_overlapping_indels;
                        }
                        else if (cvariant->type==VT_INDEL)
                        {    
                            if (bcf_get_n_filter(cvariant->v) == 0)
                            {
                                //make a new a multiallelic

                            }
                                          
                            bcf_add_filter(odw->hdr, variant->v, overlap_snp_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_indel_id);
                            ++cvariant->no_overlapping_indels;
                        }
                        else if (cvariant->type==VT_VNTR)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
                            ++variant->no_overlapping_vntrs;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_indel_id);
                            ++cvariant->no_overlapping_indels;
                        }
                        else if (cvariant->type==VT_UNDEFINED)
                        {    

                        }                        
                    }
                    ////////
                    //VNTR//
                    ////////
                    else if (variant->type==VT_VNTR)
                    {
                        if (cvariant->type==VT_SNP)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_snp_id);
                            ++variant->no_overlapping_snps;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_vntr_id);                        
                            ++cvariant->no_overlapping_vntrs;
                        }
                        else if (cvariant->type==VT_INDEL)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_indel_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_vntr_id);
                            ++cvariant->no_overlapping_vntrs;
                        }
                        else if (cvariant->type==VT_VNTR)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_vntr_id);
                            ++variant->no_overlapping_indels;
                            bcf_add_filter(odw->hdr, cvariant->v, overlap_vntr_id);
                            ++cvariant->no_overlapping_vntrs;
                            
                            //mark a better VNTR?? by score and by length?
                            
                        }
                    }
                                        
                    ++i;
                }
            }
            else //variant.rid < cvariant.rid is impossible if input file is ordered.
            {
                
                fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
                exit(1);
            }
        }

        variant_buffer.push_front(variant);
    }

    /**
     * Flush variant buffer.
     */
    void flush_variant_buffer(bcf1_t* v)
    {
        if (variant_buffer.empty())
        {
            return;
        }

        int32_t rid = bcf_get_rid(v);
        int32_t pos1 = bcf_get_pos1(v);

        while (!variant_buffer.empty())
        {
            Variant* variant = variant_buffer.back();
            
            if (variant->rid < rid) 
            {
                if (variant->type==VT_UNDEFINED)
                {
                    if (consolidate_multiallelic(variant))
                    {
                        odw->write(variant->v);
                        variant->v = NULL;
                        delete variant;
                        variant_buffer.pop_back();
                    }
                }    
                else
                {    
                    odw->write(variant->v);
                    variant->v = NULL;
                    delete variant;
                    variant_buffer.pop_back();
                }
            }
            else if (variant->rid == rid) 
            {
                if (variant->pos1 < pos1-1000)
                {
                    odw->write(variant->v);
                    variant->v = NULL;
                    delete variant;
                    variant_buffer.pop_back();
                }
                else
                {
                    break;
                }
            }
        }
    }

    /**
     * Consolidate multiallelic variant based on associated biallelic records
     * stored in vs.  Updates v which is to be the consolidated multiallelic
     * variant.
     *
     * returns true if the multiallelic variant is good to go.
     */
    bool consolidate_multiallelic(Variant* variant)
    {
        if (variant->no_overlapping_snps !=0 &&
            variant->no_overlapping_indels ==0 &&
            variant->no_overlapping_vntrs ==0)
        {
            bcf1_t* v = bcf_init1();
            
            std::vector<bcf1_t*>& vs = variant->vs;
            
            bcf_set_rid(v, bcf_get_rid(vs[0]));
            bcf_set_pos1(v, bcf_get_pos1(vs[0]));
            
            kstring_t s = {0,0,0};
            kputc(bcf_get_snp_ref(vs[0]), &s);
            kputc(',', &s);
            int32_t no_alleles = vs.size();
            char alts[no_alleles];
            for (uint32_t i=0; i<no_alleles; ++i)
            {
                alts[i] = bcf_get_snp_alt(vs[i]);
            }
            //selection sort
            for (uint32_t i=0; i<no_alleles-1; ++i)
            {
                for (uint32_t j=i+1; j<no_alleles; ++j)
                {
                    if (alts[j]<alts[i])
                    {
                        char tmp = alts[j];
                        alts[j] = alts[i];
                        alts[i] = tmp;
                    }    
                    alts[i] = bcf_get_snp_alt(vs[i]);
                }
                
                kputc(alts[i], &s);
                kputc(',', &s);
            }
            kputc(alts[no_alleles-1], &s);
            bcf_update_alleles_str(odw->hdr, v, s.s);
               
            ++no_new_multiallelic_snps;
            return true;
        }
        else if (variant->no_overlapping_snps ==0 &&
                 variant->no_overlapping_indels !=0 &&
                 variant->no_overlapping_vntrs ==0)
        {
            
//                //make multiallelic
//                
//                ++no_new_multiallelic_indels;
//                return true;                
//            
        }    
        
        return false;
    }

    /**
     * Flush variant buffer.
     */
    void flush_variant_buffer()
    {
        while (!variant_buffer.empty())
        {
            Variant* variant = variant_buffer.back();
            odw->write(variant->v);
            delete variant;
            variant_buffer.pop_back();
        }
    }

    void consolidate()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();

        Variant* variant;
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            flush_variant_buffer(v);
                       
            variant = new Variant(odw->hdr, v);
            vm->classify_variant(odw->hdr, v, *variant);
            insert_variant_record_into_buffer(variant);

            v = odw->get_bcf1_from_pool();

            ++no_total_variants;
        }

        flush_variant_buffer();

        odr->close();
        odw->close();
    };

    void print_options()
    {
        std::clog << "consolidate_variants v" << version << "\n\n";

        std::clog << "options:     input VCF file        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        if (intervals.size()!=0)
        {
            std::clog << "         [i] intervals             " << intervals.size() <<  " intervals\n";
        }
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: Total number of observed variants    " << no_total_variants << "\n";
        std::clog << "       Total number of nonoverlap variants  " << no_nonoverlap_variants << "\n";
        std::clog << "       Total number of overlap variants     " << no_overlap_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

void consolidate(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.consolidate();
    igor.print_stats();
};
