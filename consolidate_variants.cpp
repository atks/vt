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

#include "consolidate_variants.h"

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
    int32_t overlap_indel;
    int32_t overlap_vntr;
    
    /////////
    //stats//
    /////////
    int32_t no_total_variants;
    int32_t no_nonoverlap_variants;
    int32_t no_overlap_variants;

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
            std::string desc = "Consolidates variants"
                "Annotates overlapping Indels with STRs";

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
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_vntr,Description=\"Overlaps with VNTR\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_indel,Description=\"Overlaps with indel\">");
        odw->write_hdr();

        overlap_indel = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_indel");
        overlap_vntr = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_vntr");

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_total_variants = 0;
        no_nonoverlap_variants = 0;
        no_overlap_variants = 0;

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
        if (bcf_get_info_string(odw->hdr, variant->v, "TR", &tr, &n)>0)
        {
            bcf_add_filter(odw->hdr, variant->v, overlap_vntr);
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
                    if (variant->type==VT_SNP)
                    {
                        if (cvariant->type==VT_INDEL)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_indel);
                        }
                        else if (cvariant->type==VT_VNTR)
                        {    
                            bcf_add_filter(odw->hdr, variant->v, overlap_vntr);
                        }
                    }
                    else if (variant->type==VT_INDEL)
                    {
                        bcf_add_filter(odw->hdr, cvariant->v, overlap_indel);
                    }
                    else if (variant->type==VT_VNTR)
                    {
                        bcf_add_filter(odw->hdr, cvariant->v, overlap_vntr);
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
                odw->write(variant->v);
                delete variant;
                variant_buffer.pop_back();
            }
            else if (variant->rid == rid) 
            {
                if (variant->pos1 < pos1-1000)
                {
                    odw->write(variant->v);
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

    void consolidate_variants()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();

        Variant* variant;
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);

            
            flush_variant_buffer(v);
            
                       
            variant = new Variant(v);
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

void consolidate_variants(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.consolidate_variants();
    igor.print_stats();
};
