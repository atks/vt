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
    std::list<Variant> variant_buffer; //front is most recent
        
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
        odw->write_hdr();

        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap-vntr,Description=\"Overlaps with VNTR\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap-indel,Description=\"Overlaps with indel\">");
        
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=olap,Description=\"Overlapping Alleles\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=snpstr,Description=\"SNP in STR\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=badmotif,Description=\"Poorly defined motif\">");


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

//    /**
//     * Inserts a VNTR record.
//     * Returns true if successful.
//     */
//    bool insert_variant_record_into_buffer(Variant& variant)
//    {
//        std::list<VNTR>::iterator i = variant_buffer.begin();
//        while(i!=variant_buffer.end())
//        {
//            VNTR& cvntr = *i;
//
//            if (vntr.rid > cvntr.rid)
//            {
//                variant_buffer.insert(i, vntr);
//                return true;
//            }
//            else if (vntr.rid == cvntr.rid)
//            {
//                if (vntr.rbeg1 > cvntr.rbeg1)
//                {
//                    variant_buffer.insert(i, vntr);
//                    return true;
//                }
//                else if (vntr.rbeg1 == cvntr.rbeg1)
//                {
//                    if (vntr.rend1 > cvntr.rend1)
//                    {
//                        variant_buffer.insert(i, vntr);
//                        return true;
//                    }
//                    else if (cvntr.rend1 == vntr.rend1)
//                    {
//                        if (cvntr.motif > vntr.motif)
//                        {
//                            variant_buffer.insert(i, vntr);
//                            return true;
//                        }
//                        else if (cvntr.motif == vntr.motif)
//                        {
//                            //do not insert
//                         
////                            std::cerr << "NEVER inseert\n";
//                            return false;
//                        }
//                        else // cvntr.motif > vntr.motif
//                        {
//                            ++i;
//                        }
//                    }
//                    else // cvntr.rend1 > vntr.rend1
//                    {
//                        ++i;
//                    }
//                }
//                else //vntr.rbeg1 < cvntr.rbeg1
//                {
//                    ++i;
//                }
//            }
//            else //vntr.rid < cvntr.rid is impossible if input file is ordered.
//            {
//                fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
//                exit(1);
//            }
//        }
//
//        variant_buffer.push_back(vntr);
//        return true;
//    }
//
//    /**
//     * Flush variant buffer.
//     */
//    void flush_variant_buffer(bcf1_t* v)
//    {
//        if (variant_buffer.empty())
//        {
//            return;
//        }
//
//        int32_t rid = bcf_get_rid(v);
//        int32_t pos1 = bcf_get_pos1(v);
//
//        //search for vntr to start deleting from.
//        std::list<Variant>::iterator i = variant_buffer.begin();
//        while(i!=variant_buffer.end())
//        {
//            VNTR& vntr = *i;
//
////            std::cerr << vntr.rid << " " << rid << "\n";
//
//            if (vntr.rid < rid)
//            {
//                break;
//            }
//            else if (vntr.rid == rid)
//            {
//                if (vntr.rend1 < pos1-1000)
//                {
//                    break;
//                }
//            }
//            else //rid < vntr.rid is impossible
//            {
//                fprintf(stderr, "[%s:%d %s] File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
//                exit(1);
//            }
//
//            ++i;
//        }
//
//        while (i!=variant_buffer.end())
//        {
//            i = variant_buffer.erase(i);
//        }
//    }
    
    void consolidate_variants()
    {
        bcf1_t *v = odw->get_bcf1_from_pool();

        Variant* variant;
        while (odr->read(v))
        {
            bcf_unpack(v, BCF_UN_STR);
            
            variant = new Variant();
            
            

            
            ++no_total_variants;
        }


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
