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

#include "sort.h"

namespace
{

//sort and write out
int bcf_compare (const void * a, const void * b)
{
    bcf1_t *u = *(bcf1_t**) a;
    bcf1_t *v = *(bcf1_t**) b;

    if (bcf_get_rid(u)==bcf_get_rid(v))
    {
        return (bcf_get_pos1(u)-bcf_get_pos1(v));
    }

    return (bcf_get_rid(u)-bcf_get_rid(v));
}

class Igor : Program
{
    public:

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    uint32_t sort_window_size;
    std::string sort_mode;
    bool print;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    uint32_t no_variants;

    /////////
    //tools//
    /////////

    Igor(int argc, char **argv)
    {
        version = "0.5";

        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Sorts a VCF or BCF or VCF.GZ file.\n";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my;
            cmd.setOutput(&my);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::ValueArg<uint32_t> arg_sort_window_size("w", "w", "local sorting window size, set by default to 1000 under local mode. [0]", false, 0, "int", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF/VCF.GZ/BCF file [-]"
                   , false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_sort_mode("m", "m", ""
                               "sorting modes. [full]\n"
                 "              local : locally sort within a 1000bp window.  Window size may be set by -w.\n"
                 "              chrom : sort chromosomes based on order of contigs in header.\n"
                 "                      input must be indexed\n"
                 "              full  : full sort with no assumptions",
                                 false, "full", "str", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            sort_mode = arg_sort_mode.getValue();
            print = arg_print.getValue();
            sort_window_size = arg_sort_window_size.getValue();

            if (sort_mode=="local")
            {
                sort_window_size = sort_window_size ? sort_window_size : 1000;
            }
            else if (sort_mode=="chrom")
            {
                //do nothing.  A check that the index file is present is made later in sort()
            }
            else if (sort_mode=="full")
            {
                //do nothing
            }
            else
            {
                fprintf(stderr, "[%s:%d %s] Sort mode can only be full, local or chrom.\n", __FILE__,__LINE__,__FUNCTION__);
            }
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

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_variants = 0;

        ///////////////////////
        //tool initialization//
        ///////////////////////
    }

    void sort()
    {
        if (sort_mode=="local")
        {
            odr = new BCFOrderedReader(input_vcf_file, intervals);

            odw = new BCFOrderedWriter(output_vcf_file, sort_window_size);
            odw->link_hdr(odr->hdr);
            odw->write_hdr();

            bcf1_t *v = odw->get_bcf1_from_pool();
            bcf_hdr_t *h = odr->hdr;

            while (odr->read(v))
            {
                odw->write(v);
                v = odw->get_bcf1_from_pool();
                ++no_variants;

            }

            odw->close();
            odr->close();
        }
        else if (sort_mode=="chrom")
        {
            odr = new BCFOrderedReader(input_vcf_file, intervals);
            if (!odr->is_index_loaded())
            {
                fprintf(stderr, "[%s:%d %s] Chromosome sort mode requires that %s is an indexed vcf.gz file.\n", __FILE__,__LINE__,__FUNCTION__, input_vcf_file.c_str());
                exit(1);
            }

            odw = new BCFOrderedWriter(output_vcf_file, sort_window_size);
            odw->link_hdr(odr->hdr);
            odw->write_hdr();

            int32_t nseqs;
            const char ** seqs = bcf_hdr_seqnames(odr->hdr, &nseqs);
            bcf1_t *v = bcf_init1();
            bcf_hdr_t *h = odr->hdr;

            for (size_t i=0; i<nseqs; ++i)
            {
                std::string interval(seqs[i]);
                GenomeInterval ginterval(interval);
                odr->jump_to_interval(ginterval);

                while (odr->read(v))
                {
                    odw->write(v);
                    ++no_variants;
                }
            }

            bcf_destroy(v);
            odw->close();
            odr->close();
        }
        else if (sort_mode=="full")
        {
            //read into buffer 10000 records
            odr = new BCFOrderedReader(input_vcf_file, intervals);

            //need to change this to measure the amount of data read.  This can be
            size_t buffer_size = 1000000;

            bcf1_t *buffer[buffer_size];
            for (size_t i=0; i<buffer_size; ++i)
            {
                buffer[i] = bcf_init1();
            }

            size_t bptr = 0;
            std::vector<std::string> sorted_file_names;

            bcf1_t *v = buffer[bptr];
            while (odr->read(v))
            {
                if (bptr<buffer_size-1)
                {
                    ++bptr;
                }
                else
                {
                    ++bptr;

                    //sort and write out
                    qsort (buffer, bptr, sizeof(bcf1_t*), bcf_compare);

                    kstring_t s = {0,0,0};
                    if (output_vcf_file!="-")
                    {
                        kputs(output_vcf_file.c_str(), &s);
                        kputs(".", &s);
                    }
                    kputw(sorted_file_names.size()+1, &s);
                    kputs(".bcf", &s);
                    std::string file_name(s.s);
                    if (s.m) free(s.s);
                    sorted_file_names.push_back(file_name);

                    if (sorted_file_names.size()!=1)
                    {
                        odw->close();
                        delete odw;
                    }
                    odw = new BCFOrderedWriter(file_name);
                    odw->link_hdr(odr->hdr);
                    odw->write_hdr();

                    for (size_t i=0; i<bptr; ++i)
                    {
                        odw->write(buffer[i]);
                    }

                    bptr = 0;
                }

                v = buffer[bptr];
                ++no_variants;
            }

            if (bptr)
            {
                qsort (buffer, bptr, sizeof(bcf1_t*), bcf_compare);

                kstring_t s = {0,0,0};
                if (output_vcf_file!="-")
                {
                    kputs(output_vcf_file.c_str(), &s);
                    kputs(".", &s);
                }
                kputw(sorted_file_names.size()+1, &s);
                kputs(".bcf", &s);
                std::string file_name(s.s);
                if (s.m) free(s.s);
                sorted_file_names.push_back(file_name);

                if (sorted_file_names.size()!=1)
                {
                    odw->close();
                    delete odw;
                }
                odw = new BCFOrderedWriter(file_name);
                odw->link_hdr(odr->hdr);
                odw->write_hdr();

                for (size_t i=0; i<bptr; ++i)
                {
                    odw->write(buffer[i]);
                }

                odw->close();
                delete odw;

                bptr = 0;
            }

            //merge records from temporary files
            intervals.clear();
            BCFSyncedReader sr(sorted_file_names, intervals);

            std::vector<bcfptr*> current_recs;
            odw = new BCFOrderedWriter(output_vcf_file);
            odw->link_hdr(odr->hdr);
            odw->write_hdr();

            int32_t no = 0;
            while(sr.read_next_position(current_recs))
            {
                for (size_t i=0; i<current_recs.size(); ++i)
                {
                    odw->write(current_recs[i]->v);
                }

                ++no;
            }

            sr.close();

            for (size_t i=0; i< sorted_file_names.size(); ++i)
            {
                std::remove(sorted_file_names[i].c_str());
            }

            odr->close();

            odw->close();
            delete odw;

            for (size_t i=0; i<buffer_size; ++i)
            {
                bcf_destroy(buffer[i]);
            }
        }
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "sort v" << version << "\n\n";

        std::clog << "options:     input VCF file              " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file             " << output_vcf_file << "\n";
        std::clog << "         [w] sort window size            " << sort_window_size << "\n";
        std::clog << "         [m] sorting mode                " << sort_mode << "\n";
        std::clog << "         [p] print options and stats     " << (print ? "yes" : "no") << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::clog << "stats: no. variants  : " << no_variants << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool sort(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.sort();
    igor.print_stats();
    return igor.print;
};
