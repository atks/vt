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

#include "paste.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_vcf_files;
    std::string input_vcf_file_list;
    std::string output_vcf_file;
    bool print;

    ///////
    //i/o//
    ///////
    std::vector<BCFOrderedReader *> odrs;
    BCFOrderedWriter *odw;

    ///////////////
    //general use//
    ///////////////

    /////////
    //stats//
    /////////

    /////////
    //tools//
    /////////

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Pastes VCF files like the unix paste functions.\n"
                 "              This is used after the per sample genotyping step in vt.\n"
                 "              Input requirements and assumptions:\n"
                 "              1. Same variants are represented in the same order for each file (required)\n"
                 "              2. Genotype fields are the same for corresponding records (required)\n"
                 "              3. Sample names are different in all the files (warning will be given if not)\n"
                 "              4. Headers (not including the samples) are the same for all the files (unchecked assumption, will fail if output is BCF)\n"
                 "              Outputs:\n"
                 "              1. INFO fields output will be that of the first file\n"
                 "              2. Genotype fields are the same for corresponding records\n";


            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);

            cmd.parse(argc, argv);

            input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            const std::vector<std::string>& v = arg_input_vcf_files.getValue();
            print = arg_print.getValue();

            if (input_vcf_files.size()==0)
            {
                fprintf(stderr, "[E:%s:%d %s] no input vcf files.\n", __FILE__, __LINE__, __FUNCTION__);
                exit(1);
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
        std::vector<GenomeInterval> intervals;
        for (size_t i=0; i<input_vcf_files.size(); ++i)
        {
            odrs.push_back(new BCFOrderedReader(input_vcf_files[i], intervals));
        }
        odw = new BCFOrderedWriter(output_vcf_file, 0);
        odw->set_hdr(odrs[0]->hdr);

        ///////////////
        //general use//
        ///////////////

        ////////////////////////
        //stats initialization//
        ////////////////////////

        /////////
        //tools//
        /////////
    }

    void paste()
    {
        int32_t nfiles = odrs.size();

        //add all sample names to output vcf header and warn if there are more than one occurence of a sample
        int32_t no_samples = bcf_hdr_nsamples(odw->hdr);
        for (size_t i=1; i<nfiles; ++i)
        {
            for (size_t j=0; j<bcf_hdr_nsamples(odrs[i]->hdr); ++j)
            {
                if (bcf_hdr_id2int(odw->hdr, BCF_DT_SAMPLE, bcf_hdr_get_sample_name(odrs[i]->hdr, j))==-1)
                {
                    bcf_hdr_add_sample(odw->hdr, bcf_hdr_get_sample_name(odrs[i]->hdr, j));
                    ++no_samples;
                }
                else
                {
                    fprintf(stderr, "[E:%s:%d %s] %s is present more than once.\n", __FILE__, __LINE__, __FUNCTION__, bcf_hdr_get_sample_name(odrs[i]->hdr, j));
                }
            }
        }
        bcf_hdr_add_sample(odw->hdr, NULL);

        odw->write_hdr();

        std::vector<bcfptr*> current_recs;
        int ncount =0;

        std::vector<bcfptr*> sample2record(no_samples, NULL);
        std::vector<int32_t> sample2index(no_samples, -1);

        bcf1_t* vs[nfiles];
        bcf1_t *nv = bcf_init();

        for (size_t i =0; i<nfiles; ++i)
        {
            vs[i]= bcf_init();
        }

        std::vector<std::string> genotype_fields;

        while(true)
        {
            size_t read_no = 0;
            for (size_t i =0; i<nfiles; ++i)
            {
                if (odrs[i]->read(vs[i]))
                {
                    bcf_unpack(vs[i], BCF_UN_FMT);
                    ++read_no;
                }
            }

            if (read_no!=nfiles)
            {
                break;
            }

            bcf_copy(nv, vs[0]);
            int32_t nval_per_sample[vs[0]->n_fmt];
            int32_t size[vs[0]->n_fmt];
            //check consistency of FORMAT
            for (size_t i=1; i<nfiles; ++i)
            {
                if (vs[0]->n_fmt==vs[i]->n_fmt)
                {
                    for (size_t j=0; j<vs[0]->n_fmt; ++j)
                    {
                        const char* a = bcf_get_format(odrs[0]->hdr, vs[0], j);
                        const char* b = bcf_get_format(odrs[i]->hdr, vs[i], j);
                        if (i==1)
                        {
                            nval_per_sample[j] = vs[0]->d.fmt[j].n > vs[1]->d.fmt[j].n ? vs[0]->d.fmt[j].n : vs[1]->d.fmt[j].n;
                        }
                        else
                        {
                            nval_per_sample[j] = nval_per_sample[j] > vs[i]->d.fmt[j].n ? nval_per_sample[j] : vs[i]->d.fmt[j].n;
                        }
                        if (i==1)
                        {
                            size[j] = vs[0]->d.fmt[j].size > vs[1]->d.fmt[j].size ? vs[0]->d.fmt[j].size : vs[1]->d.fmt[j].size;
                        }
                        else
                        {
                            size[j] = size[j] > vs[i]->d.fmt[j].size ? size[j] : vs[i]->d.fmt[j].size;
                        }
                        if (strcmp(a,b))
                        {
                            fprintf(stderr, "[E:%s:%d %s] FORMAT not consistent between files %s %s.\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files[0].c_str(),  input_vcf_files[i].c_str());
                            exit(1);
                        }
                    }
                }
                else
                {
                    fprintf(stderr, "[E:%s:%d %s] FORMAT not consistent between files %s %s.\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files[0].c_str(),  input_vcf_files[i].c_str());
                    exit(1);
                }
            }

            //for each format, construct array
            for (size_t i=0; i<vs[0]->n_fmt; ++i)
            {
                const char* genotype_field = bcf_get_format(odrs[0]->hdr, vs[0], i);

                if (vs[0]->d.fmt[i].type == BCF_BT_INT8 ||
                    vs[0]->d.fmt[i].type == BCF_BT_INT16 ||
                    vs[0]->d.fmt[i].type == BCF_BT_INT32)
                {
                    int32_t ndst = no_samples * nval_per_sample[i] * sizeof(int32_t);

                    int32_t *data = (int32_t *) malloc(ndst);
                    int32_t *p = data;
                    int32_t np = no_samples * nval_per_sample[i];

                    for (size_t j=0; j<nfiles; ++j)
                    {
                        int32_t b = bcf_get_format_int32(odrs[j]->hdr, vs[j], genotype_field, &p, &np);
                        p = p+bcf_hdr_nsamples(odrs[j]->hdr)*vs[j]->d.fmt[i].n;
                        np = np-bcf_hdr_nsamples(odrs[j]->hdr)*nval_per_sample[i];

                        if (vs[j]->d.fmt[i].n<nval_per_sample[i])
                        {
                            int32_t steps = bcf_hdr_nsamples(odrs[j]->hdr)*(nval_per_sample[i] - vs[j]->d.fmt[i].n);
                            while (steps)
                            {
                                p[0] = *(p-1);
                                ++p;
                                --steps;
                            }
                        }
                    }

                    bcf_update_format_int32(odw->hdr, nv, genotype_field, data, no_samples * vs[0]->d.fmt[i].n);
                    free(data);
                }
                else if (vs[0]->d.fmt[i].type == BCF_BT_FLOAT)
                {
                    int32_t ndst = no_samples * nval_per_sample[i] * sizeof(float);

                    float *data = (float *) malloc(ndst);
                    float *p = data;
                    int32_t np = no_samples * nval_per_sample[i];

                    for (size_t j=0; j<nfiles; ++j)
                    {
                        int32_t b = bcf_get_format_float(odrs[j]->hdr, vs[j], genotype_field, &p, &np);
                        p = p+bcf_hdr_nsamples(odrs[j]->hdr)*vs[j]->d.fmt[i].n;
                        np = np-bcf_hdr_nsamples(odrs[j]->hdr)*nval_per_sample[i];

                        if (vs[j]->d.fmt[i].n<nval_per_sample[i])
                        {
                            int32_t steps = bcf_hdr_nsamples(odrs[j]->hdr)*(nval_per_sample[i] - vs[j]->d.fmt[i].n);
                            while (steps)
                            {
                                p[0] = *(p-1);
                                ++p;
                                --steps;
                            }
                        }
                    }

                    bcf_update_format_float(odw->hdr, nv, genotype_field, data, no_samples * vs[0]->d.fmt[i].n);
                    free(data);
                }
                else if (vs[0]->d.fmt[i].type == BCF_BT_CHAR)
                {
                    int32_t ndst = no_samples * sizeof(char *);

                    char **data = (char **) malloc(ndst);
                    int32_t cp = 0;

                    for (size_t j=0; j<nfiles; ++j)
                    {
                        char** p = NULL;
                        int32_t np = 0;
                        int32_t b = bcf_get_format_char(odrs[j]->hdr, vs[j], genotype_field, &p, &np);
                        
                        for (size_t k=0; k<bcf_hdr_nsamples(odrs[j]->hdr); ++k)
                        {
                            data[cp+k] = p[k];
                            ++cp;    
                        }
                        free(p);
                    }

                    bcf_update_format_string(odw->hdr, nv, genotype_field, const_cast<const char**>(data), no_samples * vs[0]->d.fmt[i].n);
                    
                    cp = 0;
                    for (size_t j=0; j<nfiles; ++j)
                    {
                        free(data[cp]);
                        cp += bcf_hdr_nsamples(odrs[j]->hdr);
                    }                    
                    free(data);
                }
            }

            odw->write(nv);
        }

        for (size_t i=0; i<nfiles; ++i)
        {
            odrs[i]->close();
        }

        odw->close();
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "paste v" << version << "\n\n";
        print_ifiles("options:     input VCF file        ", input_vcf_files);
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        std::cerr << "stats: Total number of files pasted  " << input_vcf_files.size() << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool paste(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.paste();
    igor.print_stats();
    return igor.print;
}

