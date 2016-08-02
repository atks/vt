/* The MIT License

   Copyright (c) 2016 Hyun Min Kang <hmkang.umich.edu> and Adrian Tan <atks@umich.edu>

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

#include "genotype.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string output_vcf_file;
    std::vector<std::string> input_files;
    std::string input_files_list;
    std::vector<GenomeInterval> intervals;

    int32_t maxBQ;
    double minContam;

    std::string sex_map_file;
    std::map<std::string,int> mSex;
    std::string xLabel;
    std::string yLabel;
    std::string mtLabel;
    int32_t xStart;
    int32_t xStop;

    std::string sample_id;
    std::string ref_fasta_file;
    bool ignore_md;
    int32_t debug;

    /**
     * Joint         - genotypes multiple samples and outputs with features
     *                 genotypic information of samples are stored in memory,
     *                 this should be used with short genomic intervals and
     *                 maybe parallelized thus.
     * Hierarchical  - genotypes by subset of individuals and stores features
     *                 temporarily in INFO fields in VCF files.  Resultant
     *                 VCF files may be merged with vt merge_genotypes.  This 
     *                 should applied on subsets of individuals and by large 
     *                 genomic sequence segments - e.g. chromosome or whole 
     *                 genome.
     */
    std::string mode;

    //read filters
    uint32_t read_mapq_cutoff;
    uint16_t read_exclude_flag;
    bool ignore_overlapping_read;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;  // for anchor VCF
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;
    uint32_t no_vntrs_genotyped;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    Estimator *est;
    ReadFilter *rf;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Genotype SNPs, Indels, VNTRs across multiple samples and combines them\n"
                 "              Input requirements and assumptions:\n"
                 "              1. Same variants are represented in the same order for each file (required)\n"
                 "              2. Genotype fields are the same for corresponding records (required)\n"
                 "              3. Sample names are different in all the files (warning will be given if not)\n"
                 "              4. Headers (not including the samples) are the same for all the files (unchecked assumption, will fail if output is BCF)\n"
                 "              Outputs:\n"
                 "              1. INFO fields output will be that of the first file\n"
                 "              2. Genotype fields are the same for corresponding records\n";


            version = "0.1";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_files_list("L", "L", "file containing list of input SAM/BAM/CRAM/BCF/VCF files", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_files("<in1.vcf>...", "Multiple SAM/BAM/CRAM/BCF/VCF files",false, "files", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i","i","Intervals[]", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_sex_map_file("g", "sex-map", "file containing sex map of each individual. ID first, X ploidy second", false, "", "file", cmd);
            TCLAP::ValueArg<int32_t> arg_max_bq("q", "q", "Maximum base quality to cap []", false, 40, "int", cmd);
            TCLAP::ValueArg<double> arg_min_contam("c", "c", "Minimum contamiation value to floor []", false, 0.01, "double", cmd);
            //TCLAP::ValueArg<int32_t> arg_n_thread("t", "t", "Number of threads", false, 1, "int", cmd);
            TCLAP::ValueArg<std::string> arg_xLabel("", "xLabel", "Contig name for X chromosome", false, "X", "str", cmd);
            TCLAP::ValueArg<std::string> arg_yLabel("", "yLabel", "Contig name for Y chromosome", false, "Y", "str", cmd);
            TCLAP::ValueArg<std::string> arg_mtLabel("", "mtLabel", "Contig name for mitochondrial chromosome", false, "MT", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_xStart("", "xStart", "Start base position of non-PAR region in X chromosome", false, 2699520, "int", cmd);
            TCLAP::ValueArg<int32_t> arg_xStop("", "xStop", "End base position of non-PAR region in X chromosome", false, 154931044, "int", cmd);

            TCLAP::ValueArg<uint32_t> arg_read_mapq_cutoff("m", "m", "MAPQ cutoff for alignments (>=) [0]", false, 0, "int", cmd);
            TCLAP::SwitchArg arg_ignore_overlapping_read("l", "l", "ignore overlapping reads [false]", cmd, false);
            TCLAP::ValueArg<uint32_t> arg_read_exclude_flag("a", "a", "read exclude flag [0x0704]", false, 0x0704, "int", cmd);
            TCLAP::ValueArg<std::string> arg_mode("m", "m", "mode [d]\n"
                          "              d : parallel access across files (up to hundreds, less memory)\n"
                          "              s : sequential access across many files (thousands or more, with more memory)\n"
                          "              s : sequential access across many files (thousands or more, with more memory)\n",
                          false, "s", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference FASTA file []", true, "", "string", cmd);
            TCLAP::ValueArg<uint32_t> arg_debug("d", "d", "debug [0]", false, 0, "int", cmd);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);

            cmd.parse(argc, argv);
            ref_fasta_file = arg_ref_fasta_file.getValue();
            debug = arg_debug.getValue();

            read_mapq_cutoff = arg_read_mapq_cutoff.getValue();
            ignore_overlapping_read = arg_ignore_overlapping_read.getValue();
            read_exclude_flag = arg_read_exclude_flag.getValue();

            input_vcf_file = arg_input_vcf_file.getValue();
            input_files_list = arg_input_files_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_files(input_files, arg_input_files.getValue(), arg_input_files_list.getValue());
            const std::vector<std::string>& v = arg_input_files.getValue();

            mode = arg_mode.getValue();

            maxBQ = arg_max_bq.getValue();
            minContam = arg_min_contam.getValue();

            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            sex_map_file = arg_sex_map_file.getValue();
            xLabel = arg_xLabel.getValue();
            yLabel = arg_yLabel.getValue();
            mtLabel = arg_mtLabel.getValue();
            xStart = arg_xStart.getValue();
            xStop = arg_xStop.getValue();

            if (input_files.size()==0)
            {
                fprintf(stderr, "[E:%s:%d %s] no input files.\n", __FILE__, __LINE__, __FUNCTION__);
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
        ////////////////////////
        //Read sex map if needed
        ////////////////////////
        //todo:  process as a pedigree file?
        if ( !sex_map_file.empty() )
        {
            htsFile *file = hts_open(sex_map_file.c_str(),"r");

            if ( file == NULL )
            {
                fprintf(stderr,"ERROR: Cannot open %s\n",sex_map_file.c_str());
                exit(1);
            }
            kstring_t *s = &file->line;
            while( hts_getline(file,'\n',s) >= 0 )
            {
                std::string ss = std::string(s->s);
                size_t idx = ss.find_first_of("\t ");
                if ( idx == std::string::npos )
                {
                    fprintf(stderr,"ERROR: Cannot parse line %s in %s\n",ss.c_str(), sex_map_file.c_str());
                    exit(1);
                }
                std::string id = ss.substr(0, idx);
                int32_t sex = atoi(ss.substr(idx+1).c_str());

                if ( mSex.find(id) != mSex.end() )
                {
                    fprintf(stderr,"ERROR: Duplicate ID %s in %s\n",id.c_str(), sex_map_file.c_str());
                    exit(1);
                }

                if ( sex == 0 )
                {
                    fprintf(stderr,"WARNING: Unknown sex for individual %s, assuming females\n",id.c_str());
                    sex = 2;
                }
                else if ( sex > 2 )
                {
                    fprintf(stderr,"ERROR: Invalid sex %d for individual %s\n",sex,id.c_str());
                    exit(1);
                }

                mSex[id] = sex;
            }
        }

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_snps_genotyped = 0;
        no_indels_genotyped = 0;
        no_vntrs_genotyped = 0;

        /////////
        //tools//
        /////////
        rf = new ReadFilter(read_mapq_cutoff, read_exclude_flag, ignore_overlapping_read);
        vm = new VariantManip(ref_fasta_file);
        est = new Estimator();
    }

    /**
     * Genotypes individuals.
     */
    void genotype()
    {
        // assume that the following features are available
        // 1. BQSUM, DPF, DPR
        // 2. GT, PL, ADF, ADR, DP, CY, ST, AL, NM
        //
        // calculate the following per-variant statistics
        // DP - total depth -- \sum_i (DPF+DPR) + \sum_i DP
        // MLEAF - allele frequency estimated by best guess genotypes
        // HWEAF - allele frequency estimated by genotype likelihood under HWE
        // HWDAF - allele frequency estimated by genotype likelihood under non-HWD
        // SLRT - signed LRT test statistics for departure from HWE
        // IBC - [Pr(Het|HWD)-Pr(Het|HWE)]/Pr(Het|HWE)
        // AB - Allele balance -- \sum_i \sqrt{#A + #R} (#A)/(#A+#R) | GT = HET
        // STR - Strand bias   -- \sum_i \sqrt{#A + #R} cor(FR,FA,RR,RA)
        // TBR - Tail distance bias  -- \sum_i \sqrt{#A + #R} cor(TD,AL)
        // MQR - Mapping quality bias -- \sum_i \sqrt{#A + #R} cor(MQ,AL)
        // BQR - Base quality bias -- \sum_i \sqrt{#A + #R} cor(BQ,AL)
        // NMR - Number of mismatch bias -- \sum \sqrt{#A + #R} cor(NM,AL - A?) | BQ > 20
        // IOR - Inflated fraction of "other" biases compared to base quality -- \sum_i \sqrt{#A + #R} {Pr(O) - E[O]} | BQ > 20
        //
        // calculate the following per-sample statistics
        // GT - best guess genotypes
        // GQ - genotype quality
        // AD - allele depth
        // PL - likelihoods
        //
        // Assume that either
        // BQSUM : DPF : DPR or
        // GT, PL, DP, CY, ST, NM exists

        // classify input files by SAM/BAM/CRAM vs BCF/VCF


        // first, check the file format
        // format is delmited
        // 1. SAMPLE
        // 2. SAM/BAM/CRAM path
        // 3. contamination fraction estimate (optional)
        std::vector<std::string> input_vcf_files;
        std::vector<std::string> input_sam_files;
        std::vector<std::string> input_sam_sample_names;
        std::vector<double>      input_sam_contams;

        std::vector<std::string> tokens;

        for(int32_t i=0; i < (int)input_files.size(); ++i)
        {
            split(tokens, "\t\r\n ", input_files[i]);

            std::string fname;
            std::string sname;
            double contam = 0;

            if (tokens.size()!=2 || tokens.size()!=3)
            {
                error("Cannot parse input line %s, 2 or 3 columns expected", input_files[i].c_str());
            }

            sname = tokens[0];
            int32_t file_type = hts_filename_type(tokens[1]);
            if (file_type==sam || file_type==bam || file_type==cram)
            {
                error("SAM/CRAM/BAM file expected in second column: %s",input_files[i].c_str());
            }
            fname = tokens[1];
            if (tokens.size()==3)
            {
                if (!str2double(tokens[2], contam))
                {
                    error("number expected in third column: %s",input_files[i].c_str());
                }

                if ( contam < minContam )
                {
                    contam = minContam;
                }
            }

            input_sam_sample_names.push_back(sname);
            input_sam_files.push_back(fname);
            input_sam_contams.push_back(contam);
        }

        notice("Identified a total of %u VCF/BCF files and %u SAM/CRAM/CRAM files as input to integrate", input_vcf_files.size(), input_sam_files.size());

        if ( input_vcf_files.size() == 0 )
        {
            error("At least one BCF/VCF files are required as input");
        }

        if ( ( input_vcf_files.size() > 1 ) && ( input_sam_files.size() > 0 ) )
        {
            error("BCF/VCF and SAM/BAM/CRAM files cannot be merged at the same time");
        }

        if ( input_sam_files.empty() )
        { // BCF/VCF only
            abort();
        }
        else
        {
            if (mode=="d")
            {
                BAMOrderedReader odr(input_sam_files[0], intervals);
                bam_hdr_t *h = odr.hdr;
                bam1_t * s = bam_init1();
                BCFSingleGenotypingBufferedReader bsgr(input_vcf_file, intervals, output_vcf_file);

                while (odr.read(s))
                {
                    ++rf->no_reads;

                    if (!rf->filter_read(h, s))
                    {
                        continue;
                    }

                    bsgr.process_read(h, s);

                    ++rf->no_passed_reads;
                    if ((rf->no_reads & 0x0000FFFF) == 0)
                    {
                        std::cerr << bam_get_chrom(h,s) << ":" << bam_get_pos1(s) << " ("  << bsgr.buffer.size() << ")\n";
                    }
                }

                bsgr.flush(h, s, true);

                no_snps_genotyped = bsgr.no_snps_genotyped;
                no_indels_genotyped = bsgr.no_indels_genotyped;
                no_vntrs_genotyped = bsgr.no_vntrs_genotyped;

                odw->close();
            }
            else
            {
                int32_t nsamples = (int32_t)input_sam_files.size();

                notice("Loading input VCF file %s in region %s", input_vcf_files[0].c_str(), intervals[0].to_string().c_str());
                BCFGenotypingBufferedReader jgbr(input_vcf_files[0], intervals, output_vcf_file, nsamples);

                BCFOrderedWriter* odw = new BCFOrderedWriter(output_vcf_file);
                bcf_hdr_transfer_contigs(jgbr.odr->hdr, odw->hdr);

                bam1_t* s = bam_init1();

                std::vector<std::string> snames(nsamples);

                for(int32_t i=0; i < nsamples; ++i)
                {
                    BAMOrderedReader odr(input_sam_files[i], intervals);
                    bam_hdr_t* h = odr.hdr;
                    int64_t no_reads = 0;

                    snames[i] = input_sam_sample_names[i].empty() ? bam_hdr_get_sample_name(odr.hdr) : input_sam_sample_names[i];
                    jgbr.set_sample(i, snames[i].c_str(), input_sam_contams[i]);

                    if ( i % 100 == 0 )
                      notice("Processing %d-th sample %s", i+1, snames[i].c_str());

                    if ( jgbr.numVariants() > 0 )
                    {
                        while( odr.read(s) )
                        {
                            ++no_reads;
                            if ( !rf->filter_read(odr.hdr, s) ) continue;

                            //jgbr.flush( h, s );
                            jgbr.process_read(h, s, i);  // process the next reads
                        }
                    }

                    if ( no_reads == 0 )
                    {
                        notice("WARNING: No read found in %d-th sample %s", i+1, snames[i].c_str());
                    }

                    jgbr.flush_sample(i);
                    odr.close();
                }

                bam_destroy1(s);
                for(int32_t i=0; i < nsamples; ++i)
                {
                    bcf_hdr_add_sample(odw->hdr, snames[i].c_str());
                }

                bcf_hdr_add_sample(odw->hdr, NULL);

                int32_t nvariants = jgbr.numVariants();
                jgbr.write_header(odw);

                for(int32_t i=0; i < nvariants; ++i)
                {
                    bcf1_t* nv = jgbr.flush_variant(i, odw->hdr);
                    odw->write(nv);
                    bcf_destroy(nv);
                }

                odw->close();
                delete odw;
            }
        }
    }

    void print_options()
    {
        std::clog << "genotype v" << version << "\n\n";
        std::clog << "options:    # input Files                       " << input_files.size() << "\n";
        std::clog << "         [o] output VCF File                      " << output_vcf_file << "\n";
        std::clog << "         [r] reference FASTA File                 " << ref_fasta_file << "\n";
        std::clog << "         [z] ignore MD tags                       " << (ignore_md ? "true": "false") << "\n";
        //std::clog << "         [m] mode of genotyping                   " << mode << "\n";
        print_int_op("         [i] intervals                            ", intervals);
        std::clog << "\n";
        std::clog << "         [t] read mapping quality cutoff          " << read_mapq_cutoff << "\n";
        std::clog << "         [l] ignore overlapping read              " << (ignore_overlapping_read ? "true" : "false") << "\n";
        std::clog << "         [a] read flag filter                     " << std::showbase << std::hex << read_exclude_flag << std::dec << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "genotype v" << version << "\n\n";

        std::clog << "\n";
        std::clog << "stats: no. reads                    : " << rf->no_reads << "\n";
        std::clog << "       no. overlapping reads        : " << rf->no_overlapping_reads << "\n";
        std::clog << "       no. low mapq reads           : " << rf->no_low_mapq_reads << "\n";
        std::clog << "       no. passed reads             : " << rf->no_passed_reads << "\n";
        std::clog << "       no. exclude flag reads       : " << rf->no_exclude_flag_reads << "\n";
        std::clog << "\n";
        std::clog << "       no. unaligned cigars         : " << rf->no_unaligned_cigars << "\n";
        std::clog << "       no. malformed del cigars     : " << rf->no_malformed_del_cigars << "\n";
        std::clog << "       no. malformed ins cigars     : " << rf->no_malformed_ins_cigars << "\n";
        std::clog << "       no. salvageable ins cigars   : " << rf->no_salvageable_ins_cigars << "\n";
        std::clog << "\n";
        std::clog << "       no. SNPs genotyped           : " << no_snps_genotyped<< "\n";
        std::clog << "       no. Indels genotyped         : " << no_indels_genotyped << "\n";
        std::clog << "       no. VNTRs genotyped          : " << no_vntrs_genotyped << "\n";
        std::clog << "\n";
    }

    ~Igor()
    {
    }

    private:
};

}

void genotype(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.genotype();
    igor.print_stats();
}

