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

#include "joint_genotype_sequential.h"

#include <errno.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

namespace
{

/**
 * For detecting overlapping reads.
 */
typedef struct
{
  int32_t start1, end1;
} interval_t;

KHASH_MAP_INIT_STR(rdict, interval_t)

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> input_files;
    std::string input_files_list;
    std::string output_vcf_file;
    bool print;
    std::vector<GenomeInterval> intervals;

    int32_t maxBQ;
    double minContam;
    int32_t nThreads;

    std::string sex_map_file;
    std::string xLabel;
    std::string yLabel;
    std::string mtLabel;
    int32_t xStart;
    int32_t xStop;

    std::string sample_id;
    std::string ref_fasta_file;
  //std::string mode;
    bool ignore_md;
    int32_t debug;

    //variables for keeping track of chromosome
    std::string chrom; //current chromosome
    int32_t tid; // current sequence id in bam
    int32_t rid; // current sequence id in bcf

    //read filters
    uint32_t read_mapq_cutoff;
    uint16_t read_exclude_flag;
    bool ignore_overlapping_read;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;  // for anchor VCF
    BCFOrderedWriter *odw;
    Estimator *est;

    std::map<std::string,int> mSex;
    khash_t(rdict) *reads;

    ///////////////
    //general use//
    ///////////////

    /////////
    //stats//
    /////////

    uint32_t no_reads;
    uint32_t no_overlapping_reads;
    uint32_t no_passed_reads;
    uint32_t no_exclude_flag_reads;
    uint32_t no_low_mapq_reads;
    uint32_t no_unaligned_cigars;
    uint32_t no_malformed_del_cigars;
    uint32_t no_malformed_ins_cigars;
    uint32_t no_salvageable_ins_cigars;

    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;
    uint32_t no_vntrs_genotyped;

    /////////
    //tools//
    /////////
    VariantManip *vm;

    // The function assumes a particular list of features in the BCF files
    //

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
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
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
        /*            TCLAP::ValueArg<std::string> arg_mode("m", "m", "mode [d]\n"
                          "              p : parallel access across files (up to hundreds, less memory)\n"
                          "              s : sequential access across many files (thousands or more, with more memory)\n",
                          false, "s", "str", cmd);*/
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference FASTA file []", true, "", "string", cmd);
            TCLAP::ValueArg<uint32_t> arg_debug("d", "d", "debug [0]", false, 0, "int", cmd);

            cmd.parse(argc, argv);
            ref_fasta_file = arg_ref_fasta_file.getValue();
            debug = arg_debug.getValue();

            read_mapq_cutoff = arg_read_mapq_cutoff.getValue();
            ignore_overlapping_read = arg_ignore_overlapping_read.getValue();
            read_exclude_flag = arg_read_exclude_flag.getValue();


            //mode = arg_mode.getValue();
            input_files_list = arg_input_files_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_files(input_files, arg_input_files.getValue(), arg_input_files_list.getValue());
            const std::vector<std::string>& v = arg_input_files.getValue();

            print = arg_print.getValue();
        maxBQ = arg_max_bq.getValue();
        minContam = arg_min_contam.getValue();
        nThreads = 1; //arg_n_thread.getValue();

#ifdef _OPENMP
        omp_set_num_threads(nThreads);
#endif

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

        est = new Estimator();
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

    void initialize()
    {
        ///////////////
        //general use//
        ///////////////

        // Read sex map if needed
        if ( !sex_map_file.empty() ) {
          htsFile *file = hts_open(sex_map_file.c_str(),"r");
          if ( file == NULL ) {
            fprintf(stderr,"ERROR: Cannot open %s\n",sex_map_file.c_str());
            exit(1);
          }
          kstring_t *s = &file->line;
          while( hts_getline(file,'\n',s) >= 0 ) {
            std::string ss = std::string(s->s);
            size_t idx = ss.find_first_of("\t ");
            if ( idx == std::string::npos ) {
              fprintf(stderr,"ERROR: Cannot parse line %s in %s\n",ss.c_str(), sex_map_file.c_str());
              exit(1);
            }
            std::string id = ss.substr(0, idx);
            int32_t sex = atoi(ss.substr(idx+1).c_str());
        
            if ( mSex.find(id) != mSex.end() ) {
              fprintf(stderr,"ERROR: Duplicate ID %s in %s\n",id.c_str(), sex_map_file.c_str());
              exit(1);
            }
        
            if ( sex == 0 ) {
              fprintf(stderr,"WARNING: Unknown sex for individual %s, assuming females\n",id.c_str());
              sex = 2;
            }
            else if ( sex > 2 ) {
              fprintf(stderr,"ERROR: Invalid sex %d for individual %s\n",sex,id.c_str());
              exit(1);
            }
            mSex[id] = sex;
          }
        }

        ////////////////////////
        //stats initialization//
        ////////////////////////

        /////////
        //tools//
        /////////
        vm = new VariantManip(ref_fasta_file);

        no_reads = 0;
        no_overlapping_reads = 0;
        no_passed_reads = 0;
        no_exclude_flag_reads = 0;
        no_low_mapq_reads = 0;
        no_unaligned_cigars = 0;
        no_malformed_del_cigars = 0;
        no_malformed_ins_cigars = 0;
        no_salvageable_ins_cigars = 0;

        no_snps_genotyped = 0;
        no_indels_genotyped = 0;
        no_vntrs_genotyped = 0;

        //for tracking overlapping reads
        reads = kh_init(rdict);

        chrom = "";
        tid = -1;
        rid = -1;
    }

    /**
     * Filter reads.
     *
     * Returns true if read is failed.
     */
  bool filter_read(bam_hdr_t* h, bam1_t *s) {
    khiter_t k;
    int32_t ret;

    if (ignore_overlapping_read) {
      //this read is part of a mate pair on the same contig
      if (bam_get_mpos1(s) && (bam_get_tid(s)==bam_get_mtid(s))) {
    //first mate
    if (bam_get_mpos1(s)>bam_get_pos1(s)) {
      //overlapping
      if (bam_get_mpos1(s)<=(bam_get_pos1(s) + bam_get_l_qseq(s) - 1)) {
        //add read that has overlapping
        //duplicate the record and perform the stitching later
        char* qname = strdup(bam_get_qname(s));
        k = kh_put(rdict, reads, qname, &ret);
        if (!ret) {
          //already present
          free(qname);
        }
        kh_val(reads, k) = {bam_get_pos1(s), bam_get_pos1(s)+bam_get_l_qseq(s)-1};
      }
    }
    else {
      //check overlap
      if((k = kh_get(rdict, reads, bam_get_qname(s)))!=kh_end(reads)) {
        if (kh_exist(reads, k)) {
          free((char*)kh_key(reads, k));
          kh_del(rdict, reads, k);
          ++no_overlapping_reads;
        }
        //set this on to remove overlapping reads.
        return false;
      }
    }
      }
    }

    if(bam_get_flag(s) & read_exclude_flag) {
      //1. unmapped
      //2. secondary alignment
      //3. not passing QC
      //4. PCR or optical duplicate
      ++no_exclude_flag_reads;
      return false;
    }

    if (bam_get_mapq(s) < read_mapq_cutoff) {
      //filter short aligments and those with too many indels (?)
      ++no_low_mapq_reads;
      return false;
    }

    //*****************************************************************
    //should we have an assertion on the correctness of the bam record?
    //Is, Ds not sandwiched in M
    //leading and trailing Is - convert to S
    //no Ms!!!!!
    //*****************************************************************
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    if (n_cigar_op) {
      uint32_t *cigar = bam_get_cigar(s);
      bool seenM = false;
      int32_t last_opchr = '^';

      for (int32_t i = 0; i < n_cigar_op; ++i) {
    int32_t opchr = bam_cigar_opchr(cigar[i]);
    int32_t oplen = bam_cigar_oplen(cigar[i]);
    if (opchr=='S') {
      if (i!=0 && i!=n_cigar_op-1) {
        std::cerr << "S issue\n";
        bam_print_key_values(h, s);
        //++malformed_cigar;
      }
    }
    else if (opchr=='M') {
      seenM = true;
    }
    else if (opchr=='D') {
      if (last_opchr!='M' || (i<=n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M')) {
        std::cerr << "D issue\n";
        ++no_malformed_del_cigars;
        bam_print_key_values(h, s);
      }
    }
    else if (opchr=='I') {
      if (last_opchr!='M' || (i<n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M')) {
        if (last_opchr!='M') {
          if (last_opchr!='^' && last_opchr!='S') {
        std::cerr << "leading I issue\n";
        bam_print_key_values(h, s);
        ++no_malformed_ins_cigars;
          }
          else {
        ++no_salvageable_ins_cigars;
          }
        }
        else if (i==n_cigar_op-1) {
          ++no_salvageable_ins_cigars;
        }
        else if (i==n_cigar_op-2 && (bam_cigar_opchr(cigar[i+1])=='S')) {
          ++no_salvageable_ins_cigars;
        }
        else {
          std::cerr << "trailing I issue\n";
          bam_print_key_values(h, s);
          ++no_malformed_ins_cigars;
        }
      }
    }

    last_opchr = opchr;
      }

      if (!seenM) {
    std::cerr << "NO! M issue\n";
    bam_print_key_values(h, s);
    ++no_unaligned_cigars;
      }
    }

    //check to see that hash should be cleared when encountering new contig.
    //some bams may not be properly formed and contain orphaned sequences that
    //can be retained in the hash
    if (bam_get_tid(s)!=tid) {
      for (k = kh_begin(reads); k != kh_end(reads); ++k) {
    if (kh_exist(reads, k)) {
      free((char*)kh_key(reads, k));
      kh_del(rdict, reads, k);
    }
      }

      tid = bam_get_tid(s);
    }

    return true;
  }

  inline bool ends_with(std::string const & value, std::string const & ending)
  {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
  }

  int32_t file_type(std::string const& value) {
    if ( ends_with(value,".vcf") || ends_with(value,".bcf") || ends_with(value,".vcf.gz") ) {
      return 2; // BCF/VCF
    }
    else if ( ends_with(value,".bam") || ends_with(value,".cram") || ends_with(value,".sam") ) {
      return 1;
    }
    else {
      errno = 0;
      char* endptr = NULL;
      double converted = strtod(value.c_str(), &endptr);

      if ( ( *endptr == 0 ) && ( errno == 0 ) ) {
    return 0;
      }
      else
    return -1;
    }
  }

  void joint_genotype_sequential()
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
    std::vector<std::string> input_vcf_files;
    std::vector<std::string> input_sam_files;
    std::vector<std::string> input_sam_sample_names;
    std::vector<double>      input_sam_contams;

    std::vector<std::string> tokens;

    for(int32_t i=0; i < (int)input_files.size(); ++i) {

      split(tokens, "\t\r\n ", input_files[i]);

      std::string fname;
      std::string sname;
      double contam = minContam;

      if ( tokens.size() == 1 ) {
    int32_t t0 = file_type(tokens[0]);
    if ( t0 > 0 ) {
      fname = tokens[0];
    }
    else {
      error("Cannot parse input line %s",input_files[i].c_str());
    }
      }
      else if ( tokens.size() == 2 ) {
    int32_t t0 = file_type(tokens[0]);
    int32_t t1 = file_type(tokens[1]);

    if ( ( t0 == 1 ) && ( t1 == 0 ) ) { // SAM - NUMERIC
      fname = tokens[0];
      contam = ::atof(tokens[1].c_str());
    }
    else if ( t1 == 1 ) {  // ID - SAM
      sname = tokens[0];
      fname = tokens[1];
    }
    else
      error("Cannot parse input line %s",input_files[i].c_str());
      }
      else if ( tokens.size() == 3 ) {
    int32_t t0 = file_type(tokens[0]);
    int32_t t1 = file_type(tokens[1]);
    int32_t t2 = file_type(tokens[2]);

    if ( ( t1 == 1 ) && ( t2 == 0 ) ) { // ID - SAM - NUMERIC
      sname = tokens[0];
      fname = tokens[1];
      contam = ::atof(tokens[2].c_str());
    }
    else
      error("Cannot parse input line %s",input_files[i].c_str());
      }

      if ( contam < minContam )
    contam = minContam;

      if ( file_type(fname) == 2 ) {
    input_vcf_files.push_back(fname);
      }
      else {
    input_sam_sample_names.push_back(sname);
    input_sam_files.push_back(fname);
    input_sam_contams.push_back(contam);
      }
    }

    notice("Identified a total of %u VCF/BCF files and %u SAM/CRAM/CRAM files as input to integrate",input_vcf_files.size(), input_sam_files.size());

    if ( input_vcf_files.size() == 0 )
      error("At least one BCF/VCF files are required as input");
    if ( ( input_vcf_files.size() > 1 ) && ( input_sam_files.size() > 0 ) )
      error("BCF/VCF and SAM/BAM/CRAM files cannot be merged at the same time");

    if ( input_sam_files.empty() ) { // BCF/VCF only
      abort();
    }
    else {
      int32_t nsamples = (int32_t)input_sam_files.size();

      notice("Loading input VCF file %s in region %s", input_vcf_files[0].c_str(), intervals[0].to_string().c_str());
      JointGenotypingBufferedReader jgbr(input_vcf_files[0], intervals, output_vcf_file, nsamples);

      BCFOrderedWriter* odw = new BCFOrderedWriter(output_vcf_file);
      bcf_hdr_transfer_contigs(jgbr.odr->hdr, odw->hdr);

      bam1_t* s = bam_init1();

#ifdef _OPENMP
      notice("Running parallel jobs across %d CPUs", nThreads);
#endif

      std::vector<std::string> snames(nsamples);

#pragma omp parallel for schedule(dynamic,1)
      for(int32_t i=0; i < nsamples; ++i) {
    BAMOrderedReader odr(input_sam_files[i], intervals);
    bam_hdr_t* h = odr.hdr;
    int64_t no_reads = 0;

    snames[i] = input_sam_sample_names[i].empty() ? bam_hdr_get_sample_name(odr.hdr) : input_sam_sample_names[i];
    jgbr.set_sample(i, snames[i].c_str(), input_sam_contams[i]);

    if ( i % 100 == 0 )
      notice("Processing %d-th sample %s", i+1, snames[i].c_str());

    if ( jgbr.numVariants() > 0 ) {
      while( odr.read(s) ) {
        ++no_reads;
        if ( !filter_read(odr.hdr, s) ) continue;

        //jgbr.flush( h, s );
        jgbr.process_read(h, s, i);  // process the next reads
      }
    }

    if ( no_reads == 0 ) {
      notice("WARNING: No read found in %d-th sample %s", i+1, snames[i].c_str());
    }

    jgbr.flush_sample(i);
    odr.close();
      }

      bam_destroy1(s);
      for(int32_t i=0; i < nsamples; ++i) {
    bcf_hdr_add_sample(odw->hdr, snames[i].c_str());
      }

      bcf_hdr_add_sample(odw->hdr, NULL);

      int32_t nvariants = jgbr.numVariants();
      jgbr.write_header(odw);

#pragma omp parallel for ordered schedule(static,1)
      for(int32_t i=0; i < nvariants; ++i) {
    bcf1_t* nv = jgbr.flush_variant(i, odw->hdr);

#pragma omp ordered
    {
      odw->write(nv);
    }
    bcf_destroy(nv);
      }
      odw->close();
      delete odw;
    }
  }

  /**
   * Print BAM for debugging purposes.
   */
  void bam_print_key_values(bam_hdr_t *h, bam1_t *s)
  {
        const char* chrom = bam_get_chrom(h, s);
        uint32_t pos1 = bam_get_pos1(s);
        kstring_t seq = {0,0,0};
        bam_get_seq_string(s, &seq);
        uint32_t len = bam_get_l_qseq(s);
        kstring_t qual = {0,0,0};
        bam_get_qual_string(s, &qual);
        kstring_t cigar_string = {0,0,0};
        bam_get_cigar_string(s, &cigar_string);
        kstring_t cigar_expanded_string = {0,0,0};
        bam_get_cigar_expanded_string(s, &cigar_expanded_string);
        uint16_t flag = bam_get_flag(s);
        uint32_t mapq = bam_get_mapq(s);

        uint8_t *aux;
        char* md = NULL;
        (aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(aux));

        std::cerr << "##################" << "\n";
        std::cerr << "chrom:pos: " << chrom << ":" << pos1 << "\n";
        std::cerr << "read     : " << seq.s << "\n";
        std::cerr << "qual     : " << qual.s << "\n";
        std::cerr << "cigar_str: " << cigar_string.s << "\n";
        std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
        std::cerr << "len      : " << len << "\n";
        std::cerr << "mapq     : " << mapq << "\n";
        std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
        std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";
        std::cerr << "md       : " << (aux?md:"") << "\n";
        std::cerr << "##################" << "\n";

        if (seq.m) free(seq.s);
        if (qual.m) free(qual.s);
        if (cigar_string.m) free(cigar_string.s);
        if (cigar_expanded_string.m) free(cigar_expanded_string.s);
  }

  void print_options() {
        std::clog << "genotype2 v" << version << "\n\n";
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
        if (!print) return;

        std::clog << "\n";
        std::cerr << "stats: Total number of files pasted  " << input_files.size() << "\n";
        std::clog << "\n";
    }

    ~Igor() {
        kh_destroy(rdict, reads);
    }

    private:
};

}

bool joint_genotype_sequential(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.joint_genotype_sequential();
    igor.print_stats();
    return igor.print;
}

