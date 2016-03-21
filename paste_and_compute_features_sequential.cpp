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

#include "paste_and_compute_features_sequential.h"

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
    std::vector<GenomeInterval> intervals;
    int32_t maxBQ;
    std::string sex_map_file;  
    std::string xLabel;
    std::string yLabel;
    std::string mtLabel;
    int32_t xStart;
    int32_t xStop;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;
    Estimator *est;

    std::map<std::string,int> mSex;

    ///////////////
    //general use//
    ///////////////

    /////////
    //stats//
    /////////

    /////////
    //tools//
    /////////

    // The function assumes a particular list of features in the BCF files
    //

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


            version = "0.1";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);
	    TCLAP::ValueArg<std::string> arg_intervals("i","i","Intervals[]", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_sex_map_file("g", "sex-map", "file containing sex map of each individual. ID fitst, X ploidy second", false, "", "file", cmd);	    
            TCLAP::ValueArg<int32_t> arg_max_bq("q", "q", "Maximum base quality to cap []", false, 30, "int", cmd);
            TCLAP::ValueArg<std::string> arg_xLabel("", "xLabel", "Contig name for X chromosome", false, "X", "str", cmd);
            TCLAP::ValueArg<std::string> arg_yLabel("", "yLabel", "Contig name for Y chromosome", false, "Y", "str", cmd);
            TCLAP::ValueArg<std::string> arg_mtLabel("", "mtLabel", "Contig name for mitochondrial chromosome", false, "MT", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_xStart("", "xStart", "Start base position of non-PAR region in X chromosome", false, 2699520, "int", cmd);
            TCLAP::ValueArg<int32_t> arg_xStop("", "xStop", "End base position of non-PAR region in X chromosome", false, 154931044, "int", cmd); 	    
	
            cmd.parse(argc, argv);

            input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            const std::vector<std::string>& v = arg_input_vcf_files.getValue();
            print = arg_print.getValue();
	    maxBQ = arg_max_bq.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
	    sex_map_file = arg_sex_map_file.getValue();
	    xLabel = arg_xLabel.getValue();
	    yLabel = arg_yLabel.getValue();
	    mtLabel = arg_mtLabel.getValue();	    
	    xStart = arg_xStart.getValue();
	    xStop = arg_xStop.getValue();

            if (input_vcf_files.size()==0)
            {
                fprintf(stderr, "[E:%s:%d %s] no input vcf files.\n", __FILE__, __LINE__, __FUNCTION__);
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
        //////////////////////
        //i/o initialization//
        //////////////////////
        // read only the first file
        odr = new BCFOrderedReader(input_vcf_files[0], intervals);
	if ( intervals.size() > 1 ) {
	  fprintf(stderr, "[E:%s:%d %s] Multiple intervals are not allowed\n", __FILE__, __LINE__, __FUNCTION__);
	  exit(1);	  
	}
      
        //for (size_t i=0; i<input_vcf_files.size(); ++i)
        //{
        //    odrs.push_back(new BCFOrderedReader(input_vcf_files[i], intervals));
        //}
        odw = new BCFOrderedWriter(output_vcf_file, 0);
	
        odw->set_hdr(odr->hdr);

        bcf_hdr_append(odw->hdr, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average Depth per Sample\">\n");	
        bcf_hdr_append(odw->hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">\n");
        //bcf_hdr_append(odw->hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Reads\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency from Best-guess Genotypes\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=GN,Number=1,Type=Integer,Description=\"Total Number of Genotypes\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWEAF,Number=A,Type=Float,Description=\"Genotype likelihood based Allele Frequency assuming HWE\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWDGF,Number=G,Type=Float,Description=\"Genotype likelihood based Genotype Frequency ignoring HWE\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=IBC,Number=1,Type=Float,Description=\"Inbreeding Coefficients calculated from genotype likelihoods\">\n");	
        bcf_hdr_append(odw->hdr, "##INFO=<ID=HWE_SLP,Number=1,Type=Float,Description=\"Signed log p-values testing  statistics based Hardy Weinberg ln(Likelihood Ratio)\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=ABE,Number=1,Type=Float,Description=\"Expected allele Balance towards Reference Allele on Heterozygous Sites\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=ABZ,Number=1,Type=Float,Description=\"Average Z-scores of Allele Balance towards Reference Allele on Heterozygous Sites\">\n");	
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NS_NREF,Number=1,Type=Integer,Description=\"Number of samples with non-reference reads\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=BQZ,Number=1,Type=Float,Description=\"Correlation between base quality and alleles\">\n");
        //bcf_hdr_append(odw->hdr, "##INFO=<ID=MQZ,Number=1,Type=Float,Description=\"Correlation between mapping quality and alleles\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=CYZ,Number=1,Type=Float,Description=\"Correlation between cycle and alleles\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=STZ,Number=1,Type=Float,Description=\"Correlation between strand and alleles\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NMZ,Number=1,Type=Float,Description=\"Correlation between mismatch counts per read and alleles\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=IOR,Number=1,Type=Float,Description=\"Inflated rate of observing of other alleles in log10 scale\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NM0,Number=1,Type=Float,Description=\"Average number of mismatches in the reads with ref alleles\">\n");	
        bcf_hdr_append(odw->hdr, "##INFO=<ID=NM1,Number=1,Type=Float,Description=\"Average number of mismatches in the reads with non-ref alleles\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depth\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=OD,Number=1,Type=Integer,Description=\"Other Allele Depth\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scale Genotype Likelihoods\">\n");	

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
    }


    float compute_correlation(int32_t n, int32_t xy, int32_t x1, int32_t x2, int32_t y1, int32_t y2, float buffer) {
      if ( n == 0 ) return 0;
      float xsd = x2/(float)n - (x1/(float)n)*(x1/(float)n);
      float ysd = y2/(float)n - (y1/(float)n)*(y1/(float)n);
      return ( ( xy/(float)n - x1 * y1 / (float) n / (float) n ) / sqrt( xsd * ysd + buffer ) );
    }

    float compute_correlation_f(int32_t n, float xy, float x1, float x2, float y1, float y2, float buffer) {
      if ( n == 0 ) return 0;
      float xsd = x2/(float)n - (x1/(float)n)*(x1/(float)n);
      float ysd = y2/(float)n - (y1/(float)n)*(y1/(float)n);
      return ( ( xy/(float)n - x1 * y1 / (float) n / (float) n ) / sqrt( xsd * ysd + buffer ) );
    }

    void paste_and_compute_features_sequential()
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

	// store everything into memory
	std::vector<int32_t> vn_alleles;
	std::vector<int32_t> vn_genos;
	std::vector<int32_t*> v_pls;
	std::vector<int32_t*> v_gts;
	std::vector<int32_t*> v_gqs;
	std::vector<int32_t*> v_ads;
	std::vector<int32_t*> v_ods;

	// check the type from FORMAT, first field must be either BQSUM or GT
	int32_t* p_bqsum = NULL;
	int32_t np_bqsum = 0;
	int32_t* p_dp = NULL;
	int32_t np_dp = 0;
	int32_t* p_gt = NULL;
	int32_t np_gt = 0;
	int32_t* p_pl = NULL;
	int32_t np_pl = 0;
	int32_t* p_bq = NULL;
	int32_t np_bq = 0;
	int32_t* p_mq = NULL;
	int32_t np_mq = 0;
	int32_t* p_cy = NULL;
	int32_t np_cy = 0;
	char**   p_st = NULL;
	int32_t np_st = 0;
	int32_t*  p_al = NULL;
	int32_t np_al = 0;
	int32_t* p_nm = NULL;
	int32_t np_nm = 0;

	// read the first genotype file to populate
	int32_t nfiles = input_vcf_files.size();
	std::vector<bool> v_skips;
	std::vector<bool> v_is_snps;	
	std::vector<int32_t> v_rids;
	std::vector<int32_t> v_poss;
	std::vector<int32_t> v_rlens;
	std::vector< std::vector<std::string> > v_d_alleles;
	std::vector< std::vector<int32_t> > v_filts;

	// determine rids for X chromosome
	int32_t x_rid = bcf_hdr_name2id(odr->hdr,xLabel.c_str());
	bool x_found = false;

        bcf1_t* v = bcf_init();
	while( odr->read(v) ) {
	  // skip multi-allelics
	  bool skip = false;
	  bcf_unpack(v, BCF_UN_ALL);
	  
	  if ( v->n_allele > 2 ) skip = true;
	  else if ( ( !intervals.empty() ) && ( ( v->pos < intervals[0].start1 ) || ( v->pos > intervals[0].end1 ) ) ) skip = true;
	  else {
	    bool is_vntr = false;
	    for(size_t i=0; i < v->n_allele; ++i) {
	      if ( strcmp(v->d.allele[i],"<VNTR>") == 0 )
		is_vntr = true;
	    }
	    if ( is_vntr ) skip = true;
	  }

	  // determine whether to skip the marker or not
	  v_skips.push_back(skip); 	  
	  if ( skip ) continue;	  

	  // populate marker information
	  v_is_snps.push_back(bcf_is_snp(v));
	  v_rids.push_back(v->rid);
	  v_poss.push_back(v->pos);
	  v_rlens.push_back(v->rlen);
	  v_d_alleles.push_back(std::vector<std::string>());
	  for(size_t i=0; i < v->n_allele; ++i) {
	    v_d_alleles.back().push_back(v->d.allele[i]);
	  }
	  vn_alleles.push_back(v->n_allele);
	  vn_genos.push_back(v->n_allele * (v->n_allele+1)/2);

	  if ( ( x_rid == v->rid ) && ( v->pos >= xStart ) && ( v->pos <= xStop ) )
	    x_found = true;

	  v_filts.push_back(std::vector<int32_t>());
	  for(size_t i=0; i < v->d.n_flt; ++i) {
	    v_filts.back().push_back(v->d.flt[i]);
	  }
	  
	  // initialize genotype fields
	  v_pls.push_back( (int32_t*)calloc( nfiles * vn_genos.back(), sizeof(int32_t) ) );
	  v_ads.push_back( (int32_t*)calloc( nfiles * vn_alleles.back(), sizeof(int32_t) ) );
	  v_gts.push_back( (int32_t*)calloc( nfiles * 2, sizeof(int32_t) ) );
	  v_gqs.push_back( (int32_t*)calloc( nfiles, sizeof(int32_t) ) );
	  v_ods.push_back( (int32_t*)calloc( nfiles, sizeof(int32_t) ) );
	}
	odr->close();
	delete odr;

	std::vector<float> v_bqr_nums(v_rids.size(), 0);
	std::vector<float> v_bqr_dens(v_rids.size(), 0);
	std::vector<float> v_mqr_nums(v_rids.size(), 0);
	std::vector<float> v_mqr_dens(v_rids.size(), 0);
	std::vector<float> v_cyr_nums(v_rids.size(), 0);
	std::vector<float> v_cyr_dens(v_rids.size(), 0);
	std::vector<float> v_str_nums(v_rids.size(), 0);
	std::vector<float> v_str_dens(v_rids.size(), 0);
	std::vector<float> v_nmr_nums(v_rids.size(), 0);
	std::vector<float> v_nmr_dens(v_rids.size(), 0);			
	std::vector<float> v_ior_nums(v_rids.size(), 0);
	std::vector<float> v_ior_dens(v_rids.size(), 0);
	std::vector<float> v_nm0_nums(v_rids.size(), 0);
	std::vector<float> v_nm0_dens(v_rids.size(), 0);
	std::vector<float> v_nm1_nums(v_rids.size(), 0);
	std::vector<float> v_nm1_dens(v_rids.size(), 0);
	std::vector<float> v_ab_nums(v_rids.size(), 0);
	std::vector<float> v_ab_dens(v_rids.size(), 0);
	std::vector<float> v_abz_nums(v_rids.size(), 0);
	std::vector<float> v_abz_dens(v_rids.size(), 0);	
	std::vector<int32_t> v_ns_nrefs(v_rids.size(), 0);
	std::vector<int32_t> v_dp_sums(v_rids.size(), 0);
	std::vector<int32_t> v_max_gqs(v_rids.size(), 0);

	//bcf_clear(v);
	bcf_destroy(v);

	for(size_t i=0; i < nfiles; ++i) {
	  //fprintf(stderr,"Reading input file %s..\n", input_vcf_files[i].c_str());
	  // set odr and v
	  odr = new BCFOrderedReader(input_vcf_files[i], intervals);
	  v = bcf_init();

	  if ( bcf_hdr_nsamples(odr->hdr) != 1 ) {
	      fprintf(stderr, "[E:%s:%d %s] The genotype file must contain exactly one sample", __FILE__, __LINE__, __FUNCTION__);
	      exit(1); 	    
	  }
	  if ( i > 0 ) 
	    bcf_hdr_add_sample(odw->hdr, bcf_hdr_get_sample_name(odr->hdr, 0));	  

	  for( size_t j=0, k=0; j < v_skips.size(); ++j) {
	    if ( ! odr->read(v) ) {
	      fprintf(stderr, "[E:%s:%d %s] Cannot read variant from genotype files. j=%zu, k=%zu, pos[k]=%d", __FILE__, __LINE__, __FUNCTION__, j, k, v_poss[k]);
	      exit(1);		
	    }
	    if ( v_skips[j] ) continue;
	    bcf_unpack(v, BCF_UN_ALL);	    

	    // check marker infor with anchor files
	    if ( ( v->rid != v_rids[k] ) || ( v->pos != v_poss[k] ) || ( v->rlen != v_rlens[k] ) || ( v->n_allele != vn_alleles[k] ) ) {
	      fprintf(stderr, "[E:%s:%d %s] Variant position or ref alleles does not match\n", __FILE__, __LINE__, __FUNCTION__);
	      exit(1);	      
	    }

	    for(size_t l=0; l < vn_alleles[k]; ++l) {
	      if ( v_d_alleles[k][l].compare(v->d.allele[l]) ) {
		  fprintf(stderr, "[E:%s:%d %s] Variant alleles does not match\n", __FILE__, __LINE__, __FUNCTION__);
		  exit(1);		
	      }
	    }

	    // extract genotype fields and calculate summary statistics
	    if ( bcf_get_format_int32(odr->hdr, v, "BQSUM", &p_bqsum, &np_bqsum) >= 0 ) { // BQSUM observed - REF-ONLY
	      if ( bcf_get_format_int32(odr->hdr, v, "DP", &p_dp, &np_dp) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- BQSUM, DP\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		
	      }
	      if ( ( np_bqsum != np_dp ) || ( np_dp != 1 ) ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not have the samme number of fields -- BQSUM, DP\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		
	      }

	      for(size_t l=0; l < vn_alleles[k]; ++l){
		v_ads[k][ i*vn_alleles[k] + l] = (l == 0 ? p_dp[0] : 0);
		v_ods[k][ i ] = 0;
		for(size_t m=0; m <= l; ++m)  {
		  if ( m == 0 ) {
		    if ( l == 0 )
		      v_pls[k][vn_genos[k] * i + l * (l+1) / 2 + m] = 0;
		    else
		      v_pls[k][vn_genos[k] * i + l * (l+1) / 2 + m] = (int32_t)floor(3.0103 * p_dp[0] + 0.5);
		  }
		  else
		    v_pls[k][vn_genos[k] * i + l * (l+1) / 2 + m] = (int32_t)floor(p_bqsum[0] + 4.771213 * p_dp[0] + 0.5);
		  }
	      }
	      v_dp_sums[k] += p_dp[0];
	    }
	    else if ( bcf_get_genotypes(odr->hdr, v, &p_gt, &np_gt) >= 0 ) {  // GT unobserved -- non-REF
	      // extract PL, DP, BQ, MQ, CY, ST, AL, NM
	      if ( bcf_get_format_int32(odr->hdr, v, "PL", &p_pl, &np_pl) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- PL\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }		
	      if ( bcf_get_format_int32(odr->hdr, v, "DP", &p_dp, &np_dp) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- DP\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }
	      if ( bcf_get_format_int32(odr->hdr, v, "BQ", &p_bq, &np_bq) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- BQ\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }
	      if ( bcf_get_format_int32(odr->hdr, v, "MQ", &p_mq, &np_mq) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- MQ\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }
	      if ( bcf_get_format_int32(odr->hdr, v, "CY", &p_cy, &np_cy) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- GT, CY\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }
	      if ( bcf_get_format_string(odr->hdr, v, "ST", &p_st, &np_st) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- GT, ST\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }		
	      if ( bcf_get_format_int32(odr->hdr, v, "AL", &p_al, &np_al) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- GT, AL\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }
	      if ( bcf_get_format_int32(odr->hdr, v, "NM", &p_nm, &np_nm) < 0 ) {
		fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- GT, NM\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  
	      }

	      // sanity checking
	      if ( np_pl != vn_genos[k] ) {
		fprintf(stderr, "[E:%s:%d %s] np_pl (%d) != n_genos (%d)\n", __FILE__, __LINE__, __FUNCTION__, np_pl, vn_genos[k]);
		exit(1);
	      }

	      if ( ( np_dp != 1 ) || ( np_gt != 2 ) ) {
		fprintf(stderr, "[E:%s:%d %s] Assertion failed in np_dp (%d) == np_gt (%d),  np_st (%d) == np_al (%d)\n", __FILE__, __LINE__, __FUNCTION__, np_dp, np_gt, np_st, np_al);
		exit(1);		
	      }
	      if ( p_dp[0] != strlen(p_st[0]) ) {
		fprintf(stderr, "[E:%s:%d %s] Assetion failed - p_dp[0] == strlen(p_st[0])\n", __FILE__, __LINE__, __FUNCTION__);
		exit(1);		  		  		  		  		
	      }

	      // currently assume everything as diploid
	      int32_t a1 = bcf_gt_allele(p_gt[0]);
	      int32_t a2 = bcf_gt_allele(p_gt[1]);
	      int32_t gt = bcf_alleles2gt(a1, a2);
	      for(size_t l=0; l < vn_genos[k]; ++l) {
		if ( p_pl[l] > 255 )
		  p_pl[l] = 255;
		v_pls[k][vn_genos[k] * i + l] = p_pl[l];
	      }
	      
	      v_dp_sums[k] += p_dp[0];

	      int32_t bq_s1 = 0;
	      int32_t bq_s2 = 0;
	      int32_t mq_s1 = 0;
	      int32_t mq_s2 = 0;
	      float   cy_s1 = 0;
	      float   cy_s2 = 0;
	      int32_t st_s1 = 0;
	      int32_t al_s1 = 0;
	      int32_t nm_s1 = 0;
	      int32_t nm_s2 = 0;
	      int32_t bq_al = 0;
	      int32_t mq_al = 0;
	      float   cy_al = 0;
	      int32_t st_al = 0;
	      int32_t nm_al = 0;
	      int32_t dp_ra = 0;
	      int32_t dp_q20 = 0;
	      double  oth_exp_q20 = 0;
	      double  oth_obs_q20 = 0;
	      
	      int32_t* ads_z = &v_ads[k][ i*vn_alleles[k] ];

	      ++v_ns_nrefs[k];
	      
	      for(size_t l=0; l < p_dp[0]; ++l) {
		if ( p_bq[l] > maxBQ ) p_bq[l] = maxBQ;
		
		if ( p_al[l] >= 0 ) {
		  if ( p_bq[l] > 20 ) {
		    oth_exp_q20 += ( est->lt->pl2prob(p_bq[l]) * 2.0 / 3.0 );
		    ++dp_q20;
		  }
		  
		  // calculate cycle-based tail distance
		  float log_td = 0-logf((float)abs(v_is_snps[k] ? p_cy[l] : (int)(rand() % 100))+1.); // temporarily ignore cycles
		  
		  ++dp_ra;
		  bq_s1 += p_bq[l];
		  bq_s2 += (p_bq[l] * p_bq[l]);
		  mq_s1 += p_mq[l];
		  mq_s2 += (p_mq[l] * p_mq[l]);
		  cy_s1 += log_td;
		  cy_s2 += (log_td * log_td);
		  st_s1 += ((p_st[0][l] == 'F') ? 0 : 1);
		  if ( p_al[l] > 0 ) {
		    ++al_s1;
		    bq_al += p_bq[l];
		    mq_al += p_mq[l];
		    cy_al += log_td; 
		    st_al += ((p_st[0][l] == 'F') ? 0 : 1);
		    nm_al += (p_nm[l]-1);
		    nm_s1 += (p_nm[l]-1);
		    nm_s2 += ((p_nm[l]-1)*(p_nm[l]-1));
		    ++ads_z[1];
		  }
		  else {
		    nm_s1 += p_nm[l];
		    nm_s2 += ( p_nm[l] * p_nm[l] );
		    ++ads_z[0];
		  }
		}
		else {
		  if ( p_bq[l] > 20 ) {
		    oth_exp_q20 += ( est->lt->pl2prob(p_bq[l]) * 2.0 / 3.0 );
		    ++oth_obs_q20;
		    ++dp_q20;
		  }
		  ++v_ods[k][i];
		}
	      }
	      
	      float sqrt_dp_ra = sqrt((float)dp_ra);
	      float ior = (float)(oth_obs_q20 / (oth_exp_q20 + 1e-6));
	      float nm1 = al_s1 == 0 ? 0 : nm_al / (float)al_s1;
	      float nm0 = (dp_ra-al_s1) == 0 ? 0 : (nm_s1-nm_al) / (float)(dp_ra-al_s1);		
	      float w_dp_ra  = log(dp_ra+1.); //sqrt(dp_ra);
	      float w_dp_q20 = log(dp_q20+1.); //sqrt(dp_q20);
	      float w_al_s1  = log(al_s1+1.); //sqrt(al_s1);
	      float w_ref_s1 = log(dp_ra-al_s1+1.);
	      
	      if ( gt == 1 ) { // het genotypes
		v_ab_nums[k] += (w_dp_ra * (dp_ra - al_s1 + 0.05) / (double)(dp_ra + 0.1));
		v_ab_dens[k] += w_dp_ra;
		
		// E(r) = 0.5(r+a) V(r) = 0.25(r+a)
		v_abz_nums[k] += w_dp_ra * (dp_ra - al_s1 - dp_ra*0.5)/sqrt(0.25 * dp_ra + 1e-3);
		v_abz_dens[k] += (w_dp_ra * w_dp_ra);
		
		float bqr = sqrt_dp_ra * compute_correlation( dp_ra, bq_al, bq_s1, bq_s2, al_s1, al_s1, .1 );
		float mqr = sqrt_dp_ra * compute_correlation( dp_ra, mq_al, mq_s1, mq_s2, al_s1, al_s1, .1 );
		float cyr = sqrt_dp_ra * compute_correlation_f( dp_ra, cy_al, cy_s1, cy_s2, (float)al_s1, (float)al_s1, .1 );
		float str = sqrt_dp_ra * compute_correlation( dp_ra, st_al, st_s1, st_s1, al_s1, al_s1, .1 );
		float nmr = sqrt_dp_ra * compute_correlation( dp_ra, nm_al, nm_s1, nm_s2, al_s1, al_s1, .1 );
		
		// Use Stouffer's method to combine the z-scores, but weighted by log of sample size
		v_bqr_nums[k] += (bqr * w_dp_ra); v_bqr_dens[k] += (w_dp_ra * w_dp_ra);
		v_mqr_nums[k] += (mqr * w_dp_ra); v_mqr_dens[k] += (w_dp_ra * w_dp_ra);
		v_cyr_nums[k] += (cyr * w_dp_ra); v_cyr_dens[k] += (w_dp_ra * w_dp_ra);
		v_str_nums[k] += (str * w_dp_ra); v_str_dens[k] += (w_dp_ra * w_dp_ra);
		v_nmr_nums[k] += (nmr * w_dp_ra); v_nmr_dens[k] += (w_dp_ra * w_dp_ra);	  
	      }
	      
	      v_ior_nums[k] += (ior * w_dp_q20); v_ior_dens[k] += w_dp_q20;
	      v_nm1_nums[k] += (nm1 * w_al_s1);  v_nm1_dens[k] += w_al_s1;
	      v_nm0_nums[k] += (nm0 * w_ref_s1); v_nm0_dens[k] += w_ref_s1;		
	    }
	    
	    ++k;
	  }
	  odr->close();
	  delete odr;
	  //bcf_clear(v);
	  bcf_destroy(v);
	}
        bcf_hdr_add_sample(odw->hdr, NULL);

        odw->write_hdr();

	std::vector<int32_t> vSex; 
	if ( x_found ) {
	  if ( mSex.empty() ) {
	    fprintf(stderr,"WARNING: No --sex-map is defined with non-PAR X chromosome markers are observed. Assuming everyone is female");
	    vSex.resize(nfiles,2);
	  }
	  else {
	    for(int i=0; i < nfiles; ++i) {
	      const char* sid = bcf_hdr_int2id(odw->hdr, BCF_DT_SAMPLE, i);
	      std::map<std::string,int>::iterator it = mSex.find(sid);
	      if ( it == mSex.end() ) { // not found
		fprintf(stderr,"WARNING: No sex information is available for %s, treating as female\n",sid);
		vSex.push_back(2);
	      }
	      else {
		vSex.push_back(it->second);
	      }
	    }
	  }
	}

	for(size_t k=0; k < v_rids.size(); ++k) {
	  bcf1_t* nv = bcf_init();
	  //bcf_clear(nv);
	  nv->rid = v_rids[k];
	  nv->pos = v_poss[k];
	  nv->rlen = v_rlens[k];
	  nv->n_sample = nfiles;

	  const char* tmp_d_alleles[vn_alleles[k]];
	  for(int l=0; l < vn_alleles[k]; ++l)
	    tmp_d_alleles[l] = v_d_alleles[k][l].c_str();
	  bcf_update_alleles(odw->hdr, nv, tmp_d_alleles, vn_alleles[k]);

	  if ( !v_filts[k].empty() ) {
	    int tmp_filts[v_filts[k].size()];
	    for(size_t l=0; l < v_filts[k].size(); ++l) {
	      tmp_filts[l] = v_filts[k][l];
	    }
	    bcf_update_filter(odw->hdr, nv, tmp_filts, (int32_t)v_filts[k].size());	    
	  }
	  
	  bcf_unpack(nv, BCF_UN_ALL);

	  bool isX = ( ( nv->rid == x_rid ) && ( nv->pos >= xStart ) && ( nv->pos <= xStop ) );
	  
	  // calculate the allele frequencies under HWE. When calculating allele frequencies, the sex information will be ignored
	  float MLE_HWE_AF[vn_alleles[k]];
	  float MLE_HWE_GF[vn_genos[k]];
	  int32_t ploidy = 2; // temporarily constant
	  int32_t n = 0;
	  est->compute_gl_af_hwe(v_pls[k], nfiles, ploidy, vn_alleles[k], MLE_HWE_AF, MLE_HWE_GF,  n, 1e-20);

	  // calculate the genotypes (diploid only)
	  double gp, gp_sum, max_gp;
	  int32_t best_gt;
	  int32_t best_a1, best_a2;
	  int32_t* pls_i;
	  int32_t an = 0;
	  int32_t acs[vn_alleles[k]];
	  int32_t gcs[vn_genos[k]];
	  float afs[vn_alleles[k]];

	  memset(acs, 0, sizeof(int32_t)*vn_alleles[k]);
	  memset(gcs, 0, sizeof(int32_t)*vn_genos[k]);	  
	  
	  for(size_t i=0; i < nfiles; ++i) {
	    pls_i = &v_pls[k][ i * vn_genos[k] ];

	    if ( isX && (vSex[i] == 1) ) { // haploid
	      max_gp = gp_sum = gp = ( est->lt->pl2prob(pls_i[0]) * MLE_HWE_AF[0] );
	      best_gt = 0; best_a1 = 0; best_a2 = 0;
	      for(size_t l=1; l < vn_alleles[k]; ++l) {
		gp = ( est->lt->pl2prob(pls_i[ l*(l+1)/2 + l]) * MLE_HWE_AF[l] );
		gp_sum += gp;
		if ( max_gp < gp ) {
		  max_gp = gp;
		  best_gt = l*(l+1)/2 + l; best_a1 = l; best_a2 = l;
		}		
	      }
	    }
	    else { // diploid
	      max_gp = gp_sum = gp = ( est->lt->pl2prob(pls_i[0]) * MLE_HWE_AF[0] * MLE_HWE_AF[0] );
	      best_gt = 0; best_a1 = 0; best_a2 = 0;
	      for(size_t l=1; l < vn_alleles[k]; ++l) {
		for(size_t m=0; m <= l; ++m) {
		  gp = ( est->lt->pl2prob(pls_i[ l*(l+1)/2 + m]) * MLE_HWE_AF[l] * MLE_HWE_AF[m] * (l == m ? 1 : 2) );
		  gp_sum += gp;
		  if ( max_gp < gp ) {
		    max_gp = gp;
		    best_gt = l*(l+1)/2 + m; best_a1 = m; best_a2 = l;
		  }
		}
	      }
	    }
	      
	    double prob = 1.-max_gp/gp_sum;
	    if ( prob <= 3.162278e-26 )
	      prob = 3.162278e-26;
	    if ( prob > 1 )
	      prob = 1;
	    
	    v_gqs[k][i] = (int32_t)est->lt->prob2pl(prob);
	    
	    if ( ( best_gt > 0 ) && ( v_max_gqs[k] < v_gqs[k][i] ) )
	      v_max_gqs[k] = v_gqs[k][i];
	    
	    v_gts[k][2*i]   = ((best_a1 + 1) << 1);
	    v_gts[k][2*i+1] = ((best_a2 + 1) << 1);	    
	    an += 2;             // still use diploid representation of chrX for now.
	    ++acs[best_a1];
	    ++acs[best_a2];
	    ++gcs[best_gt];
	  }

	  for(size_t i=0; i < vn_alleles[k]; ++i) {
	    afs[i] = acs[i]/(float)an;
	  }

	  bcf_update_format_int32(odw->hdr, nv, "GT", v_gts[k], nfiles * 2);
	  bcf_update_format_int32(odw->hdr, nv, "GQ", v_gqs[k], nfiles );	  
	  bcf_update_format_int32(odw->hdr, nv, "AD", v_ads[k], nfiles * vn_alleles[k]);
	  bcf_update_format_int32(odw->hdr, nv, "OD", v_ods[k], nfiles );	  
	  bcf_update_format_int32(odw->hdr, nv, "PL", v_pls[k], nfiles * vn_genos[k]);

	  float avgdp = (float)v_dp_sums[k]/(float)nfiles;

	  nv->qual = (float) v_max_gqs[k];
	  bcf_update_info_float(odw->hdr, nv, "AVGDP", &avgdp, 1);	  
	  bcf_update_info_int32(odw->hdr, nv, "AC", &acs[1], vn_alleles[k]-1);
	  bcf_update_info_int32(odw->hdr, nv, "AN", &an, 1);
	  bcf_update_info_float(odw->hdr, nv, "AF", &afs[1], vn_alleles[k]-1);
	  bcf_update_info_int32(odw->hdr, nv, "GC", gcs, vn_genos[k]);
	  bcf_update_info_int32(odw->hdr, nv, "GN", &nfiles, 1);

	  if (n) {
	    float* MLE_HWE_AF_PTR = &MLE_HWE_AF[1];
	    bcf_update_info_float(odw->hdr, nv, "HWEAF", MLE_HWE_AF_PTR, vn_alleles[k]-1);
	    //bcf_update_info_float(odw->hdr, nv, "HWEGF", &MLE_HWE_GF, n_genos);
	  }


	  // calculate the allele frequencies under HWD	  
	  float MLE_AF[vn_alleles[k]];
	  float MLE_GF[vn_genos[k]];
	  n = 0;
	  est->compute_gl_af(v_pls[k], nfiles, ploidy, vn_alleles[k], MLE_AF, MLE_GF,  n, 1e-20);
	  if (n) {
	    float* MLE_AF_PTR = &MLE_AF[1];
	    //bcf_update_info_float(odw->hdr, nv, "HWDAF", MLE_AF_PTR, n_alleles-1);
	    bcf_update_info_float(odw->hdr, nv, "HWDGF", &MLE_GF, vn_genos[k]);
	  }

	  if ( isX && !mSex.empty() ) { // copy only female GLs to calculate IBC and HWE_SLP
	    int32_t* p_XX_pls = (int32_t*) malloc(nfiles * vn_genos[k] * sizeof(int32_t));
	    int i, j, l;
	    for(i=0, j=0; i < nfiles; ++i) {
	      if ( vSex[i] == 2 ) {
		for(l=0; l < vn_genos[k]; ++l)  {
		  p_XX_pls[vn_genos[k] * j + l] = v_pls[k][vn_genos[k]*i + l];
		}
		++j;
	      }
	      //v_pls[k][vn_genos[k] * i + l * (l+1) / 2 + m] = 0;	      
	      //p_XX_pls[j++] = v_pls[k][i];
	      //}
	    }

	    float MLE_HWE_AF_XX[vn_alleles[k]];
	    float MLE_HWE_GF_XX[vn_genos[k]];
	    float MLE_AF_XX[vn_alleles[k]];
	    float MLE_GF_XX[vn_genos[k]];

	    // calculate allele frequencies using females
	    est->compute_gl_af_hwe(p_XX_pls, j, ploidy, vn_alleles[k], MLE_HWE_AF_XX, MLE_HWE_GF_XX,  n, 1e-20);
	    est->compute_gl_af(p_XX_pls, j, ploidy, vn_alleles[k], MLE_AF_XX, MLE_GF_XX,  n, 1e-20);


	    for(i=0; i < vn_alleles[k]; ++i) {
	      if ( MLE_HWE_AF_XX[i] < 1e-6 ) MLE_HWE_AF_XX[i] = 1e-6;
	      if ( MLE_AF_XX[i] < 1e-6 ) MLE_AF_XX[i] = 1e-6;	      
	    }
	    
	    for(i=0; i < vn_genos[k]; ++i) {
	      if ( MLE_HWE_GF_XX[i] < 1e-10 ) MLE_HWE_GF_XX[i] = 1e-10;
	      if ( MLE_GF_XX[i] < 1e-10 ) MLE_GF_XX[i] = 1e-10;	      
	    }	    
	      	    
	    float fic = 0;
	    n = 0;
	    est->compute_gl_fic(p_XX_pls, j, ploidy, MLE_HWE_AF_XX, vn_alleles[k], MLE_GF_XX, fic, n);
	    if ( isnanf(fic) ) fic = 0;	  
	    if (n) {
	      bcf_update_info_float(odw->hdr, nv, "IBC", &fic, 1);
	    }
	    
	    // calculate the LRT statistics related to HWE
	    float lrts;
	    float logp;
	    int32_t df;
	    n = 0;
	    est->compute_hwe_lrt(p_XX_pls, j, ploidy, vn_alleles[k], MLE_HWE_GF_XX, MLE_GF_XX, n, lrts, logp, df);
	    if (n) {
	      if ( fic > 0 ) logp = 0-logp;
	      bcf_update_info_float(odw->hdr, nv, "HWE_SLP", &logp, 1);
	    }
	    
	    free(p_XX_pls);
	  }
	  else {
	    float fic = 0;
	    n = 0;
	    est->compute_gl_fic(v_pls[k], nfiles, ploidy, MLE_HWE_AF, vn_alleles[k], MLE_GF, fic, n);
	    if ( isnanf(fic) ) fic = 0;	  
	    if (n) {
	      bcf_update_info_float(odw->hdr, nv, "IBC", &fic, 1);
	    }
	    
	    // calculate the LRT statistics related to HWE
	    float lrts;
	    float logp;
	    int32_t df;
	    n = 0;
	    est->compute_hwe_lrt(v_pls[k], nfiles, ploidy, vn_alleles[k], MLE_HWE_GF, MLE_GF, n, lrts, logp, df);
	    if (n) {
	      if ( fic > 0 ) logp = 0-logp;
	      bcf_update_info_float(odw->hdr, nv, "HWE_SLP", &logp, 1);
	    }
	  }

	  // add additional annotations
	  v_ns_nrefs[k] -= (nfiles - gcs[0]);
	  bcf_update_info_int32(odw->hdr, nv, "NS_NREF", &v_ns_nrefs[k], 1);
	  v_ab_nums[k] /= (v_ab_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "ABE",  &v_ab_nums[k], 1);
	  v_abz_nums[k] /= sqrt(v_abz_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "ABZ",  &v_abz_nums[k], 1);	  	  
	  v_bqr_nums[k] /= sqrt(v_bqr_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "BQZ", &v_bqr_nums[k], 1);	  
	  v_mqr_nums[k] /= sqrt(v_mqr_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "MQZ", &v_mqr_nums[k], 1);	  
	  v_cyr_nums[k] /= sqrt(v_cyr_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "CYZ", &v_cyr_nums[k], 1);	  
	  v_str_nums[k] /= sqrt(v_str_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "STZ", &v_str_nums[k], 1);
	  v_nmr_nums[k] /= sqrt(v_nmr_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "NMZ", &v_nmr_nums[k], 1);
	  v_ior_nums[k] = log(v_ior_nums[k]/v_ior_dens[k]+1e-6)/log(10.); bcf_update_info_float(odw->hdr, nv, "IOR", &v_ior_nums[k], 1);
	  v_nm1_nums[k] /= (v_nm1_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "NM1", &v_nm1_nums[k], 1);
	  v_nm0_nums[k] /= (v_nm0_dens[k]+1e-6); bcf_update_info_float(odw->hdr, nv, "NM0", &v_nm0_nums[k], 1);	  

	  //fprintf(stderr,"AC = %f, AN = %f, NS_NREF = %f, AB = %f, BQR = %f, MQR = %f, CYR = %f, STR = %f, NMR = %f, IOR = %f, NMA = %f\n", acs[1], an, ns_nref, ab_num, bqr_num, mqr_num, cyr_num, str_num, nmr_num, ior_num, nma_num);
	  
	  odw->write(nv);
	  bcf_destroy(nv);	  
	}

	//bcf_destroy(nv);

        odw->close();
	delete odw;
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "paste_and_comput_features v" << version << "\n\n";
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

bool paste_and_compute_features_sequential(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.paste_and_compute_features_sequential();
    igor.print_stats();
    return igor.print;
}

