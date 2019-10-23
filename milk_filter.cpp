/* The MIT License

   Copyright (c) 2014 Adrian Tan and Hyun Min Kang <hmkang@umich.edu>

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

#include "milk_filter.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::vector<std::string> exclude_fam_ids;
    std::string input_vcf_file;
    std::string input_ped_file;
    std::string output_vcf_file;
    bool print;
    int depth_thres;
    std::vector<GenomeInterval> intervals;
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
    NuclearPedigree  *ped;

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

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc = "Likelihood-based filtering from Pedigree.\n"
          "              This is used after the per sample genotyping step in vt.\n"
          "              Input requirements and assumptions:\n"
          "              1. PED file contains only nuclear families where\n"
          "                 a. Each nuclear family is separated by family ID\n"
          "                 b. Each sample ID is globally unique\n"
          "                 c. Duplicated samples are separated by commas\n"
          "              2. VCF files with matching ID to the PED file with GT and PL/GL fields available\n"
          "              Outputs:\n"
          "              1. INFO fields will contain the following fields\n"
          "                a. LUE - -log10 likelihood under assumption of unrelated & HWE.\n"
          "                b. LUD - -log10 likelihood under assumption of unrelated & HWD.\n"
          "                c. LRE - -log10 likelihood under assumption of related & HWE.\n"
          "                d. LRD - -log10 likelihood under assumption of related & HWD.\n"
          "                e. MILK_SCORE - Milk score as LRE-LUD\n";

            version = "0.1";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::SwitchArg arg_print("p", "p", "print options and summary []", cmd, false);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file("b", "b", "input VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_ped_file("f", "f", "input PED file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_intervals("i","i","Intervals[]", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "file", cmd);
            TCLAP::MultiArg<std::string> arg_exclude_fam_ids("x", "x", "Exclude specific family IDs []", false, "str", cmd);
            TCLAP::ValueArg<int> arg_depth_thres("d", "d", "Minimum depth for calculating high-quality Mendelian/Dup concordance (calculated based on AD or DP fields)", false, 15, "int", cmd);
            TCLAP::ValueArg<std::string> arg_sex_map_file("g", "sex-map", "file containing sex map of each individual. ID fitst, X ploidy second", false, "", "file", cmd);
            TCLAP::ValueArg<std::string> arg_xLabel("", "xLabel", "Contig name for X chromosome", false, "X", "str", cmd);
            TCLAP::ValueArg<std::string> arg_yLabel("", "yLabel", "Contig name for Y chromosome", false, "Y", "str", cmd);
            TCLAP::ValueArg<std::string> arg_mtLabel("", "mtLabel", "Contig name for mitochondrial chromosome", false, "MT", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_xStart("", "xStart", "Start base position of non-PAR region in X chromosome", false, 2699520, "int", cmd);
            TCLAP::ValueArg<int32_t> arg_xStop("", "xStop", "End base position of non-PAR region in X chromosome", false, 154931044, "int", cmd);

            cmd.parse(argc, argv);

            const std::vector<std::string>& exclude_fam_ids = arg_exclude_fam_ids.getValue();
            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            input_ped_file = arg_input_ped_file.getValue();
            print = arg_print.getValue();
            depth_thres = arg_depth_thres.getValue();

            sex_map_file = arg_sex_map_file.getValue();
            xLabel = arg_xLabel.getValue();
            yLabel = arg_yLabel.getValue();
            mtLabel = arg_mtLabel.getValue();
            xStart = arg_xStart.getValue();
            xStop = arg_xStop.getValue();
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
        odr = new BCFOrderedReader(input_vcf_file, intervals);
    if ( intervals.size() > 1 ) {
      fprintf(stderr, "[E:%s:%d %s] Multiple intervals are not allowed\n", __FILE__, __LINE__, __FUNCTION__);
      exit(1);
    }

        odw = new BCFOrderedWriter(output_vcf_file, 0);

    // copy header from the original BCF removing all records and sample info (i.e. keep ref info)
    bcf_hdr_t* null_hdr = bcf_hdr_subset(odr->hdr, 0, NULL, NULL);
    //bcf_hdr_remove(null_hdr, BCF_HL_INFO, NULL);
    //bcf_hdr_remove(null_hdr, BCF_HL_FMT, NULL);
    //bcf_hdr_remove(null_hdr, BCF_HL_FLT, NULL);
        odw->set_hdr(null_hdr);
    bcf_hdr_destroy(null_hdr);
    bcf_hdr_remove(odw->hdr, BCF_HL_FMT, NULL);

        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=LRE,Number=1,Type=Float,Description=\"Per-family log-likelihood assuming HWE and Mendelian concordance\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=LUD,Number=1,Type=Float,Description=\"Per-family log-likelihood ignoring HWE and Mendelian concordance\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=BF,Number=1,Type=Float,Description=\"Per-family Milk Bayes Factor score\">\n");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=FAMGT,Number=1,Type=String,Description=\"Genotypes for each family. 0 for missing, 1, 2, 3 for possible genotypes. Individuals are separated by underscore. Duplicated samples are concatenated without comma\">\n");

        //bcf_hdr_append(odw->hdr, "##INFO=<ID=LUE,Number=1,Type=Float,Description=\"log10 likelihood assuming unrelatedness and Hardy-Weinberg Equilibrium\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=MILK_LUD,Number=1,Type=Float,Description=\"log10 likelihood assuming unrelatedness and Hardy-Weinberg Disequilibrium\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=MILK_LRE,Number=1,Type=Float,Description=\"log10 likelihood assuming relatedness and Hardy-Weinberg Equilibrium\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=MILK_HWEAF,Number=A,Type=Float,Description=\"Maximum-likelihood allele frequency estimates using Hardy-Weinberg Equilibrium\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=MILK_HWDGF,Number=G,Type=Float,Description=\"Maximum-likelihood genotype frequency estimates ignoring Hardy-Weinberg Equilibrium\">\n");
        //bcf_hdr_append(odw->hdr, "##INFO=<ID=LRD,Number=1,Type=Float,Description=\"log10 likelihood assuming relatedness and Hardy-Weinberg Disequilibrium\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=MILK_BF,Number=1,Type=Float,Description=\"Bayes factor as LRE-LUD\">\n");
    if ( !sex_map_file.empty() ) {
      bcf_hdr_append(odw->hdr, "##INFO=<ID=MILK_CHRX,Number=1,Type=Float,Description=\"Bayes factor based on chrX heterozygosity\">\n");
    }
        bcf_hdr_append(odw->hdr, "##INFO=<ID=DUP_CONC_ALL,Number=9,Type=Integer,Description=\"Duplicated sample joint genotype across every dup pairs\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=DUP_CONC_THRES,Number=9,Type=Integer,Description=\"Duplicated sample joint genotypes across every dup pair passing the minimum depth threshold \">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=TRIO_CONC_ALL,Number=27,Type=Integer,Description=\"Trio sample joint genotypes across every offspring\">\n");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=TRIO_CONC_THRES,Number=27,Type=Integer,Description=\"Trio sample joint genotypes across every trio passing the minimum depth threshold\">\n");

    fprintf(stderr,"Loading pedigree file %s\n", input_ped_file.c_str());
    ped = new NuclearPedigree(input_ped_file.c_str());

    // match between pedigree's indivduals and BCF's individuals
    // what we need is the indices of individual IDs in each family
    int n = bcf_hdr_nsamples(odr->hdr);
    for(int i=0; i < n; ++i) {
      ped->setSampleIndex(odr->hdr->samples[i], i);
    }

    fprintf(stderr,"In the Original Pedigree: %d families, %d individuals, %d sequenced samples\n", (int)ped->famIDmap.size(), ped->numPeople(), (int)ped->smIDmap.size());

    int nRemoved = ped->removeSamplesWithoutIndex();
    fprintf(stderr,"Removed %d samples not in the VCF file from the pedigree\n", nRemoved);

    fprintf(stderr,"Overlapping Samples Only : %d families, %d individuals, %d sequenced samples\n", (int)ped->famIDmap.size(), ped->numPeople(), ped->numSamplesWithIndex());

    // write the pedigree information for every individual
    std::map<std::string, NuclearFamily*>::iterator it;
    for(it = ped->famIDmap.begin(); it != ped->famIDmap.end(); ++it) {
      bcf_hdr_add_sample(odw->hdr, it->first.c_str()); // family ID is the sample name for output

      NuclearFamily* pFam = it->second;

      if ( ( pFam->pDad ) && ( pFam->pDad->samples.size() > 1 ) ) {
        for(int i=1; i < pFam->pDad->samples.size(); ++i)
          bcf_hdr_append(odw->hdr, (std::string("##PEDIGREE=<ID=") + pFam->pDad->samples[i]->sampleID + ",Original=" + pFam->pDad->samples[0]->sampleID + ",FamilyID="+pFam->famID+">\n").c_str());
        bcf_hdr_append(odw->hdr, (std::string("##PEDIGREE=<ID=") + pFam->pDad->samples[0]->sampleID + ",Father=.,Mother=.,FamilyID="+pFam->famID+">\n").c_str());
      }

      if ( ( pFam->pMom ) && ( pFam->pMom->samples.size() > 1 ) ) {
        for(int i=1; i < pFam->pMom->samples.size(); ++i)
          bcf_hdr_append(odw->hdr, (std::string("##PEDIGREE=<ID=") + pFam->pMom->samples[i]->sampleID + ",Original=" + pFam->pMom->samples[0]->sampleID + ",FamilyID="+pFam->famID+">\n").c_str());
        bcf_hdr_append(odw->hdr, (std::string("##PEDIGREE=<ID=") + pFam->pMom->samples[0]->sampleID + ",Father=.,Mother=.,FamilyID="+pFam->famID+">\n").c_str());
      }

      for(int j=0; j < pFam->pKids.size(); ++j) {
        if ( pFam->pKids[j]->samples.size() > 1 ) {
          for(int i=1; i < pFam->pKids[j]->samples.size(); ++i)
        bcf_hdr_append(odw->hdr, (std::string("##PEDIGREE=<ID=") + pFam->pKids[j]->samples[i]->sampleID + ",Original=" + pFam->pKids[j]->samples[0]->sampleID + ",FamilyID="+pFam->famID+">\n").c_str());
          bcf_hdr_append(odw->hdr, (std::string("##PEDIGREE=<ID=") + pFam->pKids[j]->samples[0]->sampleID + ",Father=" + (pFam->pDad ? pFam->pDad->samples[0]->sampleID : ".") + ",Mother=" + (pFam->pMom ? pFam->pMom->samples[0]->sampleID : ".") + ",FamilyID=" + pFam->famID + ">\n").c_str());
        }
      }
    }
    bcf_hdr_add_sample(odw->hdr, NULL);

    odw->write_hdr();


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

    float logSumPosterior( int32_t* pls, NuclearFamilyPerson* pPerson, float* priors, bool mergeDups) {
      if ( pPerson == NULL ) {
    return 0;
      }
      else if ( mergeDups ) {
    int32_t sumPLs[3] = {0, 0, 0};
    for(int i=0; i < (int)pPerson->samples.size(); ++i) {
      int idx = pPerson->samples[i]->index;
      if ( idx >= 0 ) {
        sumPLs[0] += pls[idx * 3];
        sumPLs[1] += pls[idx * 3 + 1];
        sumPLs[2] += pls[idx * 3 + 2];
      }
    }
    return ( (float)log(priors[0] * LogTool::pl2prob(sumPLs[0]) + priors[1] * LogTool::pl2prob(sumPLs[1]) + priors[2] * LogTool::pl2prob(sumPLs[2])) );
      }
      else { // don't merge dups and treat them as independent samples
    float logSum = 0;
    for(int i=0; i < (int)pPerson->samples.size(); ++i) {
      int idx = pPerson->samples[i]->index;
      if ( idx >= 0 ) {
        logSum += log(priors[0] * LogTool::pl2prob(pls[idx * 3]) + priors[1] * LogTool::pl2prob(pls[idx * 3 + 1]) + priors[2] * LogTool::pl2prob(pls[idx * 3 + 2]) );
      }
    }
    return ( logSum );
      }
    }

    void addDupJoint( const std::vector<int>& indGTs, const std::vector<int>& indDPs, int* jointAll, int* jointThres ) {
      for(int i=1; i < (int)indGTs.size(); ++i) {
    for(int j=0; j < i; ++j) {
      if ( ( indGTs[i] > 0 ) && ( indGTs[j] > 0 ) ) {
        int cell = (indGTs[i]-1) * 3 + (indGTs[j]-1);
        ++jointAll[cell];
        if ( ( indDPs[i] >= depth_thres ) || ( indDPs[j] >= depth_thres ) ) {
          ++jointThres[cell];
        }
      }
    }
      }
    }

    void addTrioJoint( int kidGT, int kidDP, int dadGT, int dadDP, int momGT, int momDP, int* jointAll, int* jointThres ) {
      if ( ( dadGT > 0 ) && ( momGT > 0 ) && ( kidGT > 0 ) ) {
    int cell = (dadGT-1) * 9 + (momGT-1) * 3 + (kidGT-1);
    ++jointAll[cell];
    if ( ( dadDP >= depth_thres ) && ( momDP >= depth_thres ) && ( kidDP >= depth_thres ) )
      ++jointThres[cell];
      }
    }

    void getPersonGTstr( const std::vector<int>& vGTs, std::string& s ) {
      for(int i=0; i < (int)vGTs.size(); ++i)
    s += ('0' + vGTs[i]);
    }

    float logGL( int32_t* pls, NuclearFamilyPerson* pPerson, int genotype) {
      if ( pPerson == NULL ) return 0;
      else {
    int32_t sumPL = 0;
    for(int i=0; i < (int)pPerson->samples.size(); ++i) {
      int idx = pPerson->samples[i]->index;
      if ( idx >= 0 ) {
        sumPL += pls[3 * idx + genotype];
      }
    }
    return (-0.1 * sumPL);
      }
    }

    int getPersonGenoDepth( int32_t* gts, int32_t* dps, NuclearFamilyPerson* pPerson, std::vector<int>& genos, std::vector<int>& depths) {
      genos.clear();
      depths.clear();
      if ( pPerson == NULL ) return 0;
      else {
    int nSamples = 0;
    for(int i=0; i < (int)pPerson->samples.size(); ++i) {
      int idx = pPerson->samples[i]->index;
      if ( idx >= 0 ) {
        int g1 = gts[2*idx];
        int g2 = gts[2*idx+1];
        int geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          geno = 0;
        }
        else {
          geno = bcf_alleles2gt(bcf_gt_allele(g1),bcf_gt_allele(g2))+1;
        }
        //fprintf(stderr,"%s %d %d %d %d\n",pPerson->samples[i]->sampleID.c_str(), idx, g1, g2, geno);
        genos.push_back(geno);
        depths.push_back(dps[idx]);
        ++nSamples;
      }
    }
    return nSamples;
      }
    }

    void milk_filter() {
      int n = bcf_hdr_nsamples(odr->hdr);

      // load the BCF file
      bcf1_t* v = bcf_init();

      // store the PL information in a separate array
      int32_t np_gt = 0;
      int32_t np_pl = 0;
      int32_t np_dp = 0;
      int32_t* p_gt = NULL;
      int32_t* p_pl = NULL;
      int32_t* p_dp = NULL;

      int ns = ped->numSamplesWithIndex();
      int32_t* p_pl_ns = (int32_t*)calloc(ns * 3, sizeof(int32_t));
      float logTen = log(10.);

      // determine rids for X chromosome
      int32_t x_rid = bcf_hdr_name2id(odr->hdr,xLabel.c_str());
      std::vector<int32_t> vSex;

      int nskip = 0, nread = 0;
      for(nread = 0; odr->read(v); ++nread) {
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

    if ( skip ) {
      ++nskip;
      continue;
    }

    bool isX = ( v->rid == x_rid ) && ( v->pos >= xStart ) && ( v->pos <= xStop );
    if ( isX && vSex.empty() && !mSex.empty() ) { // populate vSex
      for(int i=0; i < n; ++i) {
        std::map<std::string,int>::iterator it = mSex.find(odr->hdr->samples[i]);
        if ( it == mSex.end() ) { // not found
          fprintf(stderr,"WARNING: No sex information is available for %s, treating as female\n",odr->hdr->samples[i]);
          vSex.push_back(2);
        }
        else {
          vSex.push_back(it->second);
        }
      }
    }

    bcf1_t* nv = bcf_init();

    nv->n_sample = n;

    // copy the site info
    nv->rid = v->rid;
    nv->pos = v->pos;
    nv->rlen = v->rlen;
    nv->qual = v->qual;
    //nv->n_info = v->n_info;
    nv->n_allele = v->n_allele;
    //nv->shared.m = nv->shared.l = v->shared.l;
    //nv->shared.s = (char*) malloc(nv->shared.l);
    //memcpy(nv->shared.s, v->shared.s, nv->shared.l);

    bcf_update_alleles(odw->hdr, nv, (const char**)v->d.allele, v->n_allele);
    bcf_update_filter(odw->hdr, nv, v->d.flt, v->d.n_flt);
    //nv->n_info = v->n_info;
    //nv->qual = v->qual;

    bcf_unpack(nv, BCF_UN_ALL);

    if ( (nread + 1) % 1000 == 0 )
      fprintf(stderr,"Processing %d variants at %s:%d\n", nread + 1, bcf_seqname(odw->hdr,nv), (int32_t) nv->pos);

    // get GT fields
    if ( bcf_get_genotypes(odr->hdr, v, &p_gt, &np_gt) < 0 ) {
      fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- GT\n", __FILE__, __LINE__, __FUNCTION__);
      exit(1);
    }

    // get PL fields
    if ( bcf_get_format_int32(odr->hdr, v, "PL", &p_pl, &np_pl) < 0 ) {
      fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- PL\n", __FILE__, __LINE__, __FUNCTION__);
      exit(1);
    }

    // cap the max PL value to 60
    int i;
    for(i=0; i < np_pl; ++i)
      if ( p_pl[i] > 255 ) p_pl[i] = 255;

    // get DP or AD fields
    // get PL fields
    if ( bcf_get_format_int32(odr->hdr, v, "DP", &p_dp, &np_dp) < 0 ) {
      if ( bcf_get_format_int32(odr->hdr, v, "AD", &p_dp, &np_dp) < 0 ) {
        fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- PL or AD\n", __FILE__, __LINE__, __FUNCTION__);
        exit(1);
      }

      // if AD field is available, use their sum as depth (assuming biallelics);
      for(int i=0; i < n; ++i) {
        p_dp[i] = p_dp[2*i] + p_dp[2*i+1];
      }
    }

    // copy PL fields on the subset of samples with pedigree information
    i = 0;
    std::map<std::string, NuclearFamilySample*>::iterator it;
    for(it = ped->smIDmap.begin(); it != ped->smIDmap.end(); ++it) {
      if ( it->second->index >= 0 ) {
        int j = it->second->index;
        p_pl_ns[i*3] = p_pl[j*3];
        p_pl_ns[i*3+1] = p_pl[j*3+1];
        p_pl_ns[i*3+2] = p_pl[j*3+2];
        ++i;
      }
    }

    if ( isX && !vSex.empty() ) { // if it is X chromosome
      float* hweAFs = NULL;
      int n_hweAF = 0;
      if ( bcf_get_info_float(odr->hdr, v, "HWEAF", &hweAFs, &n_hweAF) < 0 ) {
        fprintf(stderr, "[E:%s:%d %s] FORMAT field does not contain expected fields -- HWEAF\n", __FILE__, __LINE__, __FUNCTION__);
        exit(1);
      }
      float hweAF = hweAFs[0];
      free(hweAFs);
      if ( hweAF < 1e-6 ) { hweAF = 1e-6; }
      else if ( 1-hweAF < 1e-6 ) { hweAF = 1-1e-6; }

      float llkMaleDipl = 0;
      float llkMaleHapl = 0;
      // calculate the likelihoods accounting for sex vs not
      for(int i=0; i < n; ++i) {
        if ( vSex[i] == 2 ) { // diplod female
          // this does not change the overall likelihood, so just skip it
        }
        else if ( vSex[i] == 1 ) { // haploid male
          llkMaleDipl += (float)log( (1.-hweAF) * (1.-hweAF) * LogTool::pl2prob(p_pl[3*i]) + 2 * hweAF * (1.-hweAF) * LogTool::pl2prob(p_pl[3*i+1]) + hweAF * hweAF * LogTool::pl2prob(p_pl[3*i+2]) );
          llkMaleHapl += (float)log( (1.-hweAF) * LogTool::pl2prob(p_pl[3*i]) + hweAF * LogTool::pl2prob(p_pl[3*i+2]) );
        }
        else {
          fprintf(stderr, "[E:%s:%d %s] Unexpected sex %d is observed in individual index %d\n", __FILE__, __LINE__, __FUNCTION__, vSex[i], i);
          exit(1);
        }
      }

      float milkChrX = (llkMaleHapl-llkMaleDipl)/(float)log(10);
      if ( bcf_update_info_float(odw->hdr, nv, "MILK_CHRX", &milkChrX, 1) ) {
        fprintf(stderr, "[E:%s:%d %s] Error in updating MILK_CHRX FORMAT field\n", __FILE__, __LINE__, __FUNCTION__);
        abort();
      }
    }

    // estimate the allele frequencies (ignoring the ploidy info)
    float MLE_HWE_AF[2];
    float MLE_HWE_GF[3];
    float MLE_HWD_AF[2];
    float MLE_HWD_GF[3];
    float LOG_MLE_HWE_AF[2];
    float LOG_MLE_HWE_GF[3];
    float LOG_MLE_HWD_GF[3];
    int32_t ploidy = 2; // temporarily constant
    int32_t n1 = 0;
    int32_t n2 = 0;
    Estimator::compute_gl_af_hwe(p_pl_ns, ns, ploidy, 2, MLE_HWE_AF, MLE_HWE_GF, n1, 1e-20);
    Estimator::compute_gl_af(p_pl_ns, ns, ploidy, 2, MLE_HWD_AF, MLE_HWD_GF, n2, 1e-20);

    // minimum genotype frequency is set to 1e-10 to reduce outlier effects and underflow
    for(int i=0; i < 2; ++i) {
      if ( MLE_HWE_AF[i] < 1e-6 )
        MLE_HWE_AF[i] = 1e-6;
      LOG_MLE_HWE_AF[i] = logf(MLE_HWE_AF[i]);
    }

    for(int i=0; i < 3; ++i) {
      if ( MLE_HWE_GF[i] < 1e-10 ) MLE_HWE_GF[i] = 1e-10;
      if ( MLE_HWD_GF[i] < 1e-10 ) MLE_HWD_GF[i] = 1e-10;
      LOG_MLE_HWE_GF[i] = logf(MLE_HWE_GF[i]);
      LOG_MLE_HWD_GF[i] = logf(MLE_HWD_GF[i]);
    }


    // calculate the overall likelihood under HWE and HWD, ignoring the pedigree info
    float lue = 0;
    float lud = 0;

    for(i=0; i < ns; ++i) {
      lue += log( LogTool::pl2prob(p_pl_ns[3*i]) * MLE_HWE_GF[0] + LogTool::pl2prob(p_pl_ns[3*i+1]) * MLE_HWE_GF[1]  + LogTool::pl2prob(p_pl_ns[3*i+2]) * MLE_HWE_GF[2] );
      lud += log( LogTool::pl2prob(p_pl_ns[3*i]) * MLE_HWD_GF[0] + LogTool::pl2prob(p_pl_ns[3*i+1]) * MLE_HWD_GF[1]  + LogTool::pl2prob(p_pl_ns[3*i+2]) * MLE_HWD_GF[2] );
    }

    // expected probability of offspring genotypes given parents' genotypes
    float pomf[9][3] = { { 1,  0, 0}, { .5, .5,   0},   {0,  1,  0},
                 {.5, .5, 0}, {.25, .5, .25},   {0, .5, .5},
                 { 0,  1, 0}, {  0, .5,  .5},   {0,  0,  1} };

    // expected probability of offspring genotypes given parents' genotypes
    // when kid is male,
    float pomf_maleX[9][3] = { {1, 0, 0}, {.5, 0, .5}, {0, 0, 1},
                   {.5, 0, .5}, {.5, 0, .5}, {.5, 0, .5}, // dummy values - shouldn't matter
                   {1, 0, 0}, {.5, 0, .5}, {0, 0, 1} };

    int nfams = ped->famIDmap.size();
    //std::vector<double> fLUEs(0, nfams);
    float* fLUDs = (float*) calloc(nfams, sizeof(float));
    float* fLREs = (float*) calloc(nfams, sizeof(float));
    float* fBFs = (float*) calloc(nfams, sizeof(float));
      //std::vector<float> fLUDs(nfams, 0);
      //std::vector<float> fLREs(nfams, 0);
      //std::vector<float> fBFs(nfams, 0);
    std::vector<std::string> fGTs(nfams);
    float sumLUD = 0;
    float sumLRE = 0;

    int32_t dupJointAll[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int32_t dupJointThres[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int32_t trioJointAll[27]   = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int32_t trioJointThres[27] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector<std::string> famGTs(nfams);
    // calculate the genotype concordance for related and duplicated samples

    std::vector<int32_t> dadGTs;
    std::vector<int32_t> momGTs;
    std::vector<int32_t> nKids;
    std::vector< std::vector<int32_t> > kidGTs;
    std::vector<int32_t> dadDPs;
    std::vector<int32_t> momDPs;
    std::vector< std::vector<int32_t> > kidDPs;

    // calculate the likelihoods using the pedigree info
    // also calculate duplicate concordance
    std::map<std::string, NuclearFamily*>::iterator it2;
    i = 0;
    for(it2 = ped->famIDmap.begin(); it2 != ped->famIDmap.end(); ++it2, ++i) {
      NuclearFamily* pFam = it2->second;

      // calculate per-family likelihood ignoring relatedness/dups
      fLUDs[i] += logSumPosterior( p_pl, pFam->pDad, MLE_HWD_GF, false )/logTen;
      fLUDs[i] += logSumPosterior( p_pl, pFam->pMom, MLE_HWD_GF, false )/logTen;
      for(int j=0; j < (size_t)pFam->pKids.size(); ++j) {
        fLUDs[i] += logSumPosterior( p_pl, pFam->pKids[j], MLE_HWD_GF, false )/logTen;
      }


      // calculate per-family likelihood accounting for relatedness and ploidy
      // \prod_o \Pr(o|m,f)
      // calculate the sum of the likelihood
      for(int j=0; j < 3; ++j) {   // paternal genotypes
        float logGLDad = logGL(p_pl, pFam->pDad, j);

        for(int k=0; k < 3; ++k) { // maternal genotypes
          float logGLMom = logGL(p_pl, pFam->pMom, k);

          float logSumKidPosterior = ( logGLDad + logGLMom + LOG_MLE_HWE_GF[j] + LOG_MLE_HWE_GF[k] );
          if ( isX && !vSex.empty() ) {
        if ( j == 1 ) { // if dad is heterozygous, skip it
          continue;
          //logSumKidPosterior = ( logGLDad + logGLMom - 100 + LOG_MLE_HWE_GF[k] );
        }
        else {
          logSumKidPosterior = ( logGLDad + logGLMom + LOG_MLE_HWE_AF[j/2] + LOG_MLE_HWE_GF[k] );
        }
          }

          for(int l=0; l < (size_t)pFam->pKids.size(); ++l) {
        if ( isX && !vSex.empty() && pFam->pKids[l]->sex == 1 ) {
          logSumKidPosterior += logSumPosterior( p_pl, pFam->pKids[l], pomf_maleX[j*3+k], true );
        }
        else {
          logSumKidPosterior += logSumPosterior( p_pl, pFam->pKids[l], pomf[j*3+k], true );
        }
          }

          if ( j + k == 0 ) {
        //fLREs[i] = fLRDs[i] = llogSumKidPosterior;
        fLREs[i] = logSumKidPosterior;
          }
          else if ( fLREs[i] > logSumKidPosterior ) {
        fLREs[i] += (float)log(1. + exp((double)(logSumKidPosterior - fLREs[i])) );
        //fLRDs[i] += log ( 1 + exp(logSumKidPosterior - fLRDs[i]) );
          }
          else {
        fLREs[i] = logSumKidPosterior + (float)log ( 1. + (double)exp(fLREs[i] - logSumKidPosterior) );
        //fLRDs[i] = logSumKidPosterior + log ( 1 + exp(fLRDs[i] - logSumKidPosterior) );
          }
          //fLREs[i] += ( logGLDad + logGLMom + LOG_MLE_HWE_GF[j] + LOG_MLE_HWE_GF[k] );
          //fLRDs[i] += ( logGLDad + logGLMom + LOG_MLE_HWD_GF[j] + LOG_MLD_HWE_GF[k] );
        }
      }
      fLREs[i] /= logTen;
      fBFs[i] = fLREs[i] - fLUDs[i];
      //fLRDs[i] /= logTen;

      sumLUD += fLUDs[i];
      sumLRE += fLREs[i];


      // calculate genotype concordance
      // first get genotypes
      int nDad = getPersonGenoDepth( p_gt, p_dp, pFam->pDad, dadGTs, dadDPs);
      int nMom = getPersonGenoDepth( p_gt, p_dp, pFam->pMom, momGTs, momDPs);
      nKids.resize(pFam->pKids.size());
      kidGTs.resize(pFam->pKids.size());
      kidDPs.resize(pFam->pKids.size());
      for(int j=0; j < (int)pFam->pKids.size(); ++j) {
        nKids[j] = getPersonGenoDepth( p_gt, p_dp, pFam->pKids[j], kidGTs[j], kidDPs[j]);
      }

      // get the duplicate concordance and trio concordance
      if ( nDad > 1 ) addDupJoint( dadGTs, dadDPs, dupJointAll, dupJointThres );
      if ( nMom > 1 ) addDupJoint( momGTs, momDPs, dupJointAll, dupJointThres );
      for(int j=0; j < (int)pFam->pKids.size(); ++j) {
        if ( nKids[j] > 1 ) addDupJoint( kidGTs[j], kidDPs[j], dupJointAll, dupJointThres );
        if ( ( nKids[j] > 0 ) && ( nDad > 0 ) && ( nMom > 0 ) ) {
          addTrioJoint( kidGTs[j][0], kidDPs[j][0], dadGTs[0], dadDPs[0], momGTs[0], momDPs[0], trioJointAll, trioJointThres );
        }
      }

      getPersonGTstr( dadGTs, fGTs[i] );
      fGTs[i] += '_';
      getPersonGTstr( momGTs, fGTs[i] );
      for(int j=0; j < (int)pFam->pKids.size(); ++j) {
        fGTs[i] += '_';
        getPersonGTstr( kidGTs[j], fGTs[i] );
      }
    }

    float milkBF = sumLRE - sumLUD;

    //fprintf(stderr,"%f %f\n",lud/logTen,sumLUD);

    // write

        const char** s_fGTs = (const char**)malloc(sizeof(char*)*nfams);
    for(i=0; i < nfams; ++i) {
      s_fGTs[i] = fGTs[i].c_str();
    }
    //std::vector<const char*> s_fGTs;
    //for(i=0; i < (int)fGTs.size(); ++i)
    //  s_fGTs.push_back(fGTs[i].c_str());

    //fprintf(stderr,"LRE\n");
    if ( bcf_update_format_float(odw->hdr, nv, "LRE", fLREs, nfams) ) {
      fprintf(stderr, "[E:%s:%d %s] Error in updating LRE FORMAT field\n", __FILE__, __LINE__, __FUNCTION__);
      abort();
    }
    //fprintf(stderr,"LUD\n");
    if ( bcf_update_format_float(odw->hdr, nv, "LUD", fLUDs, nfams) ) {
      fprintf(stderr, "[E:%s:%d %s] Error in updating LUD FORMAT field\n", __FILE__, __LINE__, __FUNCTION__);
      abort();
    }
    //fprintf(stderr,"BF\n");
    if ( bcf_update_format_float(odw->hdr, nv, "BF", fBFs, nfams) ) {
      fprintf(stderr, "[E:%s:%d %s] Error in updating BF FORMAT field\n", __FILE__, __LINE__, __FUNCTION__);
      abort();
    }
    //fprintf(stderr,"FAMGT\n");

    if (! bcf_hdr_idinfo_exists(odw->hdr, BCF_HL_FMT, bcf_hdr_id2int(odw->hdr,BCF_DT_ID,"FAMGT")) ) {
      fprintf(stderr,"Cannot find FAMGT from header\n");
      abort();
    }

    if ( int ret = bcf_update_format_string(odw->hdr, nv, "FAMGT", s_fGTs, nfams) ) {
      fprintf(stderr, "[E:%s:%d %s] Error in updating FAMGT FORMAT field, return vale is %d\n", __FILE__, __LINE__, __FUNCTION__, ret);
      abort();
    }


    free(s_fGTs);
    free(fLUDs);
    free(fLREs);
    free(fBFs);

    // copy original INFO field to output
    //void* infoValues = NULL;
    //int nInfo = 0;
    //for(i=0; i < v->n_info; ++i) {
      //const char* key = bcf_hdr_id2hrec(odr->hdr,BCF_DT_ID,BCF_HL_INFO,v->d.info[i].key)->vals[0];
        //odr->hdr->hrec[v->d.info[i].key]->key;
      //fprintf(stderr,"Adding INFO key %s at %d (%d)\n", key, v->d.info[i].key, v->d.info[i].type);
      //bcf_get_info_values(odr->hdr, v, key, &infoValues, &nInfo, v->d.info[i].type);
      //bcf_info_t* bcf_get_info_id( v, v->d.info[i].key );
      //bcf_update_info(odw->hdr, nv, key, infoValues, nInfo, v->d.info[i].type);
    //}

    // copy the INFO fields from v to nv
    for(i=0; i < v->n_info; ++i) {
    //for(i=0; i < 2; ++i) {
      bcf_info_t& info = v->d.info[i];
      //fprintf(stderr,"info->key = %d, type = %d, len = %d, tag = %s, var_len = %d\n", info.key, info.type, info.len, odr->hdr->id[BCF_DT_ID][info.key].key, bcf_hdr_id2length(odr->hdr,BCF_HL_INFO,info.key));
      if ( info.type != BCF_BT_NULL ) {
        const char* tag = bcf_hdr_int2id(odr->hdr,BCF_DT_ID,info.key);
        int htype = bcf_hdr_id2type(odr->hdr,BCF_HL_INFO,info.key);
        int ntmp_arr = 0;
        void* tmp_arr = NULL;
        int ret = bcf_get_info_values(odr->hdr, v, tag, &tmp_arr, &ntmp_arr, htype);
        if ( ret > 0 ) {
          if ( bcf_update_info(odw->hdr, nv, tag, tmp_arr, ntmp_arr, htype) < 0 ) {
        fprintf(stderr,"Cannot write INFO field %s\n",tag);
        abort();
          }
        }
        else {
          fprintf(stderr,"Cannot retrieve INFO field %s\n",tag);
          abort();
        }
        free(tmp_arr);
      }

      /*
      if (info.key<0) {
        fprintf(stderr, "[E::%s] invalid BCF, the INFO key key=%d not present in the header.\n", __func__, info.key);
        abort();
      }
      const char* tag = bcf_hdr_int2id(odr->hdr,BCF_DT_ID,info.key); //odr->hdr->id[BCF_DT_ID][info.key].key;
      int32_t var_len = bcf_hdr_id2length(odr->hdr,BCF_HL_INFO,info.key);
      int32_t type = info.type;

      if (type==BCF_BT_INT8||type==BCF_BT_INT16||type==BCF_BT_INT32) {
        int32_t n = 0;
        int32_t* g = 0;
        int32_t ret = bcf_get_info_int32(odr->hdr, v, tag, &g, &n);
        if (ret>0) {
          bcf_update_info_int32(odw->hdr, nv, tag, g, n);
          free(g);
        }
      }
      else if (type==BCF_BT_FLOAT) {
        int32_t n = 0;
        float* g = 0;
        int32_t ret = bcf_get_info_float(odr->hdr, v, tag, &g, &n);
        if (ret>0) {
          bcf_update_info_float(odw->hdr, nv, tag, g, n);
          free(g);
        }
      }
      else if (type==BCF_BT_CHAR) {
        int32_t n = 0;
        char** g = 0;
        int32_t ret = bcf_get_info_string(odr->hdr, v, tag, &g, &n);
        if (ret>0) {
          bcf_update_info_string(odw->hdr, nv, tag, g);
          free(g);
        }
      }
      else {
        fprintf(stderr, "[E::%s] invalid BCF, Cannot recognize the INFO field type %d\n", __func__,type);
        abort();
      }
      */
    }
    //exit(1);


    bcf_update_info_float(odw->hdr, nv, "MILK_LRE", &sumLRE, 1);
    bcf_update_info_float(odw->hdr, nv, "MILK_LUD", &sumLUD, 1);
    bcf_update_info_float(odw->hdr, nv, "MILK_BF", &milkBF, 1);
    bcf_update_info_float(odw->hdr, nv, "MILK_HWEAF", &MLE_HWE_AF[1], 1);
    bcf_update_info_float(odw->hdr, nv, "MILK_HWDGF", MLE_HWD_GF, 3);
    bcf_update_info_int32(odw->hdr, nv, "DUP_CONC_ALL", dupJointAll, 9);
    bcf_update_info_int32(odw->hdr, nv, "DUP_CONC_THRES", dupJointThres, 9);
    bcf_update_info_int32(odw->hdr, nv, "TRIO_CONC_ALL", trioJointAll, 27);
    bcf_update_info_int32(odw->hdr, nv, "TRIO_CONC_THRES", trioJointThres, 27);

    //fprintf(stderr,"INFO\n");

    odw->write(nv);
    bcf_destroy(nv);
      }

      fprintf(stderr,"Procesed %d variants and skipped %d non-biallelic variants\n", nread, nskip);

      odw->close();
      delete odw;

      fprintf(stderr,"Finished writing output VCF/BCF file\n");
    };

    void print_options()
    {
        if (!print) return;

        std::clog << "milk_filter v" << version << "\n\n";
        std::clog << "         [i] input VCF file       "  << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file       " << output_vcf_file << "\n";
        std::clog << "\n";
    }

    void print_stats()
    {
        if (!print) return;

        std::clog << "\n";
        //std::cerr << "stats: Total number of files pasted  " << input_vcf_files.size() << "\n";
        std::clog << "\n";
    };

    ~Igor() {};

    private:
};

}

bool milk_filter(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.milk_filter();
    //igor.print_stats();
    return igor.print;
}

