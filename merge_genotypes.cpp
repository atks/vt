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

#include "merge_genotypes.h"

#define SINGLE     0
#define AGGREGATED 1

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
    std::string candidate_sites_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    ///////
    //i/o//
    ///////
    BCFSyncedReader *sr;
    BCFOrderedWriter *odw;
    bcf1_t *v;

    ///////////////
    //general use//
    ///////////////
    std::vector<int32_t> file_types;
    kstring_t variant;

    //////////
    //filter//
    //////////
    std::string fexp;
    Filter filter;
    bool filter_exists;
    
    /////////
    //stats//
    /////////
    uint32_t no_samples;
    uint32_t no_snps;
    uint32_t no_indels;
    uint32_t no_vntrs;

    /////////
    //tools//
    /////////
    VariantManip * vm;

    Igor(int argc, char ** argv)
    {
        //////////////////////////
        //options initialization//
        //////////////////////////
        try
        {
            std::string desc =
"Merge GT, PL, DP, ADF and ADR across samples.\n\
Extracts only the naive genotypes based on best guess genotypes.";

            version = "0.5";
            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_candidate_sites_vcf_file("c", "c", "candidate sites VCF file with annotation []", true, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_fexp("f", "f", "filter expression (applied to candidate sites file)[]", false, "", "str", cmd);
            TCLAP::UnlabeledMultiArg<std::string> arg_input_vcf_files("<in1.vcf>...", "Multiple VCF files",false, "files", cmd);

            cmd.parse(argc, argv);

            parse_files(input_vcf_files, arg_input_vcf_files.getValue(), arg_input_vcf_file_list.getValue());
            candidate_sites_vcf_file = arg_candidate_sites_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            fexp = arg_fexp.getValue();
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
        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter.parse(fexp.c_str(), false);
        filter_exists = fexp=="" ? false : true;

        //////////////////////
        //i/o initialization//
        //////////////////////
        fprintf(stderr, "[I:%s:%d %s] Initializing %zd VCF files ...", __FILE__, __LINE__, __FUNCTION__, input_vcf_files.size());
        input_vcf_files.insert(input_vcf_files.begin(), candidate_sites_vcf_file);
        sr = new BCFSyncedReader(input_vcf_files, intervals, false);
        fprintf(stderr, " done.\n");

        odw = new BCFOrderedWriter(output_vcf_file, 0);
        bcf_hdr_append(odw->hdr, "##fileformat=VCFv4.2");
        bcf_hdr_transfer_contigs(sr->hdrs[0], odw->hdr);
        bool rename = true;

//        odw->link_hdr(sr->hdrs[0]);
//        //exact alignment related statisitcs
        std::string EX_MOTIF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_MOTIF", "1", "String", "Canonical motif in a VNTR.", rename);
        std::string EX_RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_RU", "1", "String", "Repeat unit in the reference sequence.", rename);
        std::string EX_BASIS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_BASIS", "1", "String", "Basis nucleotides in the motif.", rename);
        std::string EX_MLEN = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_MLEN", "1", "Integer", "Motif length.", rename);
        std::string EX_BLEN = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_BLEN", "1", "Integer", "Basis length.", rename);
        std::string EX_REPEAT_TRACT = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_REPEAT_TRACT", "2", "Integer", "Boundary of the repeat tract detected by exact alignment.", rename);
        std::string EX_COMP = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_COMP", "4", "Integer", "Composition(%) of bases in an exact repeat tract.", rename);
        std::string EX_ENTROPY = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_ENTROPY", "1", "Float", "Entropy measure of an exact repeat tract [0,2].", rename);
        std::string EX_ENTROPY2 = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_ENTROPY2", "1", "Float", "Dinucleotide entropy measure of an exact repeat tract [0,4].", rename);
        std::string EX_KL_DIVERGENCE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_KL_DIVERGENCE", "1", "Float", "Kullback-Leibler Divergence of an exact repeat tract.", rename);
        std::string EX_KL_DIVERGENCE2 = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_KL_DIVERGENCE2", "1", "Float", "Dinucleotide Kullback-Leibler Divergence of an exact repeat tract.", rename);
        std::string EX_REF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_REF", ".", "Float", "Allele lengths in repeat units from exact alignment.", rename); 
        std::string EX_RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_RL", "1", "Integer", "Reference exact repeat tract length in bases.", rename);
        std::string EX_LL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_LL", "1", "Integer", "Longest exact repeat tract length in bases.", rename);
        std::string EX_RU_COUNTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_RU_COUNTS", "2", "Integer", "Number of exact repeat units and total number of repeat units in exact repeat tract.", rename);
        std::string EX_SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_SCORE", "1", "Float", "Score of repeat unit in exact repeat tract.", rename);
        std::string EX_TRF_SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EX_TRF_SCORE", "1", "Integer", "TRF Score for M/I/D as 2/-7/-7 in exact repeat tract.", rename);
        
        //fuzzy alignment related statisitcs
        std::string FZ_MOTIF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_MOTIF", "1", "String", "Canonical motif in a VNTR.", rename);
        std::string FZ_RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_RU", "1", "String", "Repeat unit in the reference sequence.", rename);
        std::string FZ_BASIS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_BASIS", "1", "String", "Basis nucleotides in the motif.", rename);
        std::string FZ_MLEN = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_MLEN", "1", "Integer", "Motif length.", rename);
        std::string FZ_BLEN = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_BLEN", "1", "Integer", "Basis length.", rename);
        std::string FZ_REPEAT_TRACT = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_REPEAT_TRACT", "2", "Integer", "Boundary of the repeat tract detected by fuzzy alignment.", rename);                      
        std::string FZ_COMP = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_COMP", "4", "Integer", "Composition(%) of bases in a fuzzy repeat tract.", rename);                                              
        std::string FZ_ENTROPY = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_ENTROPY", "1", "Float", "Entropy measure of a fuzzy repeat tract (0-2).", rename);                                                  
        std::string FZ_ENTROPY2 = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_ENTROPY2", "1", "Float", "Dinucleotide entropy measure of a fuzzy repeat tract (0-2).", rename);                                                  
        std::string FZ_KL_DIVERGENCE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_KL_DIVERGENCE", "1", "Float", "Kullback-Leibler Divergence of a fuzzyt repeat tract.", rename);
        std::string FZ_KL_DIVERGENCE2 = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_KL_DIVERGENCE2", "1", "Float", "Dinucleotide Kullback-Leibler Divergence of a fuzzy repeat tract.", rename);
        std::string FZ_REF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_REF", ".", "Float", "Allele lengths in repeat units from fuzzy alignment.", rename);                                                      
        std::string FZ_RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_RL", "1", "Integer", "Reference fuzzy repeat tract length in bases.", rename);                                                      
        std::string FZ_LL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_LL", "1", "Integer", "Longest fuzzy repeat tract length in bases.", rename);                                                        
        std::string FZ_RU_COUNTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_RU_COUNTS", "2", "Integer", "Number of exact repeat units and total number of repeat units in fuzzy repeat tract.", rename); 
        std::string FZ_SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_SCORE", "1", "Float", "Score of repeat unit in fuzzy repeat tract.", rename);                                                    
        std::string FZ_TRF_SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_TRF_SCORE", "1", "Integer", "TRF Score for M/I/D as 2/-7/-7 in fuzzy repeat tract.", rename);                                 
                                                                                                                                                                                                         
        std::string FLANKSEQ = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FLANKSEQ", "1", "String", "Flanking sequence 10bp on either side of REF.", rename);
        std::string EXACT_RU_AMBIGUOUS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "EXACT_RU_AMBIGUOUS", "0", "Flag", "Exact motif is ambiguous.", rename);

        std::string MOTIF = bcf_hdr_append_info_with_backup_naming(odw->hdr, "MOTIF", "1", "String", "Canonical motif in a VNTR.", rename);
        std::string RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU", "1", "String", "Repeat unit in the reference sequence.", rename);
        std::string FLANKS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FLANKS", "2", "Integer", "Exact left and right flank positions of the Indel.", rename);
        std::string FZ_FLANKS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FZ_FLANKS", "2", "Integer", "Fuzzy left and right flank positions of the Indel.", rename);
       
       
        bcf_hdr_append(odw->hdr, "##INFO=<ID=LARGE_REPEAT_REGION,Number=0,Type=Flag,Description=\"Very large repeat region, vt only detects up to 1000bp long regions.\">");
        bcf_hdr_append(odw->hdr, "##INFO=<ID=FLANKSEQ,Number=1,Type=String,Description=\"Flanking sequence 10bp on either side of detected repeat region.\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_vntr,Description=\"Overlaps with VNTR\">");
        bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_indel,Description=\"Overlaps with indel\">");

        //COMMON
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED scaled genotype likelihoods\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depth\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=ADF,Number=A,Type=Integer,Description=\"Allele Depth (Forward strand)\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=ADR,Number=A,Type=Integer,Description=\"Allele Depth (Reverse strand)\">");
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");

        //VNTR
        bcf_hdr_append(odw->hdr, "##FORMAT=<ID=CG,Number=.,Type=Float,Description=\"Repeat count genotype\">");

        odw->write_hdr();

        //add samples to output merged file
        for (int32_t i=1; i<sr->hdrs.size(); ++i)
        {
//            printf("adding sample = %s\n",  bcf_hdr_get_sample_name(sr->hdrs[i], 0));
            bcf_hdr_add_sample(odw->hdr, bcf_hdr_get_sample_name(sr->hdrs[i], 0) );
        }
        bcf_hdr_sync(odw->hdr);
       
        ///////////////
        //general use//
        ///////////////
        variant = {0,0,0};
        no_samples = sr->nfiles - 1;

        ////////////////////////
        //stats initialization//
        ////////////////////////
        no_snps = 0;
        no_indels = 0;
        no_vntrs = 0;

        /////////
        //tools//
        /////////
        vm = new VariantManip();
    }

    void merge_genotypes()
    {
        int32_t *NSAMPLES = NULL;
        int32_t no_NSAMPLES = 0;
        int32_t *GT = NULL;
        int32_t *PL = NULL;
        int32_t *AD = NULL;
        int32_t *ADF = NULL;
        int32_t *ADR = NULL;
        int32_t *BQSUM = NULL;
        int32_t *DP = NULL;
        float *CG = NULL;
        
        int32_t no_GT = 0;
        int32_t no_PL = 0;
        int32_t no_AD = 0;
        int32_t no_ADF = 0;
        int32_t no_ADR = 0;
        int32_t no_BQSUM = 0;
        int32_t no_DP = 0;
        int32_t no_CG = 0;

        std::vector<int32_t> gt;
        std::vector<int32_t> pl;
        std::vector<int32_t> ad;
        std::vector<int32_t> dp;
        std::vector<float> cg;

        bcf1_t* nv = bcf_init();
        Variant var;
        std::vector<bcfptr*> current_recs;
                
        while(sr->read_next_position(current_recs))
        {
            if (current_recs.size()!=no_samples+1)
            {
                fprintf(stderr, "[E:%s:%d %s] %d variants expected but %zd is observed.\n", __FILE__, __LINE__, __FUNCTION__, no_samples+1, current_recs.size());
                exit(1);
            }
            
            gt.resize(0);
            pl.resize(0);
            ad.resize(0);
            dp.resize(0);
            cg.resize(0);

            bcf_clear(nv);
            int32_t vtype;

            if (filter_exists)
            {
                vm->classify_variant(current_recs[0]->h, current_recs[0]->v, var);
                if (!filter.apply(current_recs[0]->h, current_recs[0]->v, &var, false))
                {
                    continue;
                }
            }
            
            if (current_recs.size())

            //for each file
            for (uint32_t i=0; i<current_recs.size(); ++i)
            {
                int32_t file_index = current_recs[i]->file_index;
                bcf1_t *v = current_recs[i]->v;
                bcf_hdr_t *h = current_recs[i]->h;


//                printf("\tfile index: %d\n", file_index);
//                bcf_print(h, v);
               

                //candidate sites file, populate info fields
                if (!file_index) 
                {
                    vtype = vm->classify_variant(h, v, var); 
                    
//                    bcf_copy(v, nv);
                    
                    bcf_set_chrom(odw->hdr, nv, bcf_get_chrom(h, v));
                    bcf_set_pos1(nv, bcf_get_pos1(v));
                    bcf_update_alleles(odw->hdr, nv, const_cast<const char**>(bcf_get_allele(v)), bcf_get_n_allele(v));
                    bcf_set_n_sample(nv, no_samples);
    
                    if (vtype==VT_SNP)
                    {
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FLANKSEQ", BCF_HT_STR);
                    }
                    else if (vtype==VT_INDEL)
                    {
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_MOTIF", BCF_HT_STR);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_MLEN", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_RU", BCF_HT_STR);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_BASIS", BCF_HT_STR);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_BLEN", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_REPEAT_TRACT", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_COMP", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_ENTROPY", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_ENTROPY2", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_KL_DIVERGENCE", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_KL_DIVERGENCE2", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_REF", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_RL", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_LL", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_RU_COUNTS", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_SCORE", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EX_TRF_SCORE", BCF_HT_INT);
                        
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_MOTIF", BCF_HT_STR);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_MLEN", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_RU", BCF_HT_STR);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_BASIS", BCF_HT_STR);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_BLEN", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_REPEAT_TRACT", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_COMP", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_ENTROPY", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_ENTROPY2", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_KL_DIVERGENCE", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_KL_DIVERGENCE2", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_REF", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_RL", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_LL", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_RU_COUNTS", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_SCORE", BCF_HT_REAL);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_TRF_SCORE", BCF_HT_INT);
                   
                        bcf_copy_info_field(h, v, odw->hdr, nv, "EXACT_RU_AMBIGUOUS", BCF_HT_FLAG);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "LARGE_REPEAT_REGION", BCF_HT_FLAG);
                    }
                    else if (vtype==VT_VNTR)
                    {
//                        printf("\t\tis a VNTR\n");
//                        printf("\t\t\tcopying info fields\n");
//                        bcf_print(h, v);
                        
//                        bcf_copy_info_field(h, v, odw->hdr, nv, "MOTIF", BCF_HT_STR);
//                        bcf_copy_info_field(h, v, odw->hdr, nv, "RU", BCF_HT_STR);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_RL", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FLANKS", BCF_HT_INT);
                        bcf_copy_info_field(h, v, odw->hdr, nv, "FZ_FLANKS", BCF_HT_INT);
                    }
                    
                    if (bcf_has_filter(h, v, const_cast<char*>("overlap_indel"))==1)
                    {
                        int32_t overlap_indel_filter_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, const_cast<char*>("overlap_indel"));
                        bcf_update_filter(odw->hdr, nv, &overlap_indel_filter_id, 1);
                    }
                    
                    if (bcf_has_filter(h, v, const_cast<char*>("overlap_vntr"))==1)
                    {
                        int32_t overlap_indel_filter_id = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, const_cast<char*>("overlap_vntr"));
                        bcf_update_filter(odw->hdr, nv, &overlap_indel_filter_id, 1);
                    }
                    
                    continue;
                }
                                
                if (vtype==VT_SNP)
                {
//                    printf("\t\tis a SNP\n");
                    
                    int32_t no_gt = bcf_get_genotypes(h, v, &GT, &no_GT);
                    int32_t no_pl = bcf_get_format_int32(h, v, "PL", &PL, &no_PL);
                    int32_t no_dp = bcf_get_format_int32(h, v, "DP", &DP, &no_DP); 
                    int32_t no_adf = bcf_get_format_int32(h, v, "ADF", &ADF, &no_ADF); 
                    int32_t no_adr = bcf_get_format_int32(h, v, "ADR", &ADR, &no_ADR);
                    int32_t no_bqsum = bcf_get_format_int32(h, v, "BQSUM", &BQSUM, &no_BQSUM); 
                    
//                    printf("GT: %d\n", no_gt);
//                    printf("PL: %d\n", no_pl);
//                    printf("DP: %d\n", no_dp);
//                    printf("ADF: %d\n", no_adf);
//                    printf("ADR: %d\n", no_adr);
//                    printf("BQSUM: %d\n", no_bqsum);
                          
                    //GT:PL:DP:AD:ADF:ADR:BQ:MQ:CY:ST:AL:NM
                    if (no_gt > 0 &&
                        no_pl > 0 &&
                        no_dp > 0 &&
                        no_adf > 0 &&
                        no_adr > 0)
                    {
//                        printf("\t\tupdating genotypes for GT,PL\n");
//                        printf("\t\t\tno_GT %d\n", no_GT);
//                        printf("\t\t\t%d %d\n", GT[0], GT[1]);
                          gt.push_back(GT[0]);
                          gt.push_back(GT[1]);
                          pl.push_back(PL[0]);
                          pl.push_back(PL[1]);
                          pl.push_back(PL[2]);
                          dp.push_back(DP[0]);
                          ad.push_back(ADF[0]+ADR[0]);
                          ad.push_back(ADF[1]+ADR[1]);

//                          printf("GT: %d\n", no_GT);
//                          printf("PL: %d\n", no_PL);
//                          printf("DP: %d\n", no_DP);
//                          printf("ADF: %d\n", no_ADF);
//                          printf("ADR: %d\n", no_ADR);

                    }
                    //GT:BQSUM:DP
                    else if (no_gt > 0 &&
                             no_bqsum > 0 &&
                             no_dp > 0)
                    {
//                          printf("\t\tupdating genotypes for GT,BQSUM\n");
//                          printf("\t\t\tno_GT %d\n", no_GT);
//                          printf("\t\t\t%d %d\n", GT[0], GT[1]);
                          gt.push_back(GT[0]);
                          gt.push_back(GT[1]);
                          pl.push_back(0);
                          pl.push_back(BQSUM[0]/3);
                          pl.push_back(BQSUM[0]);
                          dp.push_back(DP[0]);
                          ad.push_back(DP[0]);
                          ad.push_back(0);                          
                    }
                    else
                    {
                        fprintf(stderr, "[E:%s:%d %s] cannot get format values GT:PL:DP:ADF:ADR or GT:BQSUM:DP from %s\n", __FILE__, __LINE__, __FUNCTION__, sr->file_names[file_index].c_str());
                        exit(1);
                    }
                }
                else if (vtype == VT_INDEL)
                {
//                    printf("\t\tis an INDEL\n");
                    
                    int32_t no_gt = bcf_get_genotypes(h, v, &GT, &no_GT);
                    int32_t no_pl = bcf_get_format_int32(h, v, "PL", &PL, &no_PL);
                    int32_t no_dp = bcf_get_format_int32(h, v, "DP", &DP, &no_DP); 
                    int32_t no_adf = bcf_get_format_int32(h, v, "ADF", &ADF, &no_ADF); 
                    int32_t no_adr = bcf_get_format_int32(h, v, "ADR", &ADR, &no_ADR);
                    int32_t no_bqsum = bcf_get_format_int32(h, v, "BQSUM", &BQSUM, &no_BQSUM); 
                    
                    if (no_gt > 0 &&
                        no_pl > 0 &&
                        no_dp > 0 &&
                        no_adf > 0 &&
                        no_adr > 0)
                    {
//                        printf("\t\tupdating genotypes for GT,PL\n");
//                        printf("\t\t\tno_GT %d\n", no_GT);
//                        printf("\t\t\t%d %d\n", GT[0], GT[1]);
                          gt.push_back(GT[0]);
                          gt.push_back(GT[1]);
                          pl.push_back(PL[0]);
                          pl.push_back(PL[1]);
                          pl.push_back(PL[2]);
                          dp.push_back(DP[0]);
                          ad.push_back(ADF[0]+ADR[0]);
                          ad.push_back(ADF[1]+ADR[1]);

//                          printf("GT: %d\n", no_GT);
//                          printf("PL: %d\n", no_PL);
//                          printf("DP: %d\n", no_DP);
//                          printf("ADF: %d\n", no_ADF);
//                          printf("ADR: %d\n", no_ADR);

                    }
                    //GT:BQSUM:DP
                    else if (no_gt > 0 &&
                             no_bqsum > 0 &&
                             no_dp > 0)
                    {
//                          printf("\t\tupdating genotypes for GT,BQSUM\n");
//                          printf("\t\t\tno_GT %d\n", no_GT);
//                          printf("\t\t\t%d %d\n", GT[0], GT[1]);
                          gt.push_back(GT[0]);
                          gt.push_back(GT[1]);
                          pl.push_back(0);
                          pl.push_back(BQSUM[0]/3);
                          pl.push_back(BQSUM[0]);
                          dp.push_back(DP[0]);
                          ad.push_back(DP[0]);
                          ad.push_back(0);                          
                    }
                    else
                    {
                        fprintf(stderr, "[E:%s:%d %s] cannot get format values GT:PL:DP:ADF:ADR or GT:BQSUM:DP from %s\n", __FILE__, __LINE__, __FUNCTION__, sr->file_names[file_index].c_str());
                        exit(1);
                    }
                }
                else if (vtype == VT_VNTR)
                {
//                    printf("\t\tis a VNTR\n");
//                    bcf_print(h, v);
//                    
                    int32_t no_cg = bcf_get_format_float(h, v, "CG", &CG, &no_CG); 
                    
                    //CG
                    if (no_cg > 0)
                    {
//                          printf("\t\tupdating genotypes for CG\n");
//                          printf("\t\t\tno_CG %d\n", no_CG);
//                          printf("\t\t\t%f %f\n", CG[0], CG[1]);
                          cg.push_back(CG[0]);
                          cg.push_back(CG[1]);
                    }
                    else
                    {
                        fprintf(stderr, "[E:%s:%d %s] cannot get format values CG from %s\n", __FILE__, __LINE__, __FUNCTION__, sr->file_names[file_index].c_str());
                    }
                }

            }//end processing each file

            //write to merged record
            if (vtype==VT_SNP)
            {
//                printf("\tupdating genotypes %zd\n", gt.size());
//                printf("\tn_samples %zd\n", nv->n_sample);
                bcf_update_genotypes(odw->hdr, nv, &gt[0], gt.size());
                bcf_update_format_int32(odw->hdr, nv, "PL", &pl[0], pl.size());
                bcf_update_format_int32(odw->hdr, nv, "DP", &dp[0], dp.size());
                bcf_update_format_int32(odw->hdr, nv, "AD", &ad[0], ad.size());
            
                ++no_snps;
            }
            else if (vtype==VT_INDEL)
            {
//                printf("\tupdating genotypes %zd\n", gt.size());
//                printf("\tn_samples %zd\n", nv->n_sample);
                bcf_update_genotypes(odw->hdr, nv, &gt[0], gt.size());
                bcf_update_format_int32(odw->hdr, nv, "PL", &pl[0], pl.size());
                bcf_update_format_int32(odw->hdr, nv, "DP", &dp[0], dp.size());
                bcf_update_format_int32(odw->hdr, nv, "AD", &ad[0], ad.size());
                
                ++no_indels;
            }
            else if (vtype==VT_VNTR)
            {
//                printf("\tupdating genotypes %zd\n", gt.size());
//                printf("\tn_samples %zd\n", nv->n_sample);
                bcf_update_format_float(odw->hdr, nv, "CG", &cg[0], cg.size());
                
                ++no_vntrs;
                
//                bcf_print(odw->hdr, nv);
            }
            
            odw->write(nv);
//            bcf_print(odw->hdr, nv);

            //this acts as a flag to initialize a newly merged record
            vtype = VT_UNDEFINED;
            
            int32_t no_variants = no_snps+no_indels+no_vntrs;
//            if ((no_variants&0x0FFF)==0x0600)
            if ((no_variants%100)==0)
            {
                fprintf(stderr, "[I:%s:%d %s] Merged %d rows\n", __FILE__, __LINE__, __FUNCTION__, no_variants);
            }    
            
        }

        sr->close();
        odw->close();
        
        
    };

    void print_options()
    {
        std::clog << "merge_genotypes v" << version << "\n\n";
        std::clog << "options: [L] input VCF file list                      " << input_vcf_file_list << " (" << input_vcf_files.size() << " files)\n";
        std::clog << "         [c] candidate sites VCF file                 " << candidate_sites_vcf_file << "\n";
        std::clog << "         [o] output VCF file                          " << output_vcf_file << "\n";
        print_str_op("         [f] filter (applied to candidate sites file) ", fexp);
        print_int_op("         [i] intervals                                ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::clog << "stats: Total no. of SNPs                 " << no_snps << "\n";
        std::clog << "       Total no. of Indels               " << no_indels << "\n";
        std::clog << "       Total no. of VNTRs                " << no_vntrs << "\n";
        std::clog << "\n";
    };

    ~Igor()
    {
    };

    private:
};

}

void merge_genotypes(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.merge_genotypes();
    igor.print_stats();
}

