OPTFLAG = -O3
INCLUDES = -I./lib -I. -I./lib/htslib -I./lib/Rmath -I./lib/pcre2
CFLAGS = -pipe -std=c++0x $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = g++

SOURCESONLY =

SOURCES = align\
		allele\
		annotate_1000g\
		annotate_dbsnp_rsid\
		annotate_indels\
		annotate_regions\
		annotate_variants\
		augmented_bam_record\
		bam_ordered_reader\
		bcf_genotyping_buffered_reader\
		bcf_ordered_reader\
		bcf_ordered_writer\
		bcf_synced_reader\
		bed\
		candidate_motif_picker\
		candidate_region_extractor\
		cat\
		chmm\
		compute_concordance\
		compute_features\
		compute_features2\
		config\
		consolidate_vntrs\
		consolidate\
		construct_probes\
		decompose\
		decompose2\
		decompose_blocksub\
		discover\
		estimate\
		estimator\
		filter\
		flank_detector\
		gencode\
		genome_interval\
		genotype\
		genotype2\
		genotyping_record\
		hts_utils\
		hfilter\
		index\
		interval_tree\
		interval\
		lfhmm\
		lhmm\
		lhmm1\
		log_tool\
		merge\
		merge_candidate_variants\
		merge_candidate_variants2\
		motif_tree\
		motif_map\
		multi_partition\
		needle\
		normalize\
		ordered_bcf_overlap_matcher\
		ordered_region_overlap_matcher\
		partition\
		paste\
		paste_and_compute_features_sequential\
		pedigree\
		peek\
		pileup\
		pregex\
		profile_afs\
		profile_chm1\
		profile_chrom\
		profile_fic_hwe\
		profile_hwe\
		profile_indels\
		profile_len\
		profile_mendelian\
		profile_na12878\
		profile_snps\
		profile_vntrs\
		program\
		reference_sequence\
		remove_overlap\
		repeat_tract_detector\
		rfhmm\
		rminfo\
		rpartition\
		seq\
		sort\
		subset\
		sv_tree\
		svm_train\
		svm_predict\
		test\
		union_variants\
		uniq\
		utils\
		validate\
		variant\
		variant_manip\
		variant_filter\
		view\
		vntr\
		vntr_annotator\
		vntr_tree\
		vntrize\
   	    tbx_ordered_reader\
    	ahmm\

SOURCESONLY = main.cpp

TARGET = vt
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
LIBHTS = lib/htslib/libhts.a
LIBRMATH = lib/Rmath/libRmath.a
LIBPCRE2 = lib/pcre2/libpcre2.a
LIBSVM = lib/libsvm/libsvm.a

all : $(TARGET)

${LIBHTS} :
	cd lib/htslib; $(MAKE) libhts.a || exit 1; 

${LIBRMATH} :
	cd lib/Rmath; $(MAKE) libRmath.a || exit 1; 

${LIBPCRE2} :
	cd lib/pcre2; $(MAKE) libpcre2.a || exit 1; 

${LIBSVM} :
	cd lib/libsvm; $(MAKE) libsvm.a || exit 1; 
	
$(TARGET) : ${LIBHTS} ${LIBRMATH} ${LIBPCRE2}  ${LIBSVM} $(TOOLOBJ)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBHTS) $(LIBRMATH) ${LIBPCRE2} -lz -lpthread

$(TOOLOBJ): $(HEADERSONLY)

.cpp.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp

.PHONY: clean cleanvt test

clean :
	cd lib/htslib; $(MAKE) clean
	cd lib/Rmath; $(MAKE) clean
	cd lib/pcre2; $(MAKE) clean
	cd lib/libsvm; $(MAKE) clean
	-rm -rf $(TARGET) $(TOOLOBJ)

cleanvt :
	-rm -rf $(TARGET) $(TOOLOBJ)    

test : vt
	test/test.sh
	
debug : vt
	test/test.sh debug
