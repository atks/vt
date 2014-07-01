OPTFLAG = -O3 -ggdb
INCLUDES = -I./lib/include/ -I. -I./lib/include/htslib -I./lib/include/Rmath 
CFLAGS = -pipe -std=c++0x $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = g++

SOURCESONLY = genotyping_record.h

SOURCES = program\
		hts_utils\
		utils\
		bam_ordered_reader\
		bcf_ordered_reader\
		bcf_ordered_writer\
		bcf_synced_reader\
		bcf_synced_sreader\
   	    tbx_ordered_reader\
		view\
		index\
        normalize\
		merge_duplicate_variants\
		variant_manip\
		log_tool\
		interval_tree\
		genome_interval\
		compute_concordance\
		compute_features\
		partition\
		profile_indels\
		profile_snps\
		profile_na12878\
		discover\
		merge_candidate_variants\
		construct_probes\
		genotype\
		gencode\
		annotate_variants\
		annotate_str\
		lhmm\
		chmm\
		lfhmm\
		rfhmm\
		lhmm1\
		genotyping_buffer\
		lhmm_genotyping_record\
		peek\
		merge\
		profile_mendelian\
		pedigree\
		concat\
		decompose\
		remove_overlap\
		filter\
		align\
		subset\
		estimator\
		profile_afs\
		profile_hwe\
		profile_len\
		ordered_region_overlap_matcher\
		ordered_bcf_overlap_matcher\
		bed\
		profile_chrom\
		annotate_regions\
		test\
		str\
		consolidate_variants\
		genotype2\
		paste\
		annotate_dbsnp_rsid\
		config

SOURCESONLY = main.cpp

TARGET = vt
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
LIBHTS = lib/include/htslib/libhts.a
LIBRMATH = lib/include/Rmath/libRmath.a

all : ${LIBHTS} $(TARGET)

${LIBHTS} :
	cd lib/include/htslib; $(MAKE) libhts.a || exit 1; cd ..

${LIBRMATH} :
	cd lib/include/Rmath; $(MAKE) libRmath.a || exit 1; cd ..
	
$(TARGET) : ${LIBHTS} ${LIBRMATH} $(TOOLOBJ)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBHTS) $(LIBRMATH) -lz -lpthread

$(TOOLOBJ): $(HEADERSONLY)

.cpp.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp

clean :
	cd lib/include/htslib; $(MAKE) clean
	cd lib/include/Rmath; $(MAKE) clean
	-rm -rf $(TARGET) $(TOOLOBJ)

cleanvt :
	-rm -rf $(TARGET) $(TOOLOBJ)    
