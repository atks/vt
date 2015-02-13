OPTFLAG = -O3 -ggdb
INCLUDES = -I./lib/include/ -I. -I./lib/include/htslib -I./lib/include/Rmath 
CFLAGS = -pipe -std=c++0x $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = g++

SOURCESONLY = genotyping_record.h

SOURCES = align\
		allele\
		annotate_dbsnp_rsid\
		annotate_indel\
		annotate_regions\
		annotate_variants\
		bam_ordered_reader\
		bcf_ordered_reader\
		bcf_ordered_writer\
		bcf_synced_reader\
		bed\
		candidate_motif\
		cat\
		chmm\
		compute_concordance\
		compute_features\
		config\
		consolidate_variants\
		construct_probes\
		context_filter\
		decompose\
		decompose_blocksub\
		discover\
		discover2\
		estimate\
		estimator\
		filter\
		gencode\
		genome_interval\
		genotype\
		genotype2\
		genotyping_buffer\
		hts_utils\
		index\
		interval_tree\
		interval\
		lfhmm\
		lhmm\
		lhmm1\
		lhmm_genotyping_record\
		log_tool\
		merge\
		merge_candidate_variants\
		motif_suffix_tree\
		ordered_bcf_overlap_matcher\
		ordered_region_overlap_matcher\
		partition\
		paste\
		pedigree\
		peek\
		pileup\
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
		program\
		remove_overlap\
		rfhmm\
		seq\
		sort\
		str\
		subset\
		sv_tree\
		test\
		union_variants\
		uniq\
		utils\
		validate\
		variant\
		variant_manip\
		view\
		vntrize\
   	    tbx_ordered_reader\
    	ahmm\
        xcmp\
        normalize\

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
