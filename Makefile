UNAME = $(shell uname)

OPTFLAG ?= -O0 -ggdb
INCLUDES = -I./lib/include/ -I. -I./lib/include/htslib
CFLAGS = -pipe -std=c++0x $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = clang++

HEADERSONLY = 
SOURCES = program\
		filter\
		hts_utils\
		utils\
		bam_ordered_reader\
		bcf_ordered_reader\
		bcf_ordered_writer\
		bcf_synced_reader\
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
		partition\
		profile_indels\
		discover\
		merge_candidate_variants\
		construct_probes\
		genotype\
		gencode\
		annotate_variants\
		lhmm\
		genotyping_buffer\
		lhmm_genotyping_record\
		peek\
		merge\
		profile_mendel_errors\
		pedigree\
		concat

SOURCESONLY = main.cpp

TARGET = vt
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
LIBHTS = lib/include/htslib/libhts.a

all : ${LIBHTS} $(TARGET)

${LIBHTS} :
	cd lib/include/htslib; $(MAKE) libhts.a || exit 1; cd ..

$(TARGET) : ${LIBHTS} $(TOOLOBJ)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBHTS) -lz -lpthread

$(TOOLOBJ): $(HEADERSONLY)

.cpp.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp

clean :
	cd lib/include/htslib; $(MAKE) clean; cd ..
	-rm -rf $(TARGET) $(TOOLOBJ)

cleanvt :
	-rm -rf $(TARGET) $(TOOLOBJ)    