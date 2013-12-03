UNAME = $(shell uname)

OPTFLAG ?= -O0 -ggdb
INCLUDES = -I./lib/include/ -I. -I./lib/include/htslib
CFLAGS = -pipe -std=c++0x $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = clang++

HEADERSONLY =
SOURCES = program\
		variant_manip\
		log_tool\
		interval_tree\
		genome_interval\
		filter\
		rb_tree\
		lhmm\
		hts_utils\
		utils\
		bam_ordered_reader\
		bcf_ordered_reader\
		bcf_ordered_writer\
		bcf_synced_reader\
		view\
		normalize\
		merge_duplicate_variants\
		construct_probes\
		discover\
		genotype\
		profile_indels\
		compute_concordance\
		partition\
		peek

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

clean:
	-rm -rf $(TARGET) $(TOOLOBJ)
