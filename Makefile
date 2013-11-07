UNAME = $(shell uname)
STD = 
ifeq ($(UNAME), Linux)
	STD=-std=c++0x
endif

	
OPTFLAG ?= -O0 -ggdb
INCLUDES = -I./lib/include/ -I. -I./lib/include/htslib
CFLAGS = -pipe $(STD) $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = g++

HEADERSONLY =
SOURCES = program\
		variant_manip\
		log_tool\
		interval_tree\
		genome_interval\
		rb_tree\
		hts_utils\
		utils\
		ordered_reader\
		ordered_writer\
		synced_reader\
		merge_duplicate_variants\
		normalize\
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
