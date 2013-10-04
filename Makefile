OPTFLAG ?= -O3 -ggdb
INCLUDES = -I./lib/include/ -I.
CFLAGS = -pipe -std=c++0x -Wall $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = g++
CC = gcc

HEADERSONLY =
SOURCES = program\
		filter\
		variant_manip\
		log_tool\
		interval_tree\
		rb_tree\
		hts_utils\
		ordered_reader\
		ordered_writer\
		synced_reader\
		profile_indels\
		merge_duplicate_variants\
		compute_concordance\
		normalize\
		partition

SOURCESONLY = main.cpp

TARGET = vt
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
LIBHTS = lib/include/htslib/libhts.a

all : ${LIBHTS} $(TARGET)

${LIBHTS} : 
	cd lib/include/htslib; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhts.a || exit 1; cd ..

$(TARGET) : ${LIBHTS} $(TOOLOBJ)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBHTS) -lz -lpthread

$(TOOLOBJ): $(HEADERSONLY)

.cpp.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp

clean:
	-rm -rf $(TARGET) $(TOOLOBJ)

