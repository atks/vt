OPTFLAG ?= -O3 -ggdb
INCLUDES = -I./lib/include/ -I.
CFLAGS = -pipe -std=c++0x -Wall $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = g++

HEADERSONLY =
SOURCES = filter\
          compute_concordance\
          program\
          ordered_reader\
          ordered_writer\
          normalize\
          synced_reader\
          merge_duplicate_variants\
          hts_utils\
          variant_manip\
          log_tool\
          interval_tree\
          rb_tree
SOURCESONLY = main.cpp

TARGET = vt
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
LIBHTS = ./lib/include/htslib/libhts.a

all : $(TARGET)

$(TARGET) : $(TOOLOBJ)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBHTS) -lz -lpthread

#$(TOOLOBJ): $(HEADERSONLY)

.cpp.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp

clean:
	-rm -rf $(TARGET) $(TOOLOBJ)

