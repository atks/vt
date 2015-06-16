AR     = ar
RANLIB = ranlib

SOURCES = pcre2_auto_possess\
       pcre2_chartables\
       pcre2_compile\
       pcre2_config\
       pcre2_context\
       pcre2_dfa_match\
       pcre2_error\
       pcre2_jit_compile\
       pcre2_maketables\
       pcre2_match\
       pcre2_match_data\
       pcre2_newline\
       pcre2_ord2utf\
       pcre2_pattern_info\
       pcre2_string_utils\
       pcre2_study\
       pcre2_substitute\
       pcre2_substring\
       pcre2_tables\
       pcre2_ucd\
       pcre2_valid_utf\
       pcre2_xclass\

TOOLSRC = $(SOURCES:=.c)
TOOLOBJ = $(TOOLSRC:.c=.o)

.c.o :
	gcc -g -O3 -DHAVE_CONFIG_H -DPCRE2_CODE_UNIT_WIDTH=8 -I. -o $@ -c $*.c

libpcre2.a: $(TOOLOBJ) 
	@-rm -f $@
	ar -rc $@ $(TOOLOBJ)
	-ranlib $@
	
clean :
	-rm *.o *.a

