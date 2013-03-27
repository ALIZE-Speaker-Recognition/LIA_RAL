bin_PROGRAMS=GmmTokenizer
GmmTokenizer_SOURCES= \
	GmmTokenizer.cpp\
	GmmTokenizerMain.cpp

GmmTokenizer_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
GmmTokenizer_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
GmmTokenizer_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
