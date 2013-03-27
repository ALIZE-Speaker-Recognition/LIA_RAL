bin_PROGRAMS=LabelNGram
LabelNGram_SOURCES= \
	LabelNGram.cpp\
	LabelNGramMain.cpp

LabelNGram_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
LabelNGram_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
LabelNGram_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
