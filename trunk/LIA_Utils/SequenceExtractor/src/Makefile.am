bin_PROGRAMS=SequenceExtractor
SequenceExtractor_SOURCES= \
	SequenceExtractor.cpp\
	SequenceExtractorMain.cpp\
	../../BNGram/src/BNGram.cpp\
	../../LabelNGram/src/LabelNGram.cpp

SequenceExtractor_CPPFLAGS=-I../include -I../../BNGram/include -I../../LabelNGram/include/ -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
SequenceExtractor_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
SequenceExtractor_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
