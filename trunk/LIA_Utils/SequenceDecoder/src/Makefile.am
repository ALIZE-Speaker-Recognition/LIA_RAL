bin_PROGRAMS=SequenceDecoder
SequenceDecoder_SOURCES= \
	SequenceDecoderMain.cpp\
	../../SequenceExtractor/src/SequenceExtractor.cpp\
	../../BNGram/src/BNGram.cpp\
	../../LabelNGram/src/LabelNGram.cpp


SequenceDecoder_CPPFLAGS=-I../include -I../../SequenceExtractor/include -I../../LabelNGram/include -I../../BNGram/include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
SequenceDecoder_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
SequenceDecoder_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
