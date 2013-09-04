bin_PROGRAMS=ReSegmentation
ReSegmentation_SOURCES= \
	ReSegmentation.cpp\
	ReSegmentationMain.cpp

ReSegmentation_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include  -I$(top_srcdir)/include/Eigen
ReSegmentation_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
ReSegmentation_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
