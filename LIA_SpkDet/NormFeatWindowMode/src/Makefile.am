bin_PROGRAMS=NormFeatWindowMode
NormFeatWindowMode_SOURCES= \
	NormFeatWindowMode.cpp\
	NormFeatWindowModeMain.cpp

NormFeatWindowMode_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
NormFeatWindowMode_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib -L$(top_srcdir)/include/Pthread
NormFeatWindowMode_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin