SUBDIRS=SVDLIBC
bin_PROGRAMS=CovIntra
CovIntra_SOURCES= \
	CovIntra.cpp\
	CovIntraMain.cpp

CovIntra_CPPFLAGS=-ISVDLIBC/ -I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
CovIntra_LDFLAGS=-L./SVDLIBC/ -L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib -L$(top_srcdir)/include/Pthread
CovIntra_LDADD=-lsvd -lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin