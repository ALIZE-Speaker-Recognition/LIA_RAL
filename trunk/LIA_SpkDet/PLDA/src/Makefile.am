bin_PROGRAMS=PLDA
PLDA_SOURCES= \
	PLDA.cpp\
	PLDAMain.cpp

PLDA_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen -I$(LAPACK_PATH)/lapacke/include
PLDA_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib -L$(top_srcdir)/include/Pthread -L$(LAPACK_PATH)
PLDA_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD) $(LIBLAPACK)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
