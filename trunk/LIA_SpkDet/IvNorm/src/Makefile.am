bin_PROGRAMS=IvNorm
IvNorm_SOURCES= \
	IvNorm.cpp\
	IvNormMain.cpp

IvNorm_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
IvNorm_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib -L$(top_srcdir)/include/Pthread -L$(LAPACK_PATH)
IvNorm_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD) $(LIBLAPACK)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
