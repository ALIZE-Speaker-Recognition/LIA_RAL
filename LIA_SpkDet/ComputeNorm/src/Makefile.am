bin_PROGRAMS=ComputeNorm
ComputeNorm_SOURCES= \
	ComputeNorm.cpp\
	ComputeNormMain.cpp

ComputeNorm_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include
ComputeNorm_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
ComputeNorm_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
