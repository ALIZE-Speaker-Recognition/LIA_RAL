bin_PROGRAMS=ExtractParams
ExtractParams_SOURCES= \
	ExtractParams.cpp\
	ExtractParamsMain.cpp

ExtractParams_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
ExtractParams_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
ExtractParams_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
