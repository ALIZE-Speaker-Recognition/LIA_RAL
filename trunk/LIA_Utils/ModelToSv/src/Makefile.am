bin_PROGRAMS=ModelToSv
ModelToSv_SOURCES= \
	ModelToSvMain.cpp

ModelToSv_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
ModelToSv_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
ModelToSv_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
