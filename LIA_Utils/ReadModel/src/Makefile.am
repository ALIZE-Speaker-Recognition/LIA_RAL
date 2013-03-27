bin_PROGRAMS=ReadModel
ReadModel_SOURCES= \
	ReadModel.cpp\
	ReadModelMain.cpp

ReadModel_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
ReadModel_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
ReadModel_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
