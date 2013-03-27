bin_PROGRAMS=Hist
Hist_SOURCES= \
	Hist.cpp\
	HistMain.cpp

Hist_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
Hist_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
Hist_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
