bin_PROGRAMS=ShiftedDeltaFeat
ShiftedDeltaFeat_SOURCES= \
	ShiftedDeltaFeat.cpp\
	ShiftedDeltaFeatMain.cpp

ShiftedDeltaFeat_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
ShiftedDeltaFeat_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib -L$(top_srcdir)/include/Pthread
ShiftedDeltaFeat_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin