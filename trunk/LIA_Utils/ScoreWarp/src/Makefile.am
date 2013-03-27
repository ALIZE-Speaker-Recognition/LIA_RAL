bin_PROGRAMS=ScoreWarp
ScoreWarp_SOURCES= \
	ScoreWarpMain.cpp

ScoreWarp_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
ScoreWarp_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
ScoreWarp_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
