bin_PROGRAMS=FusionScore
FusionScore_SOURCES= \
	FusionScore.cpp\
	FusionScoreMain.cpp

FusionScore_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
FusionScore_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
FusionScore_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
