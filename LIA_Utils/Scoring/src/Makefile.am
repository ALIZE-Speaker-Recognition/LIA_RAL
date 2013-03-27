bin_PROGRAMS=Scoring
Scoring_SOURCES= \
	Scoring.cpp\
	ScoringMain.cpp

Scoring_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
Scoring_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
Scoring_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
