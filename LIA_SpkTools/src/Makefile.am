lib_LIBRARIES=libliatools.a

libliatools_a_SOURCES= \
	AccumulateStat.cpp\
	AccumulateJFAStat.cpp\
	AccumulateTVStat.cpp\
	FactorAnalysis.cpp\
	GeneralTools.cpp\
	IOFormat.cpp\
	PldaTools.cpp\
	SegTools.cpp\
	TopGauss.cpp\
	UnsupervisedTools.cpp\
	ClusteringCriterion.cpp\
	FileInfo.cpp\
	Hmm.cpp\
	ScoreWarp.cpp\
	SuperVectors.cpp\
	TrainTools.cpp\
	Tools.cpp\
	liatools.cpp

AM_CPPFLAGS=-I../include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen -I$(LAPACK_PATH)/lapacke/include/

all:
	mv libliatools.a $(top_srcdir)/lib/libliatools_$(OS)_$(ARCH)$(DEBUG).a

