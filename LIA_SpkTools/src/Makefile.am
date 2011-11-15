lib_LIBRARIES=libliatools.a

libliatools_a_SOURCES= \
	AccumulateStat.cpp\
	AccumulateJFAStat.cpp\
	FactorAnalysis.cpp\
	GeneralTools.cpp\
	IOFormat.cpp\
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

AM_CPPFLAGS=-I../include -I$(ALIZE_PATH)/include

all:
	mv libliatools.a $(top_srcdir)/lib/libliatools_$(OS)_$(ARCH)$(DEBUG).a

