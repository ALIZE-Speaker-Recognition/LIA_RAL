bin_PROGRAMS=TrainWorld
TrainWorld_SOURCES= \
	TrainWorld.cpp\
	TrainWorldMain.cpp

TrainWorld_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
TrainWorld_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
TrainWorld_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin