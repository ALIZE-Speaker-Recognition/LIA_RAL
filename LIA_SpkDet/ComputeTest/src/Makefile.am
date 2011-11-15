bin_PROGRAMS=ComputeTest
ComputeTest_SOURCES= \
	ComputeTest.cpp\
	ComputeTestMain.cpp

ComputeTest_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include
ComputeTest_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
ComputeTest_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin