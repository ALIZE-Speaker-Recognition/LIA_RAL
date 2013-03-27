bin_PROGRAMS=Svm
Svm_SOURCES= \
	libsvm.cpp\
	Svm.cpp\
	SvmMain.cpp

Svm_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen
Svm_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
Svm_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
