bin_PROGRAMS=EigenVoice
EigenVoice_SOURCES= \
	EigenVoice.cpp\
	EigenVoiceMain.cpp

EigenVoice_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include
EigenVoice_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
EigenVoice_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
