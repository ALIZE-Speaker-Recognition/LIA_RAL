bin_PROGRAMS=EigenChannel
EigenChannel_SOURCES= \
	EigenChannel.cpp\
	EigenChannelMain.cpp

EigenChannel_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(ALIZE_PATH)/include
EigenChannel_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib
EigenChannel_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin