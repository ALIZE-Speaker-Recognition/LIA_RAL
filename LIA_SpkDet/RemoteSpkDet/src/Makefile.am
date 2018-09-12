bin_PROGRAMS=SpkDetServer
SpkDetServer_SOURCES= \
	../../SimpleSpkDetSystem/src/SimpleSpkDetSystem.cpp\
	../../NormFeat/src/NormFeat.cpp\
	../../EnergyDetector/src/EnergyDetector.cpp\
	SpkDetServer.cpp\
	SpkDetServerMain.cpp

SpkDetServer_CPPFLAGS=-I../include -I$(top_srcdir)/LIA_SpkTools/include -I$(top_srcdir)/LIA_SpkDet/NormFeat/include -I$(top_srcdir)/LIA_SpkDet/EnergyDetector/include -I$(top_srcdir)/LIA_SpkDet/SimpleSpkDetSystem/include -I$(ALIZE_PATH)/include -I$(top_srcdir)/include/Eigen $(SPRO_INCLUDE_PATH)
SpkDetServer_LDFLAGS=-L$(top_srcdir)/lib -L$(ALIZE_PATH)/lib -L$(top_srcdir)/include/Pthread $(SPRO_LIB_PATH)
SpkDetServer_LDADD=-lliatools_${OS}_${ARCH}${DEBUG} $(lib_alize) $(LIBTHREAD) $(LIBSPRO)

all:
	mv ${bin_PROGRAMS} ${top_srcdir}/bin
