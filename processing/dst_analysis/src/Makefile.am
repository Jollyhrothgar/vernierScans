AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror -std=c++0x 

BUILT_SOURCES = \
  testexternals.C

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libVernierAnalysis.la

noinst_HEADERS = \
  VernierConstants.h 

# this can be replaced with a shell command which invokes and captures root-config --libs
libVernierAnalysis_la_LDFLAGS = \
  `root-config --evelibs` -lMinuit \
  -L$(OFFLINE_MAIN)/lib 

libVernierAnalysis_la_SOURCES = \
  BeamPositionSteps.C \
  BBCEfficiency.C \
  BBCRateSteps.C \
  BeamWidth.C \
  WcmDcctManager.C \
  VernierTreeVariables.C \
  Vernier_Dict.C

noinst_PROGRAMS = \
  testexternals 

testexternals_LDADD = \
  libVernierAnalysis.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

#rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $(INCLUDES) $(DEFS) $(CINTFLAGS) $^
Vernier_Dict.C: \
 	 BeamPositionSteps.h \
 	 BBCEfficiency.h \
 	 BBCRateSteps.h \
  	BeamWidth.h \
  	BeamSeparationData.h \
  	WcmDcctManager.h \
	VernierTreeVariables.h \
  	VernierLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict.*
