AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror -std=c++0x 

BUILT_SOURCES = \
  testexternals.C

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libVernierTimeAnalysis.la

# this can be replaced with a shell command which invokes and captures root-config --libs
libVernierTimeAnalysis_la_LDFLAGS = \
  `root-config --evelibs` -lMinuit \
  -L$(MYINSTALL)/lib -lVernierAnalysis # Link to BeamPositionSteps library, be sure to include in src files when needed.

libVernierTimeAnalysis_la_SOURCES = \
  BeamWidthTime.C \
  BBCRate.cc \
  VernierTimeAnalysis_Dict.C

noinst_PROGRAMS = \
  testexternals 

testexternals_LDADD = \
  libVernierTimeAnalysis.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

#rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $(INCLUDES) $(DEFS) $(CINTFLAGS) $^
VernierTimeAnalysis_Dict.C: \
  BBCRate.h \
  BBCRateData.h \
  BeamWidthData.h \
  BeamWidthTime.h \
  TimeCompressedLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict.*
