AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror -std=c++0x 

BUILT_SOURCES = \
  testexternals.C

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libPRDFAnalysis.la

# this can be replaced with a shell command which invokes and captures root-config --libs
# add -l RooFitCore if using RooFit is desireable
libPRDFAnalysis_la_LDFLAGS = \
  -L$(ROOTSYS)/lib `root-config --evelibs` -lMinuit -lFoam\
  -L$(MYINSTALL)/lib -lVernierAnalysis # for BeamPositionSteps

libPRDFAnalysis_la_SOURCES = \
  BeamWidthScalers.C \
  PRDFBBCRate.cc \
  PRDFAnalysis_Dict.C

noinst_PROGRAMS = \
  testexternals 

testexternals_LDADD = \
  libPRDFAnalysis.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

#rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $(INCLUDES) $(DEFS) $(CINTFLAGS) $^
PRDFAnalysis_Dict.C: \
  BeamWidthScalers.h \
  PRDFBBCRate.h \
  PRDFPacketLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict.*
