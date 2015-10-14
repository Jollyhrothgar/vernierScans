AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -Wall -Werror -std=c++0x 

BUILT_SOURCES = \
	testexternals.C

AM_CPPFLAGS = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include

lib_LTLIBRARIES = \
	libVernierHourglass.la

# this can be replaced with a shell command which invokes and captures root-config --libs
libVernierHourglass_la_LDFLAGS = \
	`root-config --evelibs` -lMinuit	\
	-L$(MYINSTALL) -lVernierAnalysis

libVernierHourglass_la_SOURCES = \
	HourglassConfiguration.cc \
	HourglassData.cc \
	HourglassSimulation.cc \
	VernierHourglass_Dict.C

noinst_PROGRAMS = \
	testexternals 

testexternals_LDADD = \
	libVernierHourglass.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "	return 0;" >> $@
	echo "}" >> $@

VernierHourglass_Dict.C: \
	HourglassConfiguration.h \
	HourglassData.h \
	HourglassSimulation.h \
	HourglassLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict.*
