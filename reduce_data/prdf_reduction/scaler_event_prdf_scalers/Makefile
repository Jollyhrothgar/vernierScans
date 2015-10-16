PACKAGE = clock

CXXFLAGS = -I.  -I$(ROOTSYS)/include -I$(OFFLINE_MAIN)/include -m32

LDFLAGS = -L$(OFFLINE_MAIN)/lib  -lpmonitor


HDRFILES = $(PACKAGE).h
LINKFILE = $(PACKAGE)LinkDef.h

ADDITIONAL_SOURCES = 


SO = lib$(PACKAGE).so

$(SO) : $(PACKAGE).cc $(PACKAGE)_dict.C $(ADDITIONAL_SOURCES) $(LINKFILE)
	$(CXX) $(CXXFLAGS) -o $@ -shared  $<  $(ADDITIONAL_SOURCES) $(PACKAGE)_dict.C $(LDFLAGS) 


$(PACKAGE)_dict.C : $(HDRFILES) $(LINKFILE)
	rootcint -f $@  -c $(CXXFLAGS) $^


.PHONY: clean

clean: 
	rm -f $(SO)  $(PACKAGE)_dict.C $(PACKAGE)_dict.h

