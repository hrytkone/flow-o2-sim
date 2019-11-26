PROGRAM       = flowO2Sim

version       = toka
CXX           = g++
CXXFLAGS      = -Ofast -Wall -g -D$(version)

LD            = g++
LDFLAGS       = -Ofast
SOFLAGS       = -shared
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)

HDRSDICT =  src/JHistos.h \
			src/JInputs.h \

ALICEPATH = $(BASEDIR)
O2INCLUDES = -I$(ALICEPATH)/O2/latest/include -I$(ALICEPATH)/O2/latest/include/AliTPCCommon -I$(ALICEPATH)/O2/latest/include/DataFormatsTOF -I$(ALICEPATH)/FairRoot/latest/include -I$(ALICEPATH)/FairLogger/latest/include -I$(ALICEPATH)/boost/latest/include -I$(ALICEPATH)/ms_gsl/latest/include -I$(ALICEPATH)/O2/latest/include/CommonConstants -I$(ALICEPATH)/O2/latest/include/GPU

HDRS    += $(HDRSDICT) nanoDict.h

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM) $(PROGRAM2)

$(PROGRAM):     $(OBJS) $(PROGRAM).C
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -lEG -lPhysics -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(FFTWINC) $(O2INCLUDES) -L$(ALICEPATH)/O2/latest/lib -o $(PROGRAM)
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS)  -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $(FFTWINC) $<


clean:
	rm -f $(PROGRAM) *.o nanoDict*

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -c -D$(version) $(HDRSDICT)
