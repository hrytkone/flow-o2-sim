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

ALICEPATH = /home/heimarry/alice/sw/slc7_x86-64
O2INCLUDES = -I$(ALICEPATH)/O2/1.0.0-1/include -I$(ALICEPATH)/O2/1.0.0-1/include/AliTPCCommon -I$(ALICEPATH)/O2/1.0.0-1/include/DataFormatsTOF -I$(ALICEPATH)/FairRoot/c672f280ec-54/include -I$(ALICEPATH)/FairLogger/v1.2.0-9/include -I$(ALICEPATH)/boost/v1.68.0-30/include -I$(ALICEPATH)/ms_gsl/1-8/include -I$(ALICEPATH)/O2/1.0.0-1/include/CommonConstants -I$(ALICEPATH)/O2/1.0.0-1/include/GPU

HDRS    += $(HDRSDICT) nanoDict.h

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM) $(PROGRAM2)

$(PROGRAM):     $(OBJS) $(PROGRAM).C
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -lEG -lPhysics -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(FFTWINC) $(O2INCLUDES) -L/home/heimarry/alice/sw/slc7_x86-64/O2/1.0.0-1/lib -o $(PROGRAM)
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
