ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
UNAME = $(shell uname -s)
ifeq ($(UNAME), Darwin)
DllSuf        = dylib
endif
ifeq ($(UNAME), Linux)
DllSuf        = so
endif

OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTLIBS_NoTMVA = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)
ifeq ($(UNAME), Linux)
ROOTLIBS       += -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA
ROOTLIBS_NoTMVA+= -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam
ROOTGLIBS      += -lMinuit -lMathMore -lMinuit2 -lRooFitCore -lRooFit -lRooStats -lFoam -lTMVA
endif
ROOTLIBS       += -L$(ROOFITSYS)/lib

# Linux with egcs
DEFINES         = -DNO_ORCA_CLASSES -I..
CXX             = g++
CXXFLAGS        = -O -Wall -fPIC $(DEFINES)  #-I./TMVA/include
ifeq ($(UNAME), Darwin)
CXXFLAGS       += -I/opt/local/include
endif
LD              = g++
LDFLAGS         = -g -O -Wall -fPIC
ifeq ($(UNAME), Darwin)
SOFLAGS         = -dynamiclib -undefined dynamic_lookup # use only -dynamiclib for TTAB
endif
ifeq ($(UNAME), Linux)
SOFLAGS         = -shared
endif

CXXFLAGS       += $(ROOTCFLAGS)
LIBS            = -I./TMVA/include -L./TMVA/lib $(ROOTLIBS) -lEG -I.. -L. -L../TopTreeProducer/src -L../TopTreeAnalysisBase
LIBS_NoTMVA     = $(ROOTLIBS_NoTMVA) -lEG -I.. -L. -L../TopTreeProducer/src -L../TopTreeAnalysisBase
ifeq ($(UNAME), Darwin)
LIBS           += -I/opt/local/include
LIBS_NoTMVA    += -I/opt/local/include
endif
GLIBS           = $(ROOTGLIBS)
#------------------------------------------------------------------------------
#SOURCES         = $(wildcard ../TopTreeAnalysisBase/Tools/src/*.cc ../TopTreeAnalysisBase/KinFitter/src/*.cc Tools/src/*.cc)
SOURCES         = $(wildcard Tools/src/*.cc)
#HEADERS         = $(wildcard ../TopTreeAnalysisBase/Tools/interface/*.h ../TopTreeAnalysisBase/Kinfitter/interface/*.h Tools/interface/*.h)
HEADERS         = $(wildcard Tools/interface/*.h)
OBJECTS	        = $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS	        = $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS        = $(SOURCES:.$(SrcSuf)=.$(DllSuf))
#for libTopTreeAnaContent.so
#SOURCESDIC      = $(wildcard ../TopTreeAnalysisBase/Reconstruction/src/Observables.cc ../TopTreeAnalysisBase/Reconstruction/src/MEzCalculator.cc ../TopTreeAnalysisBase/Content/src/*.cc ../TopTreeProducer/src/TRoot*.cc Tools/src/*.cc)
SOURCESDIC      = $(wildcard Tools/src/*.cc)
HEADERSDIC      = $(wildcard Tools/interface/*.h)
OBJECTSDIC      = $(SOURCESDIC:.$(SrcSuf)=.$(ObjSuf))


all: libTopWidthAna80.$(DllSuf)
	cp libTopWidthAna80.$(DllSuf) ~/lib/

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(OBJECTSDIC) $(DEPENDS) macros/*.exe *Dict.* *.$(DllSuf) core 

.SUFFIXES: .$(SrcSuf) .C .o .$(DllSuf)

###

Dict.$(SrcSuf): $(HEADERSDIC) ./LinkDef.h
	@echo "Generating dictionary Dict..."
	@rootcint -f Dict.$(SrcSuf) -c $(DEFINES) $(HEADERSDIC) ./LinkDef.h

libTopWidthAna80.$(DllSuf): $(OBJECTS) Dict.o
	@echo "Building libTopWidthAna..."
	$(LD) $(LIBS_NoTMVA) $(SOFLAGS) $(LDFLAGS) $+ -o $@

ADDLIBS_MACROS = -lMLP -lTreePlayer -lXMLIO

macros/%.exe: macros/%.cc $(HEADERS) libTopWidthAna80.$(DllSuf)
	$(LD) -lTopWidthAna80 $(LIBS_NoTMVA) $(ADDLIBS_MACROS) -I`root-config --incdir` $< $(LDFLAGS) -o $@


SOURCES_MACROS = $(wildcard macros/*.cc)

macros: $(SOURCES_MACROS:.cc=.exe)

