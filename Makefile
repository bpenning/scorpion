CXX=g++

CXXFLAGS=-c -fPIC -ansi -g -DLinux

OBJDIR=$(SCORPIONDIR)/obj
INCDIR=$(SCORPIONDIR)/include
LIBDIR=$(SCORPIONDIR)/lib

PYTHON_INC=$(shell python-config --includes) 
ROOT_INC=$(shell root-config --incdir)
BOOST_INC=$(BOOSTDIR)/include
ROOFIT_INC=$(ROOFITSYS)/include
LIMIT_INC=$(SCORPIONDIR)/LandS/include

BOOST_LIB=$(BOOSTDIR)/lib
LIMIT_LIB=$(SCORPIONDIR)/LandS
DELPHES_LIB=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))


ANALYSESDIR=$(SCORPIONDIR)/src/analyses
ANALYSES=$(wildcard $(ANALYSESDIR)/*.cc)
COREDIR=$(SCORPIONDIR)/src/core
CORE=$(filter-out $(COREDIR)/python.cc $(COREDIR)/fileobject_class.cc, $(wildcard $(COREDIR)/*.cc))
#filter out python.cc and fileobject_class.cc which require specialist treatment

#the below would be used if you were to e.g. make an executable instead of a library.so etc
#ANALYSESTARGETS=$(subst .cc,,$(subst $(ANALYSESDIR)/,,$(ANALYSES)))
#CORETARGETS=$(subst .cc,,$(subst $(COREDIR)/,,$(CORE))) 

ANALYSESOBJ=$(subst .cc,.o,$(subst $(ANALYSESDIR)/,$(OBJDIR)/,$(ANALYSES)))
COREOBJ=$(subst .cc,.o,$(subst $(COREDIR)/,$(OBJDIR)/,$(CORE)))



all: $(ANALYSESOBJ) $(COREOBJ) $(OBJDIR)/fileobject_class.o $(OBJDIR)/delphesdictionary.o
	@echo Building Library
	$(CXX) -I$(INCDIR) $(PYTHON_INC) -I$(BOOST_INC) -I$(ROOT_INC) -I$(LIMIT_INC) -I$(ROOFIT_INC) \
	   	-g -fPIC $(COREDIR)/python.cc -shared \
		$(ROOTLIB) -L$(shell root-config --libdir) -lMathMore -lGenVector -lMinuit \
		$(PYTHONLIB) -L$(BOOST_LIB) -lboost_python -L$(DELPHES_LIB) -lDelphes \
		-L$(ROOFITSYS)/lib -lRooFit -lRooFitCore \
		-L$(LIMIT_LIB) -llimitcode \
	   	-m64 -o $(LIBDIR)/libjad_DelphesAnalysis.so $(wildcard $(OBJDIR)/*.o)
	@echo --DONE--


$(ANALYSESOBJ): $(ANALYSES)
	@echo Compiling $(subst .o,.cc,$(subst $(OBJDIR),$(ANALYSESDIR),$@))
	$(CXX) $(CXXFLAGS) $(subst .o,.cc,$(subst $(OBJDIR),$(ANALYSESDIR),$@)) -o $@ -pthread -m64 -I$(INCDIR) -I$(ROOT_INC)

$(COREOBJ): $(CORE)
	@echo Compiling $(subst .o,.cc,$(subst $(OBJDIR),$(COREDIR),$@))
	$(CXX) $(CXXFLAGS) $(subst .o,.cc,$(subst $(OBJDIR),$(COREDIR),$@)) -o $@ -pthread -m64 -I$(INCDIR) -I$(ROOT_INC) -I$(LIMIT_INC) -I$(ROOFIT_INC) 

$(OBJDIR)/fileobject_class.o: $(COREDIR)/fileobject_class.cc
	@echo Compiling $(COREDIR)/fileobject_class.cc
	$(CXX) $(CXXFLAGS) $(COREDIR)/fileobject_class.cc -o $@ -I$(INCDIR)

$(OBJDIR)/delphesdictionary.o: $(INCDIR)/jBlockClasses.hh $(INCDIR)/jBlockClassesLinkDef.h
	@echo Compiling DELPHES dictionary
	rootcint -f DictOutput.cxx -c $(INCDIR)/jBlockClasses.hh $(INCDIR)/jBlockClassesLinkDef.h
	$(CXX) $(CXXFLAGS) DictOutput.cxx -o $(OBJDIR)/delphesdictionary.o -pthread -m64 -I$(ROOT_INC)
	rm DictOutput.cxx DictOutput.h

clean:
	@rm -f $(wildcard obj/*.o)
	@rm -f $(wildcard lib/*.so)
