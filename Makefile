PYTHON_INC=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/python/2.7.6/include/python2.7
ROOT_INC=$(ROOTSYS)/include
BOOST_INC=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.51.0/include/
#MINE
ROOFIT_INC=$(ROOFITSYS)/include/ #/vols/cms03/mc3909/root/roofit/roofitcore/inc
ROOFIT_INC2=$(ROOFITSYS)/include/ #/vols/cms03/mc3909/root/roofit/roofit/inc

LIMIT_INC=/vols/cms03/mc3909/LandS/include

PYTHON_LIB=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/python/2.7.6/lib/python2.7
ROOT_LIB=$(ROOTSYS)/lib
BOOST_LIB=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.51.0/lib/
#MINE
ROOFIT_LIB=$(ROOFITSYS)/lib/
LIMIT_LIB=/vols/cms03/mc3909/LandS
#DELPHES_LIB=/home/hep/mc3909/scorpion/
DELPHES_LIB=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

CXX=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/bin/g++

OBJDIR=obj
INCDIR=include
LIBDIR=lib

ANALYSESDIR=src/analyses
COREDIR=src/core

ANALYSES=$(wildcard $(ANALYSESDIR)/*.cc)
CORE=$(filter-out $(COREDIR)/python.cc $(COREDIR)/fileobject_class.cc, $(wildcard $(COREDIR)/*.cc))
#filter out python.cc and fileobject_class.cc which require specialist treatment

#the below would be used if you were to e.g. make an executable instead of a library.so etc
#ANALYSESTARGETS=$(subst .cc,,$(subst $(ANALYSESDIR)/,,$(ANALYSES)))
#CORETARGETS=$(subst .cc,,$(subst $(COREDIR)/,,$(CORE))) 

ANALYSESOBJ=$(subst .cc,.o,$(subst $(ANALYSESDIR)/,$(OBJDIR)/,$(ANALYSES)))
COREOBJ=$(subst .cc,.o,$(subst $(COREDIR)/,$(OBJDIR)/,$(CORE)))

CXXFLAGS=-c -fPIC -ansi -g -DLinux

all: $(ANALYSESOBJ) $(COREOBJ) $(OBJDIR)/fileobject_class.o $(OBJDIR)/delphesdictionary.o
	@echo Building Library
	$(CXX) -I$(INCDIR) -I$(PYTHON_INC) -I$(BOOST_INC) -I$(ROOT_INC) -I$(LIMIT_INC) -I$(ROOFIT_INC) \
		-L$(DELPHES_LIB) -L$(PWD) -L$(BOOST_LIB) -L$(PYTHON_LIB) -L$(ROOT_LIB) -L$(LIMIT_LIB) -L$(ROOFIT_LIB) \
	   	-g -fPIC $(COREDIR)/python.cc -shared \
		-lboost_python -lpython2.6 -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lDelphes -lRint -lPostscript \
		-lMatrix -lPhysics -lMathCore -lThread -lMathMore -lGenVector -lMinuit -lRooFit -lRooFitCore -pthread -lm -ldl -llimitcode \
	   	-rdynamic -pthread -m64 -o $(LIBDIR)/libjad_DelphesAnalysis.so $(wildcard $(OBJDIR)/*.o)
	@echo --DONE--

$(ANALYSESOBJ): $(ANALYSES)
	@echo Compiling $(subst .o,.cc,$(subst $(OBJDIR),$(ANALYSESDIR),$@))
	$(CXX) $(CXXFLAGS) $(subst .o,.cc,$(subst $(OBJDIR),$(ANALYSESDIR),$@)) -o $@ -pthread -m64 -I$(INCDIR) -I$(ROOT_INC)

$(COREOBJ): $(CORE)
	@echo Compiling $(subst .o,.cc,$(subst $(OBJDIR),$(COREDIR),$@))
	$(CXX) $(CXXFLAGS) $(subst .o,.cc,$(subst $(OBJDIR),$(COREDIR),$@)) -o $@ -pthread -m64 -I$(INCDIR) -I$(ROOT_INC) -I$(LIMIT_INC) -I$(ROOFIT_INC) -I$(ROOFIT_INC2)

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
