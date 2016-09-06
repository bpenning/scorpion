LIB  := lands.so
EXEC :=	test/CLs.exe \
	test/Bayesian.exe \
	test/lands.exe \
        test/CLs_dataCard.exe \

SOURCES := $(wildcard src/*.cc)
OBJECTS := $(patsubst src/%.cc, bin/%.o, $(SOURCES))

ifdef ROOFITSYS
    RF_CFLAGS      := -I/vols/cms03/mc3909/root/roofit/roofitcore/inc -I/vols/cms03/mc3909/root/roofit/roofit/inc
    RF_LINKERFLAGS := -L/vols/cms03/mc3909/root/lib
endif

#CXX      := g++
CXX      := /vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/gcc/4.6.2/bin/g++
CXXFLAGS := -O2 -g -fPIC -Iinclude $(shell root-config --cflags) -I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include ${RF_CFLAGS}

#LD       := g++
LD       := /vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/gcc/4.6.2/bin/g++

LDFLAGS  := $(shell root-config --libs --ldflags) -L/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/lib -lMathMore -lMinuit -lRooFit -lRooFitCore -lFoam ${RF_LINKERFLAGS}

ifeq ($(shell root-config --platform),macosx)
	MACOSXFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -Xlinker -bind_at_load -flat_namespace
endif

# Enable CPU profiling of... requires google performance tools.
# You also need to set CPUPROFILE enviroanment variable.
#PROFILE=1
ifdef PROFILE
  CXXFLAGS += -DLANDS_PROFILE
  LDFLAGS  += -lprofiler
endif


 # stop removing intermediate files
.SECONDARY:

all: $(LIB) $(EXEC)

bin/%.o: src/%.cc
	$(CXX) $(CXXFLAGS)  $< -c -o $@

${LIB}: ${OBJECTS}
	$(LD) $(LDFLAGS) $(MACOSXFLAGS) -shared  ${OBJECTS} -o $@

%.exe: %.cc ${OBJECTS}
	$(LD) $(CXXFLAGS) $(LDFLAGS) -o $@ $< ${OBJECTS}

#for MSSMA
drawMSSMA: test/drawMSSMA.exe

clean:
	@rm -v -f \
	bin/*.o \
	test/*.exe \
	*.so

cleanall: clean
	@rm -v -f \
	test/log* \
	test/*.gif \
	test/*.eps \
	test/*.root \
	*/*~ *~
