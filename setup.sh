#!/bin/bash

#need to make sure we use compatible versions of:
# 1. ROOT (with mathmore & roofit compiled against python)
# 2. PYTHON
# 3. BOOST (compiled against python)


source /cvmfs/cms.cern.ch/cmsset_default.sh;
cmsrel CMSSW_8_0_4
cd CMSSW_8_0_4/src/; cmsenv; cd ../../


#change here to your ROOTSYS directory, doesn't compile for some reason with ROOTSYS in CMSSW
export ROOTSYS='/vols/build/cms/penning/root/root/'
export PYTHONPATH=$ROOTSYS/lib/:lib/
export SCORPIONSYS='/home/hep/bpenning/build/generator/scorpion/lib'

MY_LD_LIBRARY_PATH=$ROOTSYS/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/libjpg/8b/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/python/2.7.6/lib/python2.7:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.51.0/lib/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib64/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/xz/5.0.3__5.1.2alpha/lib:/home/hep/bpenning/build/scorpion/LandS/:



export LD_LIBRARY_PATH=$MY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH:$(pwd)/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/4.8.1/lib/
export PATH=$ROOTSYS/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_7_1_5/external/slc6_amd64_gcc530/bin/python:$PATH

echo "LD_LIBRARY_PATH is:" $LD_LIBRARY_PATH
echo "PYTHONPATH is:" $PYTHONPATH
echo "ROOTSYS is " $ROOTSYS
echo "PATH is " $PATH

