#!/bin/bash

#need to make sure we use compatible versions of:
# 1. ROOT (with mathmore & roofit compiled against python)
# 2. PYTHON
# 3. BOOST (compiled against python)


#source /vols/cms/grid/setup.sh
source /home/hep/bpenning/bin/setup_cmssw.sh
source /vols/grid/cms/setup.sh
export ROOTSYS='/vols/build/cms/penning/root/root/'
export PYTHONPATH=$ROOTSYS/lib/:lib/
export SCORPIONSYS='/home/hep/bpenning/build/generator/scorpion/lib'

MY_LD_LIBRARY_PATH=$ROOTSYS/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/libjpg/8b/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/python/2.7.6/lib/python2.7:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.51.0/lib/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib64/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/xz/5.0.3__5.1.2alpha/lib:/home/hep/bpenning/build/scorpion/LandS/:

#MY_LD_LIBRARY_PATH=$ROOTSYS/lib:

#MY_LD_LIBRARY_PATH=$ROOTSYS/lib:/vols/build/cms/penning/CMSSW_8_0_4/biglib/slc6_amd64_gcc530:/vols/build/cms/penning/CMSSW_8_0_4/lib/slc6_amd64_gcc530:/vols/build/cms/penning/CMSSW_8_0_4/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_4/biglib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_4/lib/slc6_amd64_gcc530:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_4/external/slc6_amd64_gcc530/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/llvm/3.7.1/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib



export LD_LIBRARY_PATH=$MY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH:$(pwd)/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/4.8.1/lib/

#export LD_LIBRARY_PATH=$MY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH:$(pwd)/lib:$(ROOTSYS)/lib
export PATH=$ROOTSYS/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_7_1_5/external/slc6_amd64_gcc530/bin/python:$PATH

echo "LD_LIBRARY_PATH is:" $LD_LIBRARY_PATH
echo "PYTHONPATH is:" $PYTHONPATH
echo "ROOTSYS is " $ROOTSYS
echo "PATH is " $PATH

