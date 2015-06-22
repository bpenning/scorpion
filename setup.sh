#!/bin/bash

#need to make sure we use compatible versions of:
# 1. ROOT (with mathmore & roofit compiled against python)
# 2. PYTHON
# 3. BOOST (compiled against python)

#export ROOTSYS='/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/'
export ROOTSYS='/vols/cms03/mc3909/root/'

MY_LD_LIBRARY_PATH=$ROOTSYS/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/libjpg/8b/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/python/2.7.6/lib/python2.7:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.51.0/lib/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/gcc/4.8.1/lib64/:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/xz/5.0.3__5.1.2alpha/lib:/vols/cms03/mc3909/LandS/:

export LD_LIBRARY_PATH=$MY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH:$(pwd)/lib:/home/hep/mc3909/scorpion/
export PATH=$ROOTSYS/bin:/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_1_5/external/slc6_amd64_gcc481/bin/python:$PATH
#export PYTHONPATH=$ROOTSYS/lib:/home/hep/jm1103/development/delphes-stuff-rewrite-boost-makefile-gen-optimised-xsec/lib
export PYTHONPATH=$ROOTSYS/lib:$(pwd)/lib:/home/hep/mc3909/scorpion/


echo "LD_LIBRARY_PATH is:" $LD_LIBRARY_PATH
echo "PYTHONPATH is:" $PYTHONPATH
echo "ROOTSYS is " $ROOTSYS
echo "PATH is " $PATH

#snippet from https://github.com/brynmathias/AnalysisV2/blob/dev/configure.py
#iclx_root_532={"name":"IC LX (64bit ROOT 5.32)",
#      "platform" : (False,"linux"),
#      "root_sys" : (True,"/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/"),
#      "root_sys_inc" : (False,"/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/include/"),
#      "root_sys_lib" : (False,"/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/lib/"),
#      "root_sys_bin" : (False,"/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/bin/"),
#      "python_inc" : (False,"/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/python/2.6.4/include/python2.6/"),
#      "python_lib" : (False,"-lpython2.6"),
#      "boost_python_lib" : (False,"-lboost_python"),
#      "root_extra_libs" : (False,"-lMathCore -lMathMore -lGenVector"),
#      "link_search_extra" : (False,"-L/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/python/2.6.4/lib -L/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/boost/1.47.0/lib"),
#      "incdir_extra" : (False,"-I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/python/2.6.4/include/python2.6 -I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/boost/1.47.0/include"),
#      "ld_path_extra" :
#(False,"/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/gcc/4.6.2/lib:/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/libjpg/8b/lib:/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/python/2.6.4/lib:/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/boost/1.47.0//lib:/vols/grid/ui/3.2.5-0/d-cache/dcap/lib64:/vols/grid/ui/3.2.5-0/d-cache/dcap/lib:/vols/grid/ui/3.2.5-0/glite/lib:/vols/grid/ui/3.2.5-0/glite/lib64:/vols/grid/ui/3.2.5-0/globus/lib:/vols/grid/ui/3.2.5-0/lcg/lib:/vols/grid/ui/3.2.5-0/lcg/lib64:/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/gcc/4.6.2/lib64:/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/xrootd/3.1.0-cms2/lib:/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/xz/5.0.3/lib"),
#      "python_env" : (False,"/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/python/2.6.4/bin"),
#      }

#CXX+=/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/gcc/4.6.2/bin/g++

#from MC-Project/DELPHES-ANALYSIS
#MY_LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hep/jm1103/Delphes_V_1.9/lib
#export LD_LIBRARY_PATH=$MY_LD_LIBRARY_PATH
#echo "LD_LIBRARY_PATH is:" $LD_LIBRARY_PATH
