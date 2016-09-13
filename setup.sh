#!/bin/bash

#need to make sure we use compatible versions of:
# 1. ROOT (with mathmore & roofit compiled against python)
# 2. PYTHON
# 3. BOOST (compiled against python)

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc481
CMSRELEASE=CMSSW_7_1_5

if [ ! -d ${CMSRELEASE}/src ]; then 
    cmsrel ${CMSRELEASE}
fi
cd ${CMSRELEASE}/src; cmsenv; cd -

export SCORPIONDIR=`pwd`
export BOOSTDIR=/cvmfs/cms.cern.ch/${SCRAM_ARCH}/external/boost/1.51.0/include/

#install missing software  
if [ ! -d ${CMSRELEASE}/src/FWCore/Version ]; then 
    cd ${CMSRELEASE}/src
    git cms-addpkg FWCore/Version
    cd ${SCORPIONDIR} 
fi

#higgs combine tool 
if [ ! -d ${CMSRELEASE}/src/HiggsAnalysis/CombinedLimit ]; then 
    cd ${CMSRELEASE}/src
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git fetch origin
    git checkout v5.0.4   # try v5.0.1 if any issues occur
    scramv1 b clean; scramv1 b -j12 
    cd ${SCORPIONDIR} 
fi

#pheno limit 
if [ ! -d ${SCORPIONDIR}/pheno_limits ]; then 
    git clone git@github.com:bpenning/pheno_limits
fi
