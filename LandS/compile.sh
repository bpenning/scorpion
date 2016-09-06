#!/bin/bash

/vols/sl5_exp_software/cms/slc5_amd64_gcc462/external/gcc/4.6.2/bin/g++ -O2 -g -fPIC -Iinclude -pthread -m64 -I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/include -I/vols/sl5_exp_software/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include   src/BayesianBase.cc -c -o bin/BayesianBase.o

