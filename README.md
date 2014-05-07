scorpion
========

Detector simulation + searches + CLs calculation


requirements
============

* gsl
* boost
* root
* Delphes

installation
============

Compile root as follows

    ./configure  --enable-roofit --enable-python --with-python-incdir=<path> --with-python-libdir=<path> --enable-mathmore --with-gsl-incdir=<path> --with-gsl-libdir=<path>

To install:

    source setup-2.sh
    make
