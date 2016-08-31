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

Currently on slc6 with setup of CMSSW_8_0_4 (included in setup script). Note: conflics with CMSSW root version, please edid setup.sh to point to your manual root installation

To install:
```bash
    source setup.sh
    make
```
To see options run with:
```bash
./run_jaf.py --help
```
Example for runnning DM searches set up for CMS (8TeV)
```bash
   source run_example.sh
```


