import libjad_DelphesAnalysis as j
import imp
types=imp.load_source('types','/home/hep/mc3909/jaf_development/scorpion/python/core/types.py')
from types import *
DMtest_cms8xsec=1.71952298774e-12
DMtest_cms8filelist = StringVector(["/vols/cms04/mc3909/DM/1/delphes-output.root"])
DMtest_pair = j.FilePair(DMtest_cms8xsec, DMtest_cms8filelist)
DMtest_pairmap = j.jad_FilePairMap()
DMtest_pairmap['CMS8']=DMtest_pair
DMtest_map = j.FileMap('DMtest', DMtest_pairmap,1)
masslist=[DMtest_map]
