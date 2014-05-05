import libjad_DelphesAnalysis as j
import imp
types=imp.load_source('types','/home/hep/mc3909/jaf_development/scorpion/python/core/types.py')
from types import *
pmssm_219462978_cms7xsec=1.71952298774e-12
pmssm_219462978_cms7filelist = StringVector(["/home/hep/mc3909/jaf_development/scorpion/delphes-output.root"])
pmssm_219462978_pair = j.FilePair(pmssm_219462978_cms7xsec, pmssm_219462978_cms7filelist,0)
pmssm_219462978_pairmap = j.jad_FilePairMap()
pmssm_219462978_pairmap['CMS7']=pmssm_219462978_pair
pmssm_219462978_map = j.FileMap('pmssm_219462978', pmssm_219462978_pairmap)

masslist=[pmssm_219462978_map]
