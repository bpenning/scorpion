import libjad_DelphesAnalysis as j
from python.core.types import *

pmssm_219462978_cms7xsec=1.71952298774e-12
pmssm_219462978_cms7filelist = StringVector(["/vols/cms04/mc3909/FTpoint/D3output.root"])
pmssm_219462978_pair = j.FilePair(pmssm_219462978_cms7xsec, pmssm_219462978_cms7filelist)
pmssm_219462978_pairmap = j.jad_FilePairMap()
pmssm_219462978_pairmap['CMS7']=pmssm_219462978_pair
pmssm_219462978_map = j.FileMap('pmssm_219462978', pmssm_219462978_pairmap,1)

masslist=[pmssm_219462978_map]
