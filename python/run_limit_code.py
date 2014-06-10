import os
import ROOT as r
import libjad_DelphesAnalysis as j
from python.core.types import *
from samples.exp_values import *

def filemap_from_dict(filemap_dict):
    #exctract relevant parameters form dicionary
    xsec=filemap_dict['xsec']
    rootfiles=filemap_dict['rootfiles']
    experiment=filemap_dict['experiment']
    internal_name=filemap_dict['internal_name']
    #generate and return the file map
    filelist = StringVector(rootfiles)
    pair = j.FilePair(xsec, filelist)
    pairmap = j.jad_FilePairMap()
    pairmap[experiment]=pair
    return j.FileMap(internal_name, pairmap,1)

def runlim(outdir,filemap_dict,com,ss5b=0,os5b=0,lp5b=0,alphat7b=0,alphat7bb=0,
        alphat12b=0,ss11b=0,alphat40b=0,ss40b=0,lp20b=0,ge3lp20b=0,mt220b=0):
    print('XSEC: {0} barn'.format(filemap_dict['xsec']))
    #jaf needs the directory to end in '/'
    if not outdir[:-1] =='/':
        outdir+='/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    calculateRyes = True
    calculateRno = False
    
    calculateRcomboyes = True
    calculateRcombono = False
    
    geninfoyes = True
    geninfono = False
    
    docomboyes = True
    docombono = False
    
    writestatfileyes = True
    writestatfileno = False
    
    # 7 TeV analyses
    alphat7_5b = j.AlphaTb('alphaTb7_analysis5','CMS7', 32, 4.98, bg_at7b, bgunc_at7b, data_at7b, 'combined', calculateRno)
#    atlas5 = j.ATLAS5('ATLAS5 analysis','ATLAS', 11, 4.710, bg_atlas5, bgunc_atlas5, data_atlas5, 'individual', calculateRno)
#    lp5 = j.LP('LP_analysis5','CMS7', 18, 4.7, bg_lp5, bgunc_lp5, data_lp5, 'combined', calculateRno)
#    monojet5 = j.MonoJet('MonoJet_analysis75','CMS', 4, 5.0, bg_monojet5, bgunc_monojet5, data_monojet5, 'individual', calculateRno)
#    os5 = j.OS('OS_analysis5','CMS7', 6, 4.98, bg_os5, bgunc_os5, data_os5, 'combined', calculateRno)
#    ss5 = j.SS('SS_analysis5','CMS7', 3, 4.98, bg_ss5, bgunc_ss5, data_ss5, 'combined', calculateRno)
    
    # 8 TeV analyses
#    alphat8_12b = j.AlphaTb8('alphaTb8_analysis12','CMS8', 59, 11.7, bg_at8b, bgunc_at8b, data_at8b, 'combined', calculateRno)
#    ss8_11b = j.SSb8('SSb8_analysis11','CMS8', 1, 10.5, bg_ssb8, bgunc_ssb8, data_ssb8, 'combined', calculateRno)
#    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8',16,19.5,bg_lp8_20b,bgunc_lp8_20b,data_lp8_20b,'combined',calculateRno)
#    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8',16,19.5,bg_lp8_20b,bgunc_lp8_20b,data_lp8_20b,'individual',calculateRno)
    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',123,19.5,bg_zerolepmt2_8_20,
            bgunc_zerolepmt2_8_20,data_zerolepmt2_8_20,'combined',calculateRno)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',52,19.5,bg_zerolepmt2_8_20_stop_test,
#            bgunc_zerolepmt2_8_20_stop_test,data_zerolepmt2_8_20_stop_test,'combined',calculateRno)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',25,19.5,bg_zerolepmt2_8_20_gluino_test,
#            bgunc_zerolepmt2_8_20_gluino_test,data_zerolepmt2_8_20_gluino_test,'combined',calculateRno)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',25,19.5,bg_zerolepmt2_8_20_gluino_test,
#            bgunc_zerolepmt2_8_20_gluino_test,data_zerolepmt2_8_20_gluino_test,'individual',calculateRno)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',123,19.5,bg_zerolepmt2_8_20,
#            bgunc_zerolepmt2_8_20,data_zerolepmt2_8_20,'individual',calculateRno)
    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8',8,19.5,bg_lp8_20b,
            bgunc_lp8_20b,data_lp8_20b,'individual',calculateRno)
#    data_cms3l8=IntVector([int(round(x)) for x in bg_cms3l8])
    ge3lp8_20b = j.Cms3Lepton20Fb('GE3LP_analysis20','CMS8',192,19.5,bg_cms3l8,
            bgunc_cms3l8,data_cms3l8,'combined',calculateRno)

    
    # for extrapolations:
#    bg_at8b40 = DoubleVector([x * 3.4 for x in bg_at8b])
#    data_at8b40 = IntVector([int(x * 3.4) for x in data_at8b])
    
#    bg_ssb40 = DoubleVector([x * 3.8 for x in bg_ssb8])
#    data_ssb40 = IntVector([int(x * 3.8) for x in data_ssb8])
    
    # weather forecast (40/fb):
#    alphat8_40b = j.AlphaTb8('alphaTb8_analysis40','CMS8', 59, 40.0, bg_at8b40, bgunc_at8b, data_at8b40, 'combined', calculateRno)
#    ss8_40b = j.SSb8('SSb8_analysis40','CMS8', 1, 10.5, bg_ssb40, bgunc_ssb8, data_ssb40, 'combined', calculateRno)

    # 20/fb 8 TeV searches
#    ss8HighPt = j.SS8high('SS8high_analysis20','CMS8', 24, 19.5, bg_ss8HighPt, bgunc_ss8HighPt, data_ss8HighPt, 'combined', calculateRno)
#    ss8LowPt = j.SS8low('SS8low_analysis20','CMS8', 24, 19.5, bg_ss8LowPt, bgunc_ss8LowPt, data_ss8LowPt, 'combined', calculateRno)
#    zerolep8 = j.ZeroLep8('SSb8_analysis40','CMS8', 48, 19.4, bg_zerolep8, bgunc_zerolep8, data_zerolep8, 'combined', calculateRno)
    
    mgr = j.AnalysisManager(outdir, geninfono) #bool for geninfo
    
    if com == 7:
#        if ss5b:
#              mgr.add(ss5)
#        if os5b:
#             mgr.add(os5)
#        if lp5b:
#             mgr.add(lp5)
        if alphat7b:
             mgr.add(alphat7_5)
#        if alphat7bb:
#             mgr.add(alphat7_5b)
    elif com == 8:
#        if alphat12b:
#             mgr.add(alphat8_12b)
#        if ss11b:
#             mgr.add(ss8_11b)
#        if alphat40b:
#             mgr.add(alphat8_40b)
#        if ss40b:
#             mgr.add(ss8_40b)
        if lp20b:
            mgr.add(lp8_20b)
        if ge3lp20b:
            mgr.add(ge3lp8_20b)
        if mt220b:
            mgr.add(mt2_20b)
    else:
         print "WRONG ENERGY"
         sys.exit(1) 

    filemap=filemap_from_dict(filemap_dict)
    mgr.run(FileMapVector([filemap]), "")
    mgr.limit(0.20, writestatfileyes, docomboyes, calculateRcombono)
    mgr.write()
  
