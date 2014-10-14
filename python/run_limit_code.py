import os
import ROOT as r
import libjad_DelphesAnalysis as j
from python.core.types import *
from samples.exp_values import *

def filemap_from_dict(filemap_dict, internal_name):
    #exctract relevant parameters form dicionary
    pairmap = j.jad_FilePairMap()
    for experiment, filepair_dict in filemap_dict.items():
        xsec = filepair_dict['xsec']
        rootfiles = filepair_dict['rootfiles']
        #generate and return the file map
        filelist = StringVector(rootfiles)
        pair = j.FilePair(xsec, filelist)
        pairmap[experiment] = pair
    return j.FileMap(internal_name, pairmap,1)

def runlim(outdir, filemap_dict, gen_info, use_event_weights, calculate_r=False, 
        calculate_r_combo=False,  do_combo=True, write_stats_file=True, 
        alphat7bb=0, monojet20b=0, mt220b=0, alphat20b=0, lp20b=0,os5b=0, 
        ss820b=0, ge3lp20b=0):

    #jaf needs the directory to end in '/'
    if not outdir[:-1] =='/':
        outdir+='/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # 7 TeV analyses
    alphat7_5b = j.AlphaTb('alphaTb7_analysis5','CMS7', 32, 4.98, bg_at7b, 
            bgunc_at7b, data_at7b, 'combined', calculate_r)

    # Final 8 TeV searches
    monojet20 = j.MonoJet8('MonoJet20_analysis20', 'CMS8', 7, 19.5, bg_monojet20,
            bgunc_monojet20, data_monojet20, 'strongest', calculate_r)
    alphat8_20b = j.AlphaT20b('alphaT20b_analysis20', 'CMS8',75, 18.5, 
            bg_alphat20b, bgunc_alphat20b, data_alphat20b, 'combined', 
            calculate_r)
    mt2_20b = j.ZeroLepMt2('MT2_analysis20', 'CMS8', 123, 19.5, 
            bg_zerolepmt2_8_20, bgunc_zerolepmt2_8_20,data_zerolepmt2_8_20,
            'combined', calculate_r)
    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8', 16, 19.5, 
            bg_lp8_20b_all, bgunc_lp8_20b_all,data_lp8_20b_all,'strongest',
            calculate_r)
    ss8_20b = j.SS8high('SS_analysis20', 'CMS8', 24, 19.5, bg_ss8HighPt, 
            bgunc_ss8HighPt, data_ss8HighPt, 'combined', calculate_r)
    os5 = j.CmsOs5Fb('OS_analysis5','CMS7', 6, 4.98, bg_os5, bgunc_os5, 
            data_os5, 'combined', calculate_r)
    ge3lp8_20b = j.Cms3Lepton20Fb('GE3LP_analysis20', 'CMS8', 192, 19.5, 
            bg_cms3l8, bgunc_cms3l8, data_cms3l8, 'combined', calculate_r)

    
    mgr = j.AnalysisManager(outdir, gen_info, use_event_weights) #bool for geninfo
    
    #7 TeV searches
    if alphat7bb:
         mgr.add(alphat7_5b)
    #8 TeV searches
    if monojet20b:
        mgr.add(monojet20)
    if mt220b:
        mgr.add(mt2_20b)
    if alphat20b:
        mgr.add(alphat8_20b)
    if os5b:
         mgr.add(os5)
    if lp20b:
        mgr.add(lp8_20b)
    if ss820b:
        mgr.add(ss8_20b)
    if ge3lp20b:
        mgr.add(ge3lp8_20b)

    filemap = filemap_from_dict(filemap_dict, 'jaf')
    mgr.run(FileMapVector([filemap]), "")
    mgr.limit(0.20, write_stats_file, do_combo, calculate_r_combo)
    mgr.write()
  
# FIXME: for reference these commented out lines are kept... Remove them at some point
#    atlas5 = j.ATLAS5('ATLAS5 analysis','ATLAS', 11, 4.710, bg_atlas5, bgunc_atlas5, data_atlas5, 'individual', calculate_r)
#    lp5 = j.LP('LP_analysis5','CMS7', 18, 4.7, bg_lp5, bgunc_lp5, data_lp5, 'combined', calculate_r)
#    monojet5 = j.MonoJet('MonoJet_analysis75','CMS', 4, 5.0, bg_monojet5, bgunc_monojet5, data_monojet5, 'individual', calculate_r)
#    ss5 = j.SS('SS_analysis5','CMS7', 3, 4.98, bg_ss5, bgunc_ss5, data_ss5, 'combined', calculate_r)
    
    # 8 TeV analyses
#    alphat8_12b = j.AlphaTb8('alphaTb8_analysis12','CMS8', 59, 11.7, bg_at8b, bgunc_at8b, data_at8b, 'combined', calculate_r)
#    ss8_11b = j.SSb8('SSb8_analysis11','CMS8', 1, 10.5, bg_ssb8, bgunc_ssb8, data_ssb8, 'combined', calculate_r)
#    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8',16,19.5,bg_lp8_20b,bgunc_lp8_20b,data_lp8_20b,'combined',calculate_r)
#    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8',16,19.5,bg_lp8_20b,bgunc_lp8_20b,data_lp8_20b,'individual',calculate_r)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',52,19.5,bg_zerolepmt2_8_20_stop_test,
#            bgunc_zerolepmt2_8_20_stop_test,data_zerolepmt2_8_20_stop_test,'individual',calculate_r)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',25,19.5,bg_zerolepmt2_8_20_gluino_test,
#            bgunc_zerolepmt2_8_20_gluino_test,data_zerolepmt2_8_20_gluino_test,'combined',calculate_r)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',87,19.5,bg_zerolepmt2_8_20_t2qq,
#            bgunc_zerolepmt2_8_20_t2qq,data_zerolepmt2_8_20_t2qq,'combined',calculate_r)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',25,19.5,bg_zerolepmt2_8_20_gluino_test,
#            bgunc_zerolepmt2_8_20_gluino_test,data_zerolepmt2_8_20_gluino_test,'individual',calculate_r)
#    mt2_20b = j.ZeroLepMt2('MT2_analysis20','CMS8',123,19.5,bg_zerolepmt2_8_20,
#            bgunc_zerolepmt2_8_20,data_zerolepmt2_8_20,'individual',calculate_r)
#    data_lp8_20b_t2tt=IntVector([int(round(x)) for x in bg_lp8_20b_t2tt])
#    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8',8,19.5,bg_lp8_20b_t2tt,
#            bgunc_lp8_20b_t2tt,data_lp8_20b_t2tt,'individual',calculate_r)
#    data_lp8_20b_t2bbww=IntVector([int(round(x)) for x in bg_lp8_20b_t2bbww])
#    lp8_20b = j.CmsSingleLepton20Fb('LP_analysis20','CMS8',8,19.5,bg_lp8_20b_t2bbww,
#            bgunc_lp8_20b_t2bbww,data_lp8_20b_t2bbww,'individual',calculate_r)
#    ss8_20b = j.SS8high('SS_analysis20','CMS8',8,19.5,bg_ss8HighPtBtag2,
#            bgunc_ss8HighPtBtag2,data_ss8HighPtBtag2,'combined',calculate_r)
#    data_cms3l8=IntVector([int(round(x)) for x in bg_cms3l8])

    
    # for extrapolations:
#    bg_at8b40 = DoubleVector([x * 3.4 for x in bg_at8b])
#    data_at8b40 = IntVector([int(x * 3.4) for x in data_at8b])
    
#    bg_ssb40 = DoubleVector([x * 3.8 for x in bg_ssb8])
#    data_ssb40 = IntVector([int(x * 3.8) for x in data_ssb8])
    
    # weather forecast (40/fb):
#    alphat8_40b = j.AlphaTb8('alphaTb8_analysis40','CMS8', 59, 40.0, bg_at8b40, bgunc_at8b, data_at8b40, 'combined', calculate_r)
#    ss8_40b = j.SSb8('SSb8_analysis40','CMS8', 1, 10.5, bg_ssb40, bgunc_ssb8, data_ssb40, 'combined', calculate_r)

    # 20/fb 8 TeV searches
#    ss8HighPt = j.SS8high('SS8high_analysis20','CMS8', 24, 19.5, bg_ss8HighPt, bgunc_ss8HighPt, data_ss8HighPt, 'combined', calculate_r)
#    ss8LowPt = j.SS8low('SS8low_analysis20','CMS8', 24, 19.5, bg_ss8LowPt, bgunc_ss8LowPt, data_ss8LowPt, 'combined', calculate_r)
#    zerolep8 = j.ZeroLep8('SSb8_analysis40','CMS8', 48, 19.4, bg_zerolep8, bgunc_zerolep8, data_zerolep8, 'combined', calculate_r)
