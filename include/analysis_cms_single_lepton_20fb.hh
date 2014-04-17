#ifndef __ANALYSIS_CMS_SINGLE_LEPTON_20FB__
#define __ANALYSIS_CMS_SINGLE_LEPTON_20FB__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "analysis_base.hh"
#include "fileobject_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "cms_single_lepton_20fb_functions.hh"

///based on CMS-SUS-13-011; arXiv:1308.1586v2
class CmsSingleLepton20Fb : public AnalysisBase {

public:
  CmsSingleLepton20Fb(const std::string & name, 
	  const std::string & experiment,
	  const unsigned int & numBins);
  
  CmsSingleLepton20Fb(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  //const std::vector<int> & datayields,
	  const std::vector<double> & bgpred);
  
  CmsSingleLepton20Fb(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  const std::vector<double> & bgpred,
	  const std::vector<double> & bgpreduncert,
	  const std::vector<int> & datayields,
	  const std::string & fitmode,
	  const bool & calculateR);

  ~CmsSingleLepton20Fb();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:

  TSimpleArray<TRootJet> SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta);
  TSimpleArray<TRootElectron> SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta);
  TSimpleArray<TRootMuon> SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta);
  TSimpleArray<TRootETmis> makeETM(const TClonesArray *ETMISS);

//  bool checkforsecondjetphi(const TSimpleArray<TRootJet> & goodjets, const double & dphisep);
  bool gt_1_btag(const TSimpleArray<TRootJet> & jets );
  double get_mt2w(const TSimpleArray<TRootElectron> & elecs, const TSimpleArray<TRootMuon> & muons,
        const TSimpleArray<TRootJet> & jets, const TSimpleArray<TRootETmis> & etmis );

  TH1D * leadingjetpt;
  TH1D * calomet;
};

#endif
