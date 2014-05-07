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

  //miscellaneous
  bool gt_1_btag(const TSimpleArray<TRootJet> & jets );
  double dphi(double, double);
  //copies of the lepton and the leading bjet
  TRootParticle get_lepton(const TSimpleArray<TRootElectron> & elecs, const TSimpleArray<TRootMuon> & muons);
  TRootJet get_leading_bjet(const TSimpleArray<TRootJet> & jets);
  //kinematic variables
  double get_mt2w(const TRootParticle & lepton, const TSimpleArray<TRootJet> & jets, const TSimpleArray<TRootETmis> & etmis );
  double get_chi2(const TSimpleArray<TRootJet> & jets);
  double get_min_dphi(const TSimpleArray<TRootJet> & jets,const TSimpleArray<TRootETmis> & etmis);
  double get_htratio(const TSimpleArray<TRootJet> & jets,const TSimpleArray<TRootETmis> & etmis);
  double get_delta_R(const TRootParticle & lepton, const TRootJet & leading_bjet);
  double get_mt(const TRootParticle & lepton, const TSimpleArray<TRootETmis> & etmis);

//  TH1D * leadingjetpt;
//  TH1D * calomet;
  TH1D * cut_flow_hist;
  TH1D * mt_hist;
  TH1D * met_hist;
  TH1D * mt2w_hist;
  TH1D * chi2_hist;
  TH1D * htratio_hist;
  TH1D * min_dphi_hist;
  TH1D * leading_bjet_pt_hist;
  TH1D * delta_R_hist;
  TH1D * lepton_pt_hist;
};

#endif
