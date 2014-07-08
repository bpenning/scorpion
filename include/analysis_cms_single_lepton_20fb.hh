#ifndef __ANALYSIS_CMS_SINGLE_LEPTON_20FB__
#define __ANALYSIS_CMS_SINGLE_LEPTON_20FB__
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <numeric>

#include "analysis_base.hh"
#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
#include "jad_jet_functions.hh"
#include "jad_lepton_functions.hh"
#include "jBlockClasses.hh"
//#include "TSimpleArray.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
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

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:

  //kinematic variables
  double get_mt2w(const jlepton  & lepton, const std::vector<jjet> & jets, const jjet & etmis );
  double get_chi2(const std::vector<jjet> & jets);
  double get_htratio(const std::vector<jjet> & jets, const jjet & etmis);

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
  TH1D * njets;
  TH1D * nbjets;
  TH1D * lepton_flavour;
};

#endif
