#ifndef __ANALYSIS_CMS_3_LEPTON_20FB__
#define __ANALYSIS_CMS_3_LEPTON_20FB__

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
#include "TSimpleArray.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


struct ossf_bools{
    bool on_Z, above_Z, below_Z;
};
//Analysis: CMS-SUS-13-002 ; arXiv:1404.5801v1
class Cms3Lepton20Fb : public AnalysisBase {

public:
  Cms3Lepton20Fb(const std::string & name, 
	   const std::string & experiment,
	   const unsigned int & numBins);
  
  Cms3Lepton20Fb(const std::string & name, 
	   const std::string & experiment, 
	   const unsigned int & numBins,
	   const double & intlumi, 
	   const std::vector<double> & bgpred);
  
  Cms3Lepton20Fb(const std::string & name, 
	   const std::string & experiment, 
	   const unsigned int & numBins,
	   const double & intlumi, 
	   const std::vector<double> & bgpred,
	   const std::vector<double> & bgpreduncert,
	   const std::vector<int> & datayields,
	   const std::string & fitmode,
	   const bool & calculateR);
  
  ~Cms3Lepton20Fb();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:
  //  bool reject_ossf(const std::vector<jlepton> & leptons,double mlll_min, double mlll_max);
    bool reject_ossf(const std::vector<std::pair<jlepton,jlepton> >& ossf_pairs,
           std::vector<jlepton> leptons);
    //CMS-SUS-13-002 Table 2: '"On-Z" refers to events with at least one e+e- or mu+mu- (OSSF) pair with dilepton mass between
    //75 and 115 GeV, while "Off-Z" refers to events with one or two OSSF pairs, none of which fall in thes mass range'
    //Table 3: '"On-Z" refers to events with an e+e- or mu+mu- (OSSF) pair with dilepton mass between 75 and 105 GeV, while
    //"Above-Z" and "Below-Z" refer to events with an OSSF pair with mass above 105 Gev ofr 75 GeV, respectively'
    ossf_bools get_ossf_bools(const std::vector<std::pair<jlepton,jlepton> > & ossf_pairs,double mlll_min, double mlll_max);
    TH1D * cut_flow_hist;
    TH1D * njets;
    TH1D * nleptons;
    TH1D * nbjets;
    TH1D * nossfpairs;
    TH1D * ntaujets;
    TH1D * met_hist;
    TH1D * ht_hist;
};

#endif
