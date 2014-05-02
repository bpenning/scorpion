#ifndef __ANALYSIS_CMS_3_LEPTON_20FB__
#define __ANALYSIS_CMS_3_LEPTON_20FB__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "analysis_base.hh"
#include "fileobject_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


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

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:
  TSimpleArray<TRootJet> SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta);
  TSimpleArray<TRootElectron> SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta);
  TSimpleArray<TRootMuon> SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta);
  TSimpleArray<TRootJet> SubArrayBadJets(const TClonesArray *JET, const float & pt, const float & eta);
  TSimpleArray<TRootETmis> makeETM(const TClonesArray *ETMISS);
};

#endif
