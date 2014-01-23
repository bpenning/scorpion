#ifndef __ANALYSISALPHATB__
#define __ANALYSISALPHATB__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "alphat_functions.hh"
#include "analysis_base.hh"
#include "fileobject_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


class AlphaTb : public AnalysisBase {

public:
  AlphaTb(const std::string & name, 
	 const std::string & experiment,
	 const unsigned int & numBins);
  
  AlphaTb(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 //const std::vector<int> & datayields,
	 const std::vector<double> & bgpred);

  AlphaTb(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 const std::vector<double> & bgpred,
	 const std::vector<double> & bgpreduncert,
	 const std::vector<int> & datayields,
	 const std::string & fitmode,
	 const bool & calculateR);

  ~AlphaTb();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:
  TSimpleArray<TRootJet> SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta);
  TSimpleArray<TRootElectron> SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta);
  TSimpleArray<TRootMuon> SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta);
  TSimpleArray<TRootJet> SubArrayBadJets(const TClonesArray *JET, const float & pt, const float & eta);
  TSimpleArray<TRootETmis> makeETM(const TClonesArray *ETMISS);

  TH1D * leadingjetpt;
  TH1D * hthist;
  TH1D * mhthist;
  TH1D * calomethist;
  TH1D * athist;
  TH1D * njets;
  TH1D * mht_over_ht;
  TH1D * btagrate;
  TH2D * calomet_vs_mht;
  TH2D * ht_vs_mht_pre_alphaT;
  TH2D * ht_vs_mht_post_alphaT;

};

#endif
