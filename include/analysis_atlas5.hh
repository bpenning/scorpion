#ifndef __ANALYSISATLAS5__
#define __ANALYSISATLAS5__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "analysis_base.hh"
#include "fileobject_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "TH1.h"
#include "TMath.h"

class ATLAS5 : public AnalysisBase {

public:
  ATLAS5(const std::string & name, 
	 const std::string & experiment,
	 const unsigned int & numBins);
  
  ATLAS5(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 //const std::vector<int> & datayields,
	 const std::vector<double> & bgpred);

  ATLAS5(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 const std::vector<double> & bgpred,
	 const std::vector<double> & bgpreduncert,
	 const std::vector<int> & datayields,
	 const std::string & fitmode,
	 const bool & calculateR);

  ~ATLAS5();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:  

  bool jetdrsep (const TSimpleArray<TRootJet> & jets, const double & DR);
  bool deltaphiptjet(const TSimpleArray<TRootJet> & jets, const double & px, const double & py, const double & deltaphi, const int & jetcuts, double & phiret, double & dphiret);

  TSimpleArray<TRootElectron> SubArrayEl(const TClonesArray * ELEC, float pt, float eta);
  TSimpleArray<TRootMuon> SubArrayMu(const TClonesArray *MUON, float pt, float eta);
  TSimpleArray<TRootJet> SubArrayGoodJets2(const TClonesArray *JET, float pt, float eta, double DRlim, const TSimpleArray<TRootElectron> &ele);
  TSimpleArray<TRootElectron> SubArrayEl2(const TClonesArray * ELEC, float pt, float eta, double DRLmin, const TSimpleArray<TRootJet> &jet);
  TSimpleArray<TRootMuon> SubArrayMu2(const TClonesArray *MUON, float pt, float eta,double DRLmin, const TSimpleArray<TRootJet> &jet);
  TSimpleArray<TRootETmis> makeETM(const TClonesArray *ETMISS);

  TH1D * njets;
  TH1D * PHIPT;
  TH1D * DELTAPHI;
  TH1D * PHI;
  TH1D * PHI40;
  TH1D * hthist;
  TH1D * etmisshist;
  TH1D * ptdiff;
  TH1D * leadingjetpt;
  TH1D * etdiff;
  TH1D * ATLASA;
  TH1D * ATLASAd;
  TH1D * ATLASB;
  TH1D * ATLASC;
  TH1D * ATLASD;
  TH1D * ATLASE;
  
};

#endif
