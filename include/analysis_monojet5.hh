#ifndef __ANALYSISMONOJET__
#define __ANALYSISMONOJET__

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


class MonoJet : public AnalysisBase {

public:
  MonoJet(const std::string & name, 
	  const std::string & experiment,
	  const unsigned int & numBins);
  
  MonoJet(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  //const std::vector<int> & datayields,
	  const std::vector<double> & bgpred);
  
  MonoJet(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  const std::vector<double> & bgpred,
	  const std::vector<double> & bgpreduncert,
	  const std::vector<int> & datayields,
	  const std::string & fitmode,
	  const bool & calculateR);

  ~MonoJet();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:

  TSimpleArray<TRootJet> SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta);
  TSimpleArray<TRootElectron> SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta);
  TSimpleArray<TRootMuon> SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta);
  TSimpleArray<TRootETmis> makeETM(const TClonesArray *ETMISS);

  bool checkforsecondjetphi(const TSimpleArray<TRootJet> & goodjets, const double & dphisep);

  TH1D * leadingjetpt;
  TH1D * calomet;
};

#endif
