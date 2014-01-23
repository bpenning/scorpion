#ifndef __ANALYSISLP__
#define __ANALYSISLP__

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


class LP : public AnalysisBase {

public:
  LP(const std::string & name, 
     const std::string & experiment,
     const unsigned int & numBins);
  
  LP(const std::string & name, 
     const std::string & experiment, 
     const unsigned int & numBins,
     const double & intlumi, 
     //const std::vector<int> & datayields,
     const std::vector<double> & bgpred);
  
  LP(const std::string & name, 
     const std::string & experiment, 
     const unsigned int & numBins,
     const double & intlumi, 
     const std::vector<double> & bgpred,
     const std::vector<double> & bgpreduncert,
     const std::vector<int> & datayields,
     const std::string & fitmode,
     const bool & calculateR);

  ~LP();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:

  TSimpleArray<TRootJet> SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta);
  TSimpleArray<TRootElectron> SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta, const float & riso);
  TSimpleArray<TRootMuon> SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta, const float & riso);
  TSimpleArray<TRootETmis> makeETM(const TClonesArray *ETMISS);

  TH1D * leadingjetpt;
  TH1D * num_mu_veto;
  TH1D * num_ele_veto;
  TH1D * hthist;
  TH1D * lphist_el;
  TH1D * sthist_el;
  TH1D * mthist_el;
  TH1D * lphist_mu;
  TH1D * sthist_mu;
  TH1D * mthist_mu;
};

#endif
