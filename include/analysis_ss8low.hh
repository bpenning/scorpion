#ifndef __ANALYSISSSB8__
#define __ANALYSISSSB8__

#include <iostream>
#include <algorithm>
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

#include "jad_lepton_class.hh"
#include "jad_jet_class.hh"
#include "jad_particle_class.hh"
#include "pair_info.hh"

class SS8low : public AnalysisBase {

public:
  SS8low(const std::string & name, 
       const std::string & experiment,
       const unsigned int & numBins);
  
  SS8low(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       //const std::vector<int> & datayields,
       const std::vector<double> & bgpred);
  
  SS8low(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       const std::vector<double> & bgpred,
       const std::vector<double> & bgpreduncert,
       const std::vector<int> & datayields,
       const std::string & fitmode,
       const bool & calculateR);
  
  ~SS8low();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:
  double getht(const std::vector<jjet> & jets);
  double getmet(const TClonesArray *ETMISS);
  unsigned int getnbtags(const std::vector<jjet> & jets);
  std::vector<jjet> getjets(const TClonesArray *JET, const float & pt, const float & eta);
  std::vector<jlepton> getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam);
  

  TH1D * njethist;
  TH1D * hthist;
  TH1D * methist;
  TH1D * leadingjetpt;  
  TH1D * sfinvmass;
  TH1D * ptleadinglep;
  TH2D * ptchargeonevschargetwo;
  TH1D * numbjets;
  TH1D * drjetcombo;
};

#endif