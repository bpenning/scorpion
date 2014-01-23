#ifndef __ANALYSISSS__
#define __ANALYSISSS__

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

class SS : public AnalysisBase {

public:
  SS(const std::string & name, 
     const std::string & experiment,
     const unsigned int & numBins);
  
  SS(const std::string & name, 
     const std::string & experiment, 
     const unsigned int & numBins,
     const double & intlumi, 
     //const std::vector<int> & datayields,
     const std::vector<double> & bgpred);
  
  SS(const std::string & name, 
     const std::string & experiment, 
     const unsigned int & numBins,
     const double & intlumi, 
     const std::vector<double> & bgpred,
     const std::vector<double> & bgpreduncert,
     const std::vector<int> & datayields,
     const std::string & fitmode,
     const bool & calculateR);

  ~SS();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:
  double getht(const std::vector<jjet> & jets);
  double getmet(const TClonesArray *ETMISS);
  unsigned int getnbtags(const std::vector<jjet> & jets);
  std::vector<jjet> getjets(const TClonesArray *JET, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim);
  std::vector<jlepton> getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam);
  std::vector<jparticle> getparticles(const TClonesArray * GEN, const unsigned int & pidfilter=0);

  TH1D * njethist;
  TH1D * hthist;
  TH1D * methist;
  TH1D * leadingjetpt;  
  TH1D * sfinvmass;
  TH1D * ptleadinglep;
  TH2D * ptchargeonevschargetwo;
  TH1D * numbjets;
  TH1D * jetptb;
  TH1D * bjetptb;
  TH1D * jetptc;
  TH1D * bjetptc;
  TH1D * jetptm;
  TH1D * bjetptm;
  TH1D * drjetcombo;
};

#endif
