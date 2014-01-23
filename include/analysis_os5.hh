#ifndef __ANALYSISOS__
#define __ANALYSISOS__

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

class OS : public AnalysisBase {

public:
  OS(const std::string & name, 
     const std::string & experiment,
     const unsigned int & numBins);
  
  OS(const std::string & name, 
     const std::string & experiment, 
     const unsigned int & numBins,
     const double & intlumi, 
     //const std::vector<int> & datayields,
     const std::vector<double> & bgpred);
  
  OS(const std::string & name, 
     const std::string & experiment, 
     const unsigned int & numBins,
     const double & intlumi, 
     const std::vector<double> & bgpred,
     const std::vector<double> & bgpreduncert,
     const std::vector<int> & datayields,
     const std::string & fitmode,
     const bool & calculateR);

  ~OS();

  void Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight);
  void initHistos();

private:
  double getht(const std::vector<jjet> & jets);
  double getmet(const TClonesArray *ETMISS);
  bool checkDRjetlep(const std::vector<jjet> & jets, const std::vector<jlepton> & leptons, const double & DRmin);
  std::vector<jjet> getjets(const TClonesArray *JET, const float & pt, const float & eta);
  std::vector<jjet> getjets2(const TClonesArray *JET, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim);
  std::vector<jlepton> getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam, bool poscharge);

  TH1D * ptleadingleppos;
  TH1D * ptleadinglepneg;
  TH2D * ptleadinglepposvsneg;
  TH1D * leadingjetpt;
  TH1D * deltarjetlep;
  TH1D * njethist;
  TH1D * hthist;
  TH1D * methist;
  TH2D * htvsmetsr1sf;
  TH2D * htvsmetsr2sf;
  TH2D * htvsmetsr3sf;
  TH2D * htvsmetsr1of;
  TH2D * htvsmetsr2of;
  TH2D * htvsmetsr3of;
  TH1D * sfinvmass;

};

#endif
