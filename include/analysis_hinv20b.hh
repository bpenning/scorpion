#ifndef __ANALYSISHINV20B__
#define __ANALYSISHINV20B__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

//!!Review all below here
#include "alphat_functions.hh"
#include "jad_jet_functions.hh"
#include "jad_lepton_functions.hh"
#include "jad_track_functions.hh"
#include "jad_photon_functions.hh"
#include "analysis_base.hh"
#include "fileobject_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"


class Hinv20b : public AnalysisBase {

public:
  Hinv20b(const std::string & name, 
	 const std::string & experiment,
	 const unsigned int & numBins);
  
  Hinv20b(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 //const std::vector<int> & datayields,
	 const std::vector<double> & bgpred);

  Hinv20b(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 const std::vector<double> & bgpred,
	 const std::vector<double> & bgpreduncert,
	 const std::vector<int> & datayields,
	 const std::string & fitmode,
	 const bool & calculateR);

  ~Hinv20b();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:

  TH1D * njets;
  TH1D * njets_cjv;
  TH1D * nelectrons;
  TH1D * nmuons;
  TH1D * jet1pt;
  TH1D * jet2pt;
  TH1D * jet3pt;
  TH1D * jet1eta;
  TH1D * jet2eta;
  TH1D * jet3eta;
  TH1D * jet1phi;
  TH1D * jet2phi;
  TH1D * jet3phi;
  TH1D * jetmet_mindphi;
  TH1D * met;
  TH1D * metsignificance;
  TH1D * deltaphijj;
  TH1D * deltaetajj;
  TH1D * mjj;
  TH1D * event_weight;
};

#endif
