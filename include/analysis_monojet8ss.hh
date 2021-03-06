#ifndef __ANALYSISMONOJET8SS__
#define __ANALYSISMONOJET8SS__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "analysis_base.hh"
#include "jad_jet_functions.hh"
#include "jad_lepton_functions.hh"
#include "analysis_base.hh"
#include "fileobject_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


class MonoJet8ss : public AnalysisBase {

public:
  MonoJet8ss(const std::string & name, 
	  const std::string & experiment,
	  const unsigned int & numBins);
  
  MonoJet8ss(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  //const std::vector<int> & datayields,
	  const std::vector<double> & bgpred);
  
  MonoJet8ss(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  const std::vector<double> & bgpred,
	  const std::vector<double> & bgpreduncert,
	  const std::vector<int> & datayields,
	  const std::string & fitmode,
	  const bool & calculateR);

  ~MonoJet8ss();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:


  bool checkforsecondjetphi(const std::vector<jjet> & goodjets, const double & dphisep);

  TH1D * leadingjetpt;
  TH1D * calomet;
  TH1D * event_weight;

  TH1D * cutflow;
 
  TH1D * njets;

};

#endif
