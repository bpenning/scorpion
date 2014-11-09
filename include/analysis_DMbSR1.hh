#ifndef __ANALYSISDMB_SR1__
#define __ANALYSISDMB_SR1__

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


class DMbSR1 : public AnalysisBase {

public:
  DMbSR1(const std::string & name, 
	  const std::string & experiment,
	  const unsigned int & numBins);
  
  DMbSR1(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  //const std::vector<int> & datayields,
	  const std::vector<double> & bgpred);
  
  DMbSR1(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  const std::vector<double> & bgpred,
	  const std::vector<double> & bgpreduncert,
	  const std::vector<int> & datayields,
	  const std::string & fitmode,
	  const bool & calculateR);

  ~DMbSR1();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:


  bool checkforsecondjetphi(const std::vector<jjet> & goodjets, const double & dphisep);

  TH1D * leadingjetpt;
  TH1D * calomet;

  TH1D * cutflow;
  TH1D * event_weight;
  TH1D * njets;

};

#endif
