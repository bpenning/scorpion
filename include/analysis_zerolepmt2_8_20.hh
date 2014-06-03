#ifndef __ANALYSISZEROLEPMT2_8_20__
#define __ANALYSISZEROLEPMT2_8_20__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "alphat_functions.hh"
#include "jad_jet_functions.hh"
#include "jad_lepton_functions.hh"
#include "analysis_base.hh"
#include "fileobject_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "zerolepmt2_functions.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "mt2_bisect.hh"


class ZeroLepMt2 : public AnalysisBase {
  
public:
  ZeroLepMt2(const std::string & name, 
	     const std::string & experiment,
	     const unsigned int & numBins);
  
  ZeroLepMt2(const std::string & name, 
	     const std::string & experiment, 
	     const unsigned int & numBins,
	     const double & intlumi, 
	     //const std::vector<int> & datayields,
	     const std::vector<double> & bgpred);
  
  ZeroLepMt2(const std::string & name, 
	     const std::string & experiment, 
	     const unsigned int & numBins,
	     const double & intlumi, 
	     const std::vector<double> & bgpred,
	     const std::vector<double> & bgpreduncert,
	     const std::vector<int> & datayields,
	     const std::string & fitmode,
	     const bool & calculateR);

  ~ZeroLepMt2();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:

  //TH1D * leadingjetpt;
  TH1D * hthist;
  //TH1D * mhthist;
  TH1D * methist;
  //TH1D * athist;
  //TH1D * njets;
  //TH1D * mht_over_ht;
  TH1D * cut_flow;
  TH1D * btagrate;
  TH2D * met_vs_mht;
  TH2D * low_ht_met_vs_mt2;

};

#endif
