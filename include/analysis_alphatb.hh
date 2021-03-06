#ifndef __ANALYSISALPHATB__
#define __ANALYSISALPHATB__

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
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


class AlphaTb : public AnalysisBase {

public:
  AlphaTb(const std::string & name, 
	 const std::string & experiment,
	 const unsigned int & numBins);
  
  AlphaTb(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 //const std::vector<int> & datayields,
	 const std::vector<double> & bgpred);

  AlphaTb(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 const std::vector<double> & bgpred,
	 const std::vector<double> & bgpreduncert,
	 const std::vector<int> & datayields,
	 const std::string & fitmode,
	 const bool & calculateR);

  ~AlphaTb();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:

  TH1I * cut_sel;
  TH1D * leadingjetpt;
  TH1D * hthist;
  TH1D * mhthist;
  TH1D * calomethist;
  TH1D * athist;
  TH1D * athist2jets;
  TH1D * athist3jets;
  TH1D * athist4jets;
  TH1D * athist5jets;
  TH1D * njets;
  TH1D * ejets;
  TH1D * bjets;
  TH1D * mjets;
  TH1D * mht_over_ht;
  TH1D * btagrate;
  TH2D * calomet_vs_mht;
  TH2D * ht_vs_mht_pre_alphaT;
  TH2D * ht_vs_mht_post_alphaT;
  TH1D * jet1pt;
  TH1D * jet2pt;
  TH1D * deltaphi;

};

#endif
