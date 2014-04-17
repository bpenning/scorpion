#ifndef __ANALYSISALPHAT__
#define __ANALYSISALPHAT__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "alphat_functions.hh"
#include "analysis_base.hh"
#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
#include "jad_jet_functions.hh"
#include "jad_lepton_functions.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"


class AlphaT : public AnalysisBase {

public:
  AlphaT(const std::string & name, 
	 const std::string & experiment,
	 const unsigned int & numBins);
  
  AlphaT(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 //const std::vector<int> & datayields,
	 const std::vector<double> & bgpred);

  AlphaT(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 const std::vector<double> & bgpred,
	 const std::vector<double> & bgpreduncert,
	 const std::vector<int> & datayields,
	 const std::string & fitmode,
	 const bool & calculateR);

  ~AlphaT();

  void Run(const Reader * treereader, const Reader *gentreereader, const double & weight);
  void initHistos();

private:

  TH1D * leadingjetpt;
  TH1D * hthist;
  TH1D * mhthist;
  TH1D * calomethist;
  TH1D * athist;
  TH1D * njets;
  TH1D * mht_over_ht;
  TH1D * btagrate;
  TH2D * calomet_vs_mht;
  TH2D * ht_vs_mht_pre_alphaT;
  TH2D * ht_vs_mht_post_alphaT;

};

#endif
