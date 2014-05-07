#ifndef __ANALYSISZEROLEP8__
#define __ANALYSISZEROLEP8__

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


class ZeroLep8 : public AnalysisBase {
  
public:
  ZeroLep8(const std::string & name, 
	     const std::string & experiment,
	     const unsigned int & numBins);
  
  ZeroLep8(const std::string & name, 
	     const std::string & experiment, 
	     const unsigned int & numBins,
	     const double & intlumi, 
	     //const std::vector<int> & datayields,
	     const std::vector<double> & bgpred);
  
  ZeroLep8(const std::string & name, 
	     const std::string & experiment, 
	     const unsigned int & numBins,
	     const double & intlumi, 
	     const std::vector<double> & bgpred,
	     const std::vector<double> & bgpreduncert,
	     const std::vector<int> & datayields,
	     const std::string & fitmode,
	     const bool & calculateR);

  ~ZeroLep8();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
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

};

#endif
