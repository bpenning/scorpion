#ifndef __ANALYSISALPHAT20B__
#define __ANALYSISALPHAT20B__

#include <iostream>
#include <string>
#include <vector>
#include <numeric>

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
#include "EcalMap.hh"


class AlphaT20b : public AnalysisBase {

public:
  AlphaT20b(const std::string & name, 
	 const std::string & experiment,
	 const unsigned int & numBins);
  
  AlphaT20b(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 //const std::vector<int> & datayields,
	 const std::vector<double> & bgpred);

  AlphaT20b(const std::string & name, 
	 const std::string & experiment, 
	 const unsigned int & numBins,
	 const double & intlumi, 
	 const std::vector<double> & bgpred,
	 const std::vector<double> & bgpreduncert,
	 const std::vector<int> & datayields,
	 const std::string & fitmode,
	 const bool & calculateR);

  ~AlphaT20b();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();
  double biasedDPhiMin(std::vector<jjet>  inJets);

private:

  TH1D * cut_sel;
  TH1D * cut_sel23b0;
  TH1D * cut_selGe4b0;
  TH1D * cut_sel23b1;
  TH1D * cut_selGe4b1;
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
  TH1D * boostHisto;
  TH1D * ejets;
  TH1D * bDPhiHisto23At53b0EcalCut;
  TH1D * bDPhiHistoGe4At53b0EcalCut;
  TH1D * bDPhiHisto23At55b0EcalCut;
  TH1D * bDPhiHistoGe4At55b0EcalCut;
  TH1D * bDPhiHisto23At53b1EcalCut;
  TH1D * bDPhiHistoGe4At53b1EcalCut;
  TH1D * bDPhiHisto23At55b1EcalCut;
  TH1D * bDPhiHistoGe4At55b1EcalCut;
  TH1D * bDPhiHisto23At53b0;
  TH1D * bDPhiHistoGe4At53b0;
  TH1D * bDPhiHisto23At55b0;
  TH1D * bDPhiHistoGe4At55b0;
  TH1D * bDPhiHisto23At53b1;
  TH1D * bDPhiHistoGe4At53b1;
  TH1D * bDPhiHisto23At55b1;
  TH1D * bDPhiHistoGe4At55b1;

  TH1D * bDPhiHisto23At53;
  TH1D * bDPhiHistoGe4At53;
  TH1D * bDPhiHisto23At53EcalCut;
  TH1D * bDPhiHistoGe4At53EcalCut;
  TH1D * bjets;
  TH1D * cTrack;
  TH1D * mjets;
  TH1D * mht_over_ht;
  TH1D * btagrate;
  TH2D * calomet_vs_mht;
  TH2D * ht_vs_mht_pre_alphaT;
  TH2D * ht_vs_mht_post_alphaT;
  TH2D * ecal_map;
  TH1D * jet1pt;
  TH1D * jet2pt;
  TH1D * deltaphi;
  TH1D * event_weight;
};

#endif
