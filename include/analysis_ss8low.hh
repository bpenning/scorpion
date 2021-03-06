#ifndef __ANALYSISSS8LOW__
#define __ANALYSISSS8LOW__

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
#include "jad_jet_functions.hh"
#include "jad_lepton_functions.hh"
#include "pair_info.hh"


class SS8low : public AnalysisBase {

public:
  SS8low(const std::string & name, 
       const std::string & experiment,
       const unsigned int & numBins);
  
  SS8low(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       //const std::vector<int> & datayields,
       const std::vector<double> & bgpred);
  
  SS8low(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       const std::vector<double> & bgpred,
       const std::vector<double> & bgpreduncert,
       const std::vector<int> & datayields,
       const std::string & fitmode,
       const bool & calculateR);
  
  ~SS8low();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:

  TH1D * njethist;
  TH1D * hthist;
  TH1D * methist;
  TH1D * leadingjetpt;  
  TH1D * sfinvmass;
  TH1D * ptleadinglep;
  TH2D * ptchargeonevschargetwo;
  TH1D * numbjets;
  TH1D * drjetcombo;
};

#endif
