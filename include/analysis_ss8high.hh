#ifndef __ANALYSISSS8HIGH__
#define __ANALYSISSS8HIGH__

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
#include "jad_particle_functions.hh"
#include "pair_info.hh"


class SS8high : public AnalysisBase {

public:
  SS8high(const std::string & name, 
       const std::string & experiment,
       const unsigned int & numBins);
  
  SS8high(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       //const std::vector<int> & datayields,
       const std::vector<double> & bgpred);
  
  SS8high(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       const std::vector<double> & bgpred,
       const std::vector<double> & bgpreduncert,
       const std::vector<int> & datayields,
       const std::string & fitmode,
       const bool & calculateR);
  
  ~SS8high();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:

  TH1D * njethist;
  TH1D * hthist;
  TH1D * methist;
  TH1D * leadingjetpt;  
  TH1D * leadingbjetpt;  
  TH1D * sfinvmass;
  TH1D * ptleadinglep;
  TH2D * ptchargeonevschargetwo;
  TH1D * numbjets;
  TH1D * drjetcombo;
  TH1D * genelectronshist;
  TH1D * recelectronshist;
  TH1D * genmuonshist;
  TH1D * recmuonshist;
  TH1D * gentaushist;
  TH1D * rectaushist;
  TH1D * genbjetshist;
  TH1D * recbjetshist;
};

#endif
