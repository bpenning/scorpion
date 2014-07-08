#ifndef __ANALYSIS_OS5__
#define __ANALYSIS_OS5__
#include <vector>

#include "analysis_base.hh"
#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
#include "jad_jet_functions.hh"
#include "jad_lepton_functions.hh"
#include "TH1.h"
#include "TH2.h"

//CMS-SUS-11-011, arXiv:1206.3949
class CmsOs5Fb : public AnalysisBase {

public:
  CmsOs5Fb(const std::string & name, 
	  const std::string & experiment,
	  const unsigned int & numBins);
  
  CmsOs5Fb(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  //const std::vector<int> & datayields,
	  const std::vector<double> & bgpred);
  
  CmsOs5Fb(const std::string & name, 
	  const std::string & experiment, 
	  const unsigned int & numBins,
	  const double & intlumi, 
	  const std::vector<double> & bgpred,
	  const std::vector<double> & bgpreduncert,
	  const std::vector<int> & datayields,
	  const std::string & fitmode,
	  const bool & calculateR);

  ~CmsOs5Fb();

  void Run(const Reader * treereader, const Reader * gentreereader, const double & weight);
  void initHistos();

private:
  TH1D * cut_flow_hist;
  TH1D * ptleadingleppos;
  TH1D * ptleadinglepneg;
  TH2D * ptleadinglepposvsneg;
  TH1D * leadingjetpt;
  TH1D * deltarjetlep;
  TH1D * njethist;
  TH1D * hthist;
  TH1D * methist;
  TH2D * htvsmetsr1sf;
  TH2D * htvsmetsr2sf;
  TH2D * htvsmetsr3sf;
  TH2D * htvsmetsr1of;
  TH2D * htvsmetsr2of;
  TH2D * htvsmetsr3of;
  TH1D * sfinvmass;

};

#endif
