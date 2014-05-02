#include "analysis_cms_3_lepton_20fb.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
Cms3Lepton20Fb::Cms3Lepton20Fb(const std::string & name, 
		   const std::string & experiment,
		   const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
Cms3Lepton20Fb::Cms3Lepton20Fb(const std::string & name, 
		   const std::string & experiment, 
		   const unsigned int & numBins,
		   const double & intlumi, 
		   const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

Cms3Lepton20Fb::Cms3Lepton20Fb(const std::string & name, 
		   const std::string & experiment, 
		   const unsigned int & numBins,
		   const double & intlumi, 
		   const std::vector<double> & bgpred,
		   const std::vector<double> & bgpreduncert,
		   const std::vector<int> & datayields,
		   const std::string & fitmode,
		   const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}

Cms3Lepton20Fb::~Cms3Lepton20Fb() {}

void Cms3Lepton20Fb::initHistos() {
  andir->cd();
}


void Cms3Lepton20Fb::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {
  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> ele=goodleptons(treereader->GetElec(), 10.0, 2.4,2.5,2.6);         //the central isolated electrons, pt > PT_ELEC GeV
  std::vector<jlepton> mu=goodleptons(treereader->GetMuon(), 10.0, 2.1,2.5,2.6);          //the central isolated muons, pt > PT_MUON GeV
  std::vector<jjet> badjets=badjetsSkim(treereader->GetJet(),50.0,3.0);   //check for jets which we should veto
  std::vector<jjet> goodjets275=goodjetsSkim(treereader->GetJet(),37.0,3.0); //check for jets which we should analyse

  return;
}
