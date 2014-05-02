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

TSimpleArray<TRootJet> Cms3Lepton20Fb::SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta) {

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    if(jet->PT > pt && fabs(jet->Eta) < eta) {
      array.Add(jet);
    }
  }
  return array;
}

TSimpleArray<TRootElectron> Cms3Lepton20Fb::SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta) {
  TIter itElec((TCollection*)ELEC);
  TRootElectron *elec;
  itElec.Reset();
  TSimpleArray<TRootElectron> array;
  while( (elec = (TRootElectron*) itElec.Next()) )
    {
      if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
      array.Add(elec);
    }
  return array;
}

TSimpleArray<TRootMuon> Cms3Lepton20Fb::SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta) {
  TIter itMuon((TCollection*)MUON);
  TRootMuon *muon;
  itMuon.Reset();
  TSimpleArray<TRootMuon> array;
  while( (muon = (TRootMuon*) itMuon.Next()) )
    {
      if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
      array.Add(muon);
    }
  return array;
}

TSimpleArray<TRootTauJet> Cms3Lepton20Fb::SubArrayTauJet(const TClonesArray *TAUJET, const float & pt, const float & eta) {
  TIter itTauJet((TCollection*)TAUJET);
  TRootTauJet *taujet;
  itTauJet.Reset();
  TSimpleArray<TRootTauJet> array;
  while( (taujet = (TRootTauJet*) itTauJet.Next()) )
    {
      if(taujet->PT<pt || !taujet->IsolFlag || fabs(taujet->Eta) > eta) continue;
      array.Add(taujet);
    }
  return array;
}

TSimpleArray<TRootJet> Cms3Lepton20Fb::SubArrayBadJets(const TClonesArray *JET, const float & pt, const float & eta) {

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    //check if any jet has Pt>50 and |eta|>3.0
    if(jet->PT > pt && fabs(jet->Eta) > eta) {
      array.Add(jet);
    }
  }
  return array;
}

TSimpleArray<TRootETmis> Cms3Lepton20Fb::makeETM(const TClonesArray *ETMISS) {
  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();
  TSimpleArray<TRootETmis> array;
  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
    array.Add(etm);
  }
  return array;
}


void Cms3Lepton20Fb::initHistos() {
  andir->cd();
}

void Cms3Lepton20Fb::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  TSimpleArray<TRootElectron> ele=SubArrayEl(treereader.Elec(), 10.0, 2.4);         //the central isolated electrons, pt > PT_ELEC GeV
  TSimpleArray<TRootMuon>     mu=SubArrayMu(treereader.Muon(), 10.0, 2.4);          //the central isolated muons, pt > PT_MUON GeV
  TSimpleArray<TRootJet>      badjets=SubArrayBadJets(treereader.Jet(),50.0,3.0);   //check for jets which we should veto
  TSimpleArray<TRootJet>      goodjets275=SubArrayGoodJets(treereader.Jet(),37.0,3.0); //check for jets which we should analyse
  TSimpleArray<TRootJet>      goodjets325=SubArrayGoodJets(treereader.Jet(),43.0,3.0); //check for jets which we should analyse
  TSimpleArray<TRootJet>      goodjets375=SubArrayGoodJets(treereader.Jet(),50.0,3.0); //check for jets which we should analyse
  TSimpleArray<TRootETmis>    etmis=makeETM(treereader.ETMis()); //Missing transverse energy array

  return;
}
