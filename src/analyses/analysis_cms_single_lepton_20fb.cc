#include "analysis_cms_single_lepton_20fb.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
CmsSingleLepton20Fb::CmsSingleLepton20Fb(const std::string & name, 
		 const std::string & experiment,
		 const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
CmsSingleLepton20Fb::CmsSingleLepton20Fb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 //const std::vector<int> & datayields,
		 const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

CmsSingleLepton20Fb::CmsSingleLepton20Fb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 const std::vector<double> & bgpred,
		 const std::vector<double> & bgpreduncert,
		 const std::vector<int> & datayields,
		 const std::string & fitmode,
		 const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


CmsSingleLepton20Fb::~CmsSingleLepton20Fb() {}

TSimpleArray<TRootJet> CmsSingleLepton20Fb::SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta) {

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    //check if any jet has Pt>50 and |eta|<3.0
    if(jet->PT > pt && fabs(jet->Eta) < eta) {
      array.Add(jet);
    }
  }
  return array;
}

TSimpleArray<TRootElectron> CmsSingleLepton20Fb::SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta) {
  TIter itElec((TCollection*)ELEC);
  TRootElectron *elec;
  itElec.Reset();
  TSimpleArray<TRootElectron> array;
  while( (elec = (TRootElectron*) itElec.Next()) )
    {
      if (fabs(elec->Eta)<eta && elec->PT > pt && elec->IsolFlag  ){
          array.Add(elec);
      }
    }
  return array;
}

TSimpleArray<TRootMuon> CmsSingleLepton20Fb::SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta) {
  TIter itMuon((TCollection*)MUON);
  TRootMuon *muon;
  itMuon.Reset();
  TSimpleArray<TRootMuon> array;
  while( (muon = (TRootMuon*) itMuon.Next()) )
    {
      if (fabs(muon->Eta)<eta && muon->PT > pt && muon->IsolFlag  ){
        array.Add(muon);
      }
    }
  return array;
}

TSimpleArray<TRootETmis> CmsSingleLepton20Fb::makeETM(const TClonesArray *ETMISS) {
  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();
  TSimpleArray<TRootETmis> array;
  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
    array.Add(etm);
  }
  return array;
}

//bool CmsSingleLepton20Fb::checkforsecondjetphi(const TSimpleArray<TRootJet> & goodjets, const double & dphisep) {
//
//  if(goodjets.GetEntries() == 2) {
//    
//    //check if deltaphi j1, j2 < 2.5
//    double phijet1 = TMath::ATan2(goodjets[0]->Py, goodjets[0]->Px); //calculate jet phi
//    double phijet2 = TMath::ATan2(goodjets[1]->Py, goodjets[1]->Px); //calculate jet phi
//    double dphij12 = fabs(phijet1 - phijet2); //calculate the distance in phi between the jet total momentum
//    
//    //for two vectors in a 2D plane, they cannot be more than pi away from each other
//    if(dphij12 > TMath::Pi()) {
//      dphij12 = fabs(dphij12 - (2.0 * TMath::Pi()));
//    }
//    
//    if (dphij12 < dphisep) {
//      return true;
//    } else {
//      return false;
//    }
//  }
//  
//  return false;
//   
//}

void CmsSingleLepton20Fb::initHistos() {
  andir->cd();
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  calomet = new TH1D("calomet",";E_{T}^{miss} [GeV];Entries",200,-5.,1995.);
}

bool CmsSingleLepton20Fb::gt_1_btag(const TSimpleArray<TRootJet> & jets ){
  bool gt_1_btag=false;
  for(unsigned int ijet = 0; ijet < jets.GetEntries(); ijet++) {
    if (jets[ijet]->Btag){
      gt_1_btag=true;
      break;
    }
  }
  return gt_1_btag;
}

double get_mt2w(const TSimpleArray<TRootElectron> & elecs, const TSimpleArray<TRootMuon> & muons,
        const TSimpleArray<TRootJet> & jets, const TSimpleArray<TRootETmis> & etmis ){
  //inputs from jets
  int njets=jets.GetEntries();
  std::vector<LorentzVector> jet_momenta(njets);
  std::vector<bool> btags(njets);
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    jet_momenta[ijet]=LorentzVector(jets[ijet]->E,jets[ijet]->Px,jets[ijet]->Py,jets[ijet]->Pz);
    btags[ijet]=jets[ijet]->Btag;
  }
  //inputs from leptons
  double lepton_e,lepton_px,lepton_py,lepton_pz;
  if (elecs.GetEntries()==1){
    lepton_e=elecs[0]->E; lepton_px=elecs[0]->Px; lepton_py=elecs[0]->Py; lepton_pz=elecs[0]->Pz; 
  } else if (muons.GetEntries()==1){
    lepton_e=muons[0]->E; lepton_px=muons[0]->Px; lepton_py=muons[0]->Py; lepton_pz=muons[0]->Pz; 
  } else{
    std::cerr << "ERROR: LEPTON SELECTION WASN'T DONE PROPERLY" << std::endl;
    return -123456789;
  }
  LorentzVector lepton_momentum(lepton_e,lepton_px,lepton_py,lepton_pz);
  //inputs from missing ET
  double met = etmis[0]->ET;
  double metphi = etmis[0]->Phi; 
  //calculate MT2W 
  return calculateMT2w(jet_momenta, btags, lepton_momentum, met, metphi);
}

void CmsSingleLepton20Fb::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  TSimpleArray<TRootElectron> goodelecs=SubArrayEl(treereader.Elec(), 30.0, 1.4442); //isolated electrons, pt > 30 GeV, |eta| < 1.4442 (barrel ECAL)
  TSimpleArray<TRootElectron> allelecs=SubArrayEl(treereader.Elec(), 5.0, 1.4442); //isolated electrons, pt > 5 GeV, |eta| < 1.4442 (barrel ECAL)
  TSimpleArray<TRootMuon> goodmuons=SubArrayMu(treereader.Muon(), 25.0, 2.1); //the central isolated muons, pt > 25 GeV, |eta| < 2.1
  TSimpleArray<TRootMuon> allmu=SubArrayMu(treereader.Muon(), 5.0, 2.1); //the central isolated muons, pt > 25 GeV, |eta| < 2.1
  TSimpleArray<TRootJet> goodjets=SubArrayGoodJets(treereader.Jet(), 30.0, 2.4); //check pt > 30 GeV, |eta|<2.4
  TSimpleArray<TRootETmis> etmis=makeETM(treereader.ETMis()); //Missing transverse energy array

  double calo_met = etmis[0]->ET;

  //preselection
  bool preselected=false;
  // one isolated lepton
  if (goodelecs.GetEntries()+goodmuons.GetEntries()==1){
    // no second isolated lepton of pt>5GeV. FIXME: isolation algorithm should be: pTsum<0.2*pT && DeltaR<0.4
    if (allelecs.GetEntries()+allmu.GetEntries()==1){
      // FIXME: no additional isolated track of pt>10 GeV with opposite sign wrt primary lepton
      if (true){
        //FIXME: no jet of pt>20 GeV "consistent with hadronic tau lepton" 
        if (true){
          // at least four jets, at least one btag
          if (goodjets.GetEntries()>=4 && gt_1_btag(goodjets)){
            // MET>100 GeV
            if (calo_met>100){
              preselected=true;
            }
          }
        }
      }
    }
  }

  if (preselected){
    double mt2w=get_mt2w(goodelecs,goodmuons,goodjets,etmis);
  }
    

  return;
}
