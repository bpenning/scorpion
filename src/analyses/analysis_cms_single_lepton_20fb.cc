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

void CmsSingleLepton20Fb::initHistos() {
  andir->cd();
  cut_flow_hist = new TH1D("cuts",";cuts;entries",7,0,7);
  //same histograms as in Fig. 2 1308.1586
  mt_hist= new TH1D("mt",";M_{T} [GeV];Entries",10 ,0 ,300);
  met_hist= new TH1D("met",";E_{T}^{miss} [GeV];Entries",10 ,100 ,350);
  mt2w_hist= new TH1D("mt2w",";M_{T2}^W [GeV];Entries",17 ,75 ,500);
  chi2_hist= new TH1D("chi2",";Hadronic~top~#chi^2;Entries", 20,0 ,20);
  htratio_hist= new TH1D("htratio",";H_T^{ratio};Entries",25 ,0 ,1);
  min_dphi_hist= new TH1D("min_dphi",";#min#{#Delta#phi(j_1,E^{miss}_T),#Delta#phi(j_2,E^{miss}_T) #}[rad];Entries", 15,0 ,TMath::Pi());
  leading_bjet_pt_hist= new TH1D("leading_bjet_pt",";p_T(b_T);Entries",9, 30 ,300 );
  delta_R_hist= new TH1D("delta_R",";#Delta R(l,b_1);Entries",15 ,0 ,5);
  lepton_pt_hist= new TH1D("lepton_pt",";p_T(l) [GeV];Entries", 8, 20,100);
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

double CmsSingleLepton20Fb::get_mt2w(const TRootParticle & lepton, const TSimpleArray<TRootJet> & jets, const TSimpleArray<TRootETmis> & etmis ){
  //inputs from jets
  int njets=jets.GetEntries();
  std::vector<LorentzVector> jet_momenta(njets);
  std::vector<bool> btags(njets);
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    jet_momenta[ijet]=LorentzVector(jets[ijet]->E,jets[ijet]->Px,jets[ijet]->Py,jets[ijet]->Pz);
    btags[ijet]=jets[ijet]->Btag;
  }
  //inputs from leptons
  LorentzVector lepton_momentum(lepton.E,lepton.Px,lepton.Py,lepton.Pz);
  //inputs from missing ET
  double met = etmis[0]->ET;
  double metphi = etmis[0]->Phi; 
  //calculate MT2W 
  return calculateMT2w(jet_momenta, btags, lepton_momentum, met, metphi);
}

double CmsSingleLepton20Fb::get_chi2(const TSimpleArray<TRootJet> & jets){
  //inputs from jets
  int njets=jets.GetEntries();
  std::vector<LorentzVector> jet_momenta(njets);
  //FIXME: fix the jet enery resolution to 10% for now. 
  //Verena Martinez Outschoorn should tell how to do it better
  std::vector<float> jet_energy_resolutions(njets,0.1);
  std::vector<bool> btags(njets);
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    jet_momenta[ijet]=LorentzVector(jets[ijet]->E,jets[ijet]->Px,jets[ijet]->Py,jets[ijet]->Pz);
    btags[ijet]=jets[ijet]->Btag;
  }
  return calculateChi2(jet_momenta, jet_energy_resolutions, btags);
}

double CmsSingleLepton20Fb::dphi(double phi1, double phi2){
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

double CmsSingleLepton20Fb::get_min_dphi(const TSimpleArray<TRootJet> & jets,const TSimpleArray<TRootETmis> & etmis){
 double phi_jet1=jets[0]->Phi;
 double phi_jet2=jets[1]->Phi;
 double phi_etmiss=etmis[0]->Phi;
 return std::min(dphi(phi_jet1,phi_etmiss),dphi(phi_jet2,phi_etmiss));
}

//FIXME: not sure if the opposite hemisphere is defined using Phi only
double CmsSingleLepton20Fb::get_htratio(const TSimpleArray<TRootJet> & jets,const TSimpleArray<TRootETmis> & etmis){
  int njets=jets.GetEntries();
  double same_hemisphere_httot;
  double httot;
  double et;
  double phi_etmis=etmis[0]->Phi;
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    et=jets[ijet]->PT;
    httot+=et;
    if (std::abs(jets[ijet]->Phi-phi_etmis)>TMath::Pi()/2){
      same_hemisphere_httot+=et;
    }
  }
  return same_hemisphere_httot/httot;
}

TRootJet CmsSingleLepton20Fb::get_leading_bjet(const TSimpleArray<TRootJet> & jets){
  TRootJet result;
  int njets=jets.GetEntries();
  double e, px, py, pz;
  bool Btag;
  int NTracks, NCalo;
  double EHoverEE;
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    if (jets[ijet]->Btag){
      e=jets[ijet]->E; 
      px=jets[ijet]->Px; 
      py=jets[ijet]->Py; 
      pz=jets[ijet]->Pz; 
      result.Set(px,py,pz,e);
      result.Btag=true;
      result.NTracks=jets[ijet]->NTracks;
      result.NCalo=jets[ijet]->NCalo;
      result.EHoverEE=jets[ijet]->EHoverEE;
      break;
    }
  }
  return result;
}

double CmsSingleLepton20Fb::get_delta_R(const TRootParticle & lepton, const TRootJet & leading_bjet){
  double lepton_eta=lepton.Eta;
  double lepton_phi=lepton.Phi;
  double bjet_eta=leading_bjet.Eta;
  double bjet_phi=leading_bjet.Phi;
  double delta_eta = lepton_eta - bjet_eta;
  double delta_phi = dphi(lepton_phi, bjet_phi);
  return sqrt(delta_eta*delta_eta+delta_phi*delta_phi);
}

double CmsSingleLepton20Fb::get_mt(const TRootParticle & lepton, const TSimpleArray<TRootETmis> & etmis){
  double met_et=etmis[0]->ET;
  double met_phi=etmis[0]->Phi;
  double lepton_pt=lepton.PT;
  double lepton_phi=lepton.Phi;
  return sqrt(2*met_et*lepton_pt*(1-cos(dphi(lepton_phi,met_phi))));
}

TRootParticle CmsSingleLepton20Fb::get_lepton(const TSimpleArray<TRootElectron> & elecs, const TSimpleArray<TRootMuon> & muons){
  TRootParticle result;
  double e, px, py, pz;
  if (elecs.GetEntries()==1){
    e=elecs[0]->E; 
    px=elecs[0]->Px; 
    py=elecs[0]->Py; 
    pz=elecs[0]->Pz; 
  } else if (muons.GetEntries()==1){
    e=muons[0]->E; 
    px=muons[0]->Px; 
    py=muons[0]->Py; 
    pz=muons[0]->Pz; 
  } else{
    std::cerr << "ERROR: LEPTON SELECTION WASN'T DONE PROPERLY" << std::endl;
  }
  result.Set(px,py,pz,e);
  return result;
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
  cut_flow_hist->Fill(0.5);
  if (goodelecs.GetEntries()+goodmuons.GetEntries()==1){
    cut_flow_hist->Fill(1.5);
    // no second isolated lepton of pt>5GeV. FIXME: isolation algorithm should be: pTsum<0.2*pT && DeltaR<0.4
    if (allelecs.GetEntries()+allmu.GetEntries()==1){
      cut_flow_hist->Fill(2.5);
      // FIXME: no additional isolated track of pt>10 GeV with opposite sign wrt primary lepton
      if (true){
        cut_flow_hist->Fill(3.5);
        //FIXME: no jet of pt>20 GeV "consistent with hadronic tau lepton" 
        if (true){
          cut_flow_hist->Fill(4.5);
          // at least four jets, at least one btag
          if (goodjets.GetEntries()>=4 && gt_1_btag(goodjets)){
            cut_flow_hist->Fill(5.5);
            // MET>100 GeV
            if (calo_met>100){
              cut_flow_hist->Fill(6.5);
              preselected=true;
            }
          }
        }
      }
    }
  }

  if (preselected){
    //get single lepton
    TRootParticle lepton=get_lepton(goodelecs,goodmuons);
    TRootJet leading_bjet=get_leading_bjet(goodjets);
    //get 9 kinematic variables
    double mt2w=get_mt2w(lepton,goodjets,etmis);
    double chi2=get_chi2(goodjets);
    double min_dphi=get_min_dphi(goodjets,etmis);
    //FIXME: check if get_htratio has been defined correctly
    double htratio=get_htratio(goodjets,etmis);
    double leading_bjet_pt=leading_bjet.PT;
    double delta_R=get_delta_R(lepton,leading_bjet);
    double lepton_pt=lepton.PT;
    double mt=get_mt(lepton,etmis);
    double met=etmis[0]->ET;
    //fill histograms (cf. Fig. 2)
    mt_hist->Fill(mt,1000*weight);
    met_hist->Fill(met,1000*weight);
    mt2w_hist->Fill(mt2w,1000*weight);
    chi2_hist->Fill(chi2,1000*weight);
    htratio_hist->Fill(htratio,1000*weight);
    min_dphi_hist->Fill(min_dphi,1000*weight);
    leading_bjet_pt_hist->Fill(leading_bjet_pt,1000*weight);
    delta_R_hist->Fill(delta_R,1000*weight);
    lepton_pt_hist->Fill(lepton_pt,1000*weight);
    // stop->top neu; low DeltaM
    if (met>150 && min_dphi>0.8 && chi2<5 ) mSigPred.at(0)+=weight;
    if (met>200 && min_dphi>0.8 && chi2<5 ) mSigPred.at(1)+=weight;
    if (met>250 && min_dphi>0.8 && chi2<5 ) mSigPred.at(2)+=weight;
    if (met>300 && min_dphi>0.8 && chi2<5 ) mSigPred.at(3)+=weight;
    // stop->top neu; high DeltaM
    if (met>150 && mt2w>200 && min_dphi>0.8 && chi2<5 ) mSigPred.at(4)+=weight;
    if (met>200 && mt2w>200 && min_dphi>0.8 && chi2<5 ) mSigPred.at(5)+=weight;
    if (met>250 && mt2w>200 && min_dphi>0.8 && chi2<5 ) mSigPred.at(6)+=weight;
    if (met>300 && mt2w>200 && min_dphi>0.8 && chi2<5 ) mSigPred.at(7)+=weight;
    // stop->bot char; low DeltaM
    if (met>100 && min_dphi>0.8 ) mSigPred.at(8)+=weight;
    if (met>150 && min_dphi>0.8 ) mSigPred.at(9)+=weight;
    if (met>200 && min_dphi>0.8 ) mSigPred.at(10)+=weight;
    if (met>250 && min_dphi>0.8 ) mSigPred.at(11)+=weight;
    // stop->bot char; high DeltaM
    if (met>100 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(12)+=weight;
    if (met>150 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(13)+=weight;
    if (met>200 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(14)+=weight;
    if (met>250 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(15)+=weight;
  }
  return;
}
