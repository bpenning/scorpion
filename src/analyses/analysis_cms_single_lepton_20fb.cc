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


void CmsSingleLepton20Fb::initHistos() {
  andir->cd();
  cut_flow_hist = new TH1D("cuts",";cuts;entries",7,0,7);
  njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
  nbjets = new TH1D("nbjets", ";N_{bjets};Entries",10,-0.5,9.5);
  lepton_flavour = new TH1D("lepton_flavour", ";lepton flavour;Entries",2,0,2);
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

double CmsSingleLepton20Fb::get_mt2w(const jlepton  & lepton, const std::vector<jjet> & jets, const jjet & etmis ){
  //inputs from jets
  int njets=jets.size();
  std::vector<LorentzVector> jet_momenta(njets);
  std::vector<bool> btags(njets);
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    jet_momenta[ijet]=LorentzVector(jets[ijet].Px(),jets[ijet].Py(),jets[ijet].Pz(),jets[ijet].E());
    btags[ijet]=jets[ijet].Btag();
  }
  //inputs from leptons
  LorentzVector lepton_momentum(lepton.Px(),lepton.Py(),lepton.Pz(),lepton.E());
  //inputs from missing ET
  double met = etmis.Et();
  double metphi = etmis.Phi(); 
  //calculate MT2W 
  return calculateMT2w(jet_momenta, btags, lepton_momentum, met, metphi);
}

double CmsSingleLepton20Fb::get_chi2(const std::vector<jjet> & jets){
  //inputs from jets
  int njets=jets.size();
  std::vector<LorentzVector> jet_momenta(njets);
  //FIXME: fix the jet enery resolution to 10% for now. 
  //Verena Martinez Outschoorn should tell how to do it better
  std::vector<float> jet_energy_resolutions(njets,0.1);
  std::vector<bool> btags(njets);
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    jet_momenta[ijet]=LorentzVector(jets[ijet].Px(),jets[ijet].Py(),jets[ijet].Pz(),jets[ijet].E());
    btags[ijet]=jets[ijet].Btag();
  }
  return calculateChi2(jet_momenta, jet_energy_resolutions, btags);
}


//FIXME: not sure if the opposite hemisphere is defined using Phi only
double CmsSingleLepton20Fb::get_htratio(const std::vector<jjet> & jets,const jjet & etmis){
  int njets=jets.size();
  double same_hemisphere_httot;
  double httot;
  double et;
  for(unsigned int ijet = 0; ijet < njets; ijet++) {
    et=jets[ijet].Pt();
    httot+=et;
    if (std::abs(etmis.DeltaPhi(jets[ijet]))<TMath::Pi()/2){
      same_hemisphere_httot+=et;
    }
  }
  return same_hemisphere_httot/httot;
}

void CmsSingleLepton20Fb::Run(const Reader * treereader, const Reader * gentreereader, const double & weight){

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //CMS-SUS-13-011 p.3: "Events are required to have an electron (muon) with pt>30 (25) GeV. Electrons are required to lie in the barrel region of the ECAL
  //(|eta|<1.4442) while muons are considered up to |eta|<2.1."
  std::vector<jlepton> leptons=leptonSkim(treereader->GetElec(),treereader->GetMuon(), 30, 1.4442,25.0,2.1);
  std::vector<jlepton> allleptons=leptonSkim(treereader->GetElec(),treereader->GetMuon(), 5.0, 1.4442,5.0,2.1);
  //CMS-SUS-13-011 p.4: "Selected events are required to contain at least four jets with pt>30 GeV and |eta|<2.4."
  std::vector<jjet> jets=goodjetsSkim(treereader->GetJet(),30,2.4); 
  std::vector<jjet> bjets=goodbjetsSkim(treereader->GetJet(),30,2.4); 
  //CMS-SUS-13-011 p.4: "We also reject ... , as well as a jet of pt>20 GeV consitent with the hadronic decay of a tau lepton"
  std::vector<jjet> taujets=goodjetsSkim(treereader->GetTauJet(),20,2.4); 
  jjet etmis=(treereader->GetMet())[0]; 

  double met=etmis.Et();
  njets->Fill(jets.size());
  nbjets->Fill(bjets.size());
  //preselection
  bool preselected=false;
  // one isolated lepton
  cut_flow_hist->Fill(0.5);
  if (leptons.size()==1){
    cut_flow_hist->Fill(1.5);
    // no second isolated lepton of pt>5GeV. FIXME: isolation algorithm should be: pTsum<0.2*pT && DeltaR<0.4
    if (allleptons.size()==1){
      cut_flow_hist->Fill(2.5);
      // FIXME: no additional isolated track of pt>10 GeV with opposite sign wrt primary lepton
      if (true){
        cut_flow_hist->Fill(3.5);
        //no jet of pt>20 GeV "consistent with hadronic tau lepton" 
        if (taujets.size()==0){
          cut_flow_hist->Fill(4.5);
          // at least four jets, at least one btag
          if (jets.size()>=4 && bjets.size()>=1){
            cut_flow_hist->Fill(5.5);
            // MET>100 GeV
            if (met>100){
              cut_flow_hist->Fill(6.5);
              preselected=true;
            }
          }
        }
      }
    }
  }

  if (preselected){
    jlepton lepton=leptons[0];
    //get 8 kinematic variables
    double mt2w=get_mt2w(lepton,jets,etmis);
    double chi2=get_chi2(jets);
    double min_dphi=std::min(TMath::Abs(etmis.DeltaPhi(jets[0])),TMath::Abs(etmis.DeltaPhi(jets[1])));
    //FIXME: check if get_htratio has been defined correctly
    double htratio=get_htratio(jets,etmis);
    double leading_bjet_pt=bjets[0].Pt();
    double delta_R=lepton.DeltaR(bjets[0]);
    double lepton_pt=lepton.Pt();
    double mt=TMath::Sqrt(2*etmis.Et()*lepton.Pt()*(1-TMath::Cos(lepton.DeltaPhi(etmis))));
    //fill histograms (cf. Fig. 2)
    mt_hist->Fill(mt,1000*weight);
    met_hist->Fill(met,1000*weight);
    mt2w_hist->Fill(mt2w,1000*weight);
    chi2_hist->Fill(chi2,1000*weight);
    htratio_hist->Fill(htratio,1000*weight);
    min_dphi_hist->Fill(min_dphi,1000*weight);
    leading_bjet_pt_hist->Fill(leading_bjet_pt,10*weight);
    delta_R_hist->Fill(delta_R,1000*weight);
    lepton_pt_hist->Fill(lepton_pt,10*weight);
    //fill other histograms
    lepton_flavour->Fill((lepton.Flavour()=="muon")+0.5);
    //choose which signal regions to use
    enum SignalRegionsSelection {SRall,SRT2tt,SRT2bbww};
    SignalRegionsSelection signalRegions=SRT2tt;
    if (signalRegions==SRall){
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
    else if (signalRegions=SRT2tt){
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
    }
    else if (signalRegions=SRT2bbww){
      // stop->bot char; low DeltaM
      if (met>100 && min_dphi>0.8 ) mSigPred.at(0)+=weight;
      if (met>150 && min_dphi>0.8 ) mSigPred.at(1)+=weight;
      if (met>200 && min_dphi>0.8 ) mSigPred.at(2)+=weight;
      if (met>250 && min_dphi>0.8 ) mSigPred.at(3)+=weight;
      // stop->bot char; high DeltaM
      if (met>100 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(4)+=weight;
      if (met>150 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(5)+=weight;
      if (met>200 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(6)+=weight;
      if (met>250 && mt2w>200 && min_dphi>0.8 && leading_bjet_pt>100 ) mSigPred.at(7)+=weight;
    }
  }
  return;
}
