//Implementation of the 0 lepton search SUS13011
//arXiv:1305.2390 v2

#include "analysis_zerolepmt2_8_20.hh"

double lundDistance(jjet testJet, jjet refJet){
  double Ei=refJet.E();
  double Ek=testJet.E();
  double pi=refJet.P();
//  double eTot=testJet.E()+refJet.E();
  //NOTE: cos(theta)=cos(-theta)
  double dTheta = testJet.Theta()-refJet.Theta();
//  return (refJet.E()-refJet.P()*TMath::Cos(dTheta) )*refJet.E()/(eTot*eTot);
//  return (Ei-pi*TMath::Cos(dTheta))*Ei/((Ei+Ek)*(Ei+Ek));
  double brackets = testJet*refJet/Ek;
  return brackets*Ei/((Ei+Ek)*(Ei+Ek));
}

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
ZeroLepMt2::ZeroLepMt2(const std::string & name, 
    const std::string & experiment,
    const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

  //constructor for object where you also want to run a limit and so other information
  //is necessary
  ZeroLepMt2::ZeroLepMt2(const std::string & name, 
      const std::string & experiment, 
      const unsigned int & numBins,
		       const double & intlumi, 
		       //const std::vector<int> & datayields,
		       const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

ZeroLepMt2::ZeroLepMt2(const std::string & name, 
		       const std::string & experiment, 
           const unsigned int & numBins,
           const double & intlumi, 
           const std::vector<double> & bgpred,
           const std::vector<double> & bgpreduncert,
           const std::vector<int> & datayields,
           const std::string & fitmode,
           const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


  ZeroLepMt2::~ZeroLepMt2() {}


  void ZeroLepMt2::initHistos() {
    andir->cd();
//    leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
    ht20hist = new TH1D("ht20hist", ";H_{T} [GeV];Entries",250,-5.,2495.);
    ht40hist = new TH1D("ht40hist", ";H_{T} [GeV];Entries",250,-5.,2495.);
    ht50hist = new TH1D("ht50hist", ";H_{T} [GeV];Entries",250,-5.,2495.);
    mt2hist = new TH1D("mt2hist", ";M_{T2} [GeV];Entries",51, 0.0, 750.0);
//    mhthist = new TH1D("mhthist", ";Missing H_{T} [GeV];Entries",200,-5.,1995.);
    methist = new TH1D("methist", ";CALO Missing E_{T} [GeV];Entries",200,-5.,1995.);
    njets20hist = new TH1D("njets20", ";N_{jets};Entries",10,-0.5,9.5);
    njets40hist = new TH1D("njets40", ";N_{jets};Entries",10,-0.5,9.5);
    njets50hist = new TH1D("njets50", ";N_{jets};Entries",10,-0.5,9.5);
//    mht_over_ht = new TH1D("mht_over_ht",";MH_{T}/H_{T};Entries",200,-0.005, 1.995);
    btagrate20hist = new TH1D("btagrate20",";#b-tags;Entries",10,-0.5,9.5);
    btagrate40hist = new TH1D("btagrate40",";#b-tags;Entries",10,-0.5,9.5);
    btagrate50hist = new TH1D("btagrate50",";#b-tags;Entries",10,-0.5,9.5);
//    met_vs_mht = new TH2D("met_vs_mht",";MET;MHT",200,-5.,1995., 200,-5.,1995.);
    cut_flow = new TH1D("cut_flow",";;cuts",7,-0.5,6.5);
    low_ht_met_vs_mt2 = new TH2D("low_ht_met_vs_mt2",
            "low H_T;MET[GeV];M_{T2}[GeV]", 60, 0.0, 1500.0, 60, 0.0, 1500.0);
    medium_ht_met_vs_mt2 = new TH2D("medium_ht_met_vs_mt2",
            "medium H_T;MET[GeV];M_{T2}[GeV]", 60, 0.0, 1500.0, 60, 0.0, 1500.0);
    high_ht_met_vs_mt2 = new TH2D("high_ht_met_vs_mt2",
            "high H_T;MET[GeV];M_{T2}[GeV]", 60, 0.0, 1500.0, 60, 0.0, 1500.0);
  }

void ZeroLepMt2::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {
  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over
  cut_flow->Fill(0);

  //Get the required variables and fill the histograms

  std::vector<jjet> etmis = treereader->GetMet(); //Missing transverse energy array
  double met = (etmis.size() == 1) ? etmis[0].E(): -1;
  if(met < -0.5) std::cout << "MET error in 0 lepton search" << std::endl; 
  methist->Fill(met, weight);

  std::vector<jlepton> elecs = treereader->GetElec();
  std::vector<jlepton> muons = treereader->GetMuon();
  std::vector<jlepton> leptons = leptonSkim(elecs, muons, 10.0, 2.4, 10.0, 2.4);
  std::vector<jlepton> leptons20 = leptonSkim(elecs, muons, 20.0, 2.4, 20.0, 2.4);
  std::vector<jjet> taus = goodjetsSkim(treereader->GetTauJet(), 20.0, 2.3);
  //AN2013_215_v10 p.17 "L1FastJet and L2L3 corrected pf-CHS (charged hadron 
  //subracted)-jets with pT > 40 GeV and |η | < 2.4.
  std::vector<jjet> goodjets40 = goodjetsSkim(treereader->GetJet(), 40.0, 2.4);
  //For the HT calculation pf-CHS-jets with pT > 50 GeV and |η | < 3.0 are used.
  // This choice is done such that the HT trigger is fully eﬃcient at an oﬄine 
  // HT > 750 GeV.
  std::vector<jjet> goodjets50 = goodjetsSkim(treereader->GetJet(), 50.0, 3.0);
  //For the HTmiss and hemisphere calculation we use pf-CHS-jets with pT > 20 
  //GeV and |η | < 2.4, that pass the loose jet ID. With this choice we have the
  //smallest bias of the MT2 distribution in QCD events.
  std::vector<jjet> goodjets20 = goodjetsSkim(treereader->GetJet(), 20.0, 2.4);

  //HT used in the analysis
  double ht = getht(goodjets50);
  //all values of HT for sanity check
  double ht20 = getht(goodjets20);
  double ht40 = getht(goodjets40);
  double ht50 = getht(goodjets50);
  ht20hist->Fill(ht20, weight);
  ht40hist->Fill(ht40, weight);
  ht50hist->Fill(ht50, weight);
  //#b-tags used in the analysis
  int nBTags = getnbtags(goodjets40);
  //all b-tag values for sanity check
  int nBTags20 = getnbtags(goodjets20);
  int nBTags40 = getnbtags(goodjets40);
  int nBTags50 = getnbtags(goodjets50);
  btagrate20hist->Fill(nBTags20, weight);
  btagrate40hist->Fill(nBTags40, weight);
  btagrate50hist->Fill(nBTags50, weight);
  //# jets
  int nJets = goodjets40.size();
  // all #jets for sanity check
  int nJets20 = goodjets20.size();
  int nJets40 = goodjets40.size();
  int nJets50 = goodjets50.size();
  njets20hist->Fill(nJets20, weight);
  njets40hist->Fill(nJets40, weight);
  njets50hist->Fill(nJets50, weight);

  //If the general conditions are passed
  bool selected=false;
  // CMS PAS SUS-13-019 p4. and p5.
  // "At least two good jets are required ,,, "
  if(nJets >= 2){
      cut_flow->Fill(1);
      // "... the two leading jets with pt>100 GeV."
      if (goodjets40[0].Pt() > 100.0 && goodjets40[1].Pt() > 100.0){
          cut_flow->Fill(2);
          // "Events with possible contributions from beam halo processes, 
          // anomalous calorimeter or traking noise are rejected. In order to 
          // reject events with an important contribution of soft and/or forward
          // jets to the momentum imbalance, a maximum difference of 70 GeV is 
          // imposed between the MET vector and the vector sum of the pt of all
          // leptons and jets candidates with pt>20 passing the jet ID."
          TLorentzVector vectorSum = TLorentzVector();
          for(std::vector<jjet>::const_iterator j=goodjets20.begin(); 
                  j!=goodjets20.end(); j++) 
              vectorSum+=(TLorentzVector) *j;
          for(std::vector<jjet>::const_iterator t=taus.begin(); 
                  t!=taus.end(); t++) 
              vectorSum+=(TLorentzVector) *t;
          for(std::vector<jlepton>::const_iterator l=leptons.begin(); 
                  l!=leptons.end(); l++) 
              vectorSum += (TLorentzVector) *l;
          if(fabs(met-vectorSum.Pt()) < 70.0){
              cut_flow->Fill(3);
              //"It has been argued previously that MT2 is protected against jet
              // energy mismeasurement in dijet events. In multijet evens such 
              // mismeasurement can lead to hemispheres not being back-to-back 
              // and resultsing in larger values of MT2. To protect against this
              // effect a minimum difference in azimuth phi between any of the 
              // four leading jets and the MET, deltaphi >=0.3, is required."
              bool dPhiRequirement = true;
              int iMax = (nJets<4) ? nJets : 4;
              for(unsigned i=0; i<iMax; ++i){
                double dPhi=TMath::Abs(etmis[0].DeltaPhi(goodjets40[i]));
                dPhiRequirement = (dPhi>=0.3);
                if(!dPhiRequirement) break;
              }
              if(dPhiRequirement){ 
                  cut_flow->Fill(4);
                  //"Finally, events are vetoed if they contain and isolated 
                  //electron, muon or tau, to suppress the contributions from 
                  //W(l nu)+jets, Z(ll)+jets, and top backgrounds."
                  if(leptons.size() == 0 && taus.size() == 0){ 
                      cut_flow->Fill(5);
                      selected=true;
                  }
              }
          } 
      }
  }
  if (selected){

    // to calculate MT2 we first need to construct two pseudo jets
    // (see description in include/zerolepmt2_functions.hh)
    std::vector<int> pseudoJetsGrouping=getPseudoJetsGrouping(goodjets20);
    jjet jet1,jet2;
    for (int i=0;i<pseudoJetsGrouping.size();i++){
        if (pseudoJetsGrouping[i]==1)
            jet1+=goodjets20[i];
        else
            jet2+=goodjets20[i];
    }

    //Calculate mt2
    mt2_bisect::mt2 mt2_event;

    // CMS PAS SUS-13-019 p5. "Massless pseudo-jets have been used as input to 
    // MT2 and zero test mass, as this was found..."
    double pa[3],pb[3],pmiss[3];
    pa[0]=0.;
    pb[0]=0.;
    pmiss[0]=0.;
    pa[1]=jet1.Px();
    pb[1]=jet2.Px();
    pmiss[1]=etmis[0].Px();
    pa[2]=jet1.Py();
    pb[2]=jet2.Py();
    pmiss[2]=etmis[0].Py();

    mt2_event.set_momenta( pa, pb, pmiss );
    // zero test mass
    mt2_event.set_mn( 0. );

    double Mt2 = mt2_event.get_mt2();       
    mt2hist->Fill(Mt2,weight);

    //Low Ht region first
    if (450 < ht && ht < 750)
        low_ht_met_vs_mt2->Fill(met, Mt2, weight);
    if (750 < ht && ht < 1200)
        medium_ht_met_vs_mt2->Fill(met, Mt2, weight);
    if (1200 < ht )
        high_ht_met_vs_mt2->Fill(met, Mt2, weight);
    if(ht>450. && ht<750. && met>=200.){

      if(nJets==2 && nBTags==0){

        if(Mt2>=200. && Mt2 < 240.) mSigPred.at(0)+=weight;
        else if(Mt2>=240. && Mt2 < 290.) mSigPred.at(1)+=weight;
        else if(Mt2>=290. && Mt2 < 350.) mSigPred.at(2)+=weight;
        else if(Mt2>=350. && Mt2 < 420.) mSigPred.at(3)+=weight;
        else if(Mt2>=420. && Mt2 < 490.) mSigPred.at(4)+=weight;
        else if(Mt2>=490. && Mt2 < 570.) mSigPred.at(5)+=weight;
        else if(Mt2>=570. && Mt2 < 650.) mSigPred.at(6)+=weight;
        else if(Mt2>=650.) mSigPred.at(7)+=weight;

      }else if(nJets==2 && nBTags>=1){

        if(Mt2>=200. && Mt2 < 250.) mSigPred.at(8)+=weight;
        else if(Mt2>=250. && Mt2 < 310.) mSigPred.at(9)+=weight;
        else if(Mt2>=310. && Mt2 < 380.) mSigPred.at(10)+=weight;
        else if(Mt2>=380. && Mt2 < 450.) mSigPred.at(11)+=weight;
        else if(Mt2>=450. && Mt2 < 550.) mSigPred.at(12)+=weight;
        else if(Mt2>=550.) mSigPred.at(13)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==0){

        if(Mt2>=200. && Mt2 < 240.) mSigPred.at(14)+=weight;
        else if(Mt2>=240. && Mt2 < 290.) mSigPred.at(15)+=weight;
        else if(Mt2>=290. && Mt2 < 350.) mSigPred.at(16)+=weight;
        else if(Mt2>=350. && Mt2 < 420.) mSigPred.at(17)+=weight;
        else if(Mt2>=420. && Mt2 < 490.) mSigPred.at(18)+=weight;
        else if(Mt2>=490. && Mt2 < 570.) mSigPred.at(19)+=weight;
        else if(Mt2>=570. && Mt2 < 650.) mSigPred.at(20)+=weight;
        else if(Mt2>=650.) mSigPred.at(21)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==1){

        if(Mt2>=200. && Mt2 < 250.) mSigPred.at(22)+=weight;
        else if(Mt2>=250. && Mt2 < 310.) mSigPred.at(23)+=weight;
        else if(Mt2>=310. && Mt2 < 380.) mSigPred.at(24)+=weight;
        else if(Mt2>=380. && Mt2 < 450.) mSigPred.at(25)+=weight;
        else if(Mt2>=450. && Mt2 < 550.) mSigPred.at(26)+=weight;
        else if(Mt2>=550.) mSigPred.at(27)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==2){

        if(Mt2>=200. && Mt2 < 250.) mSigPred.at(28)+=weight;
        else if(Mt2>=250. && Mt2 < 325.) mSigPred.at(29)+=weight;
        else if(Mt2>=325. && Mt2 < 425.) mSigPred.at(30)+=weight;
        else if(Mt2>=425.) mSigPred.at(31)+=weight;

      }else if(nJets>=6 && nBTags==0){

        if(Mt2>=200. && Mt2 < 280.) mSigPred.at(32)+=weight;
        else if(Mt2>=280. && Mt2 < 380.) mSigPred.at(33)+=weight;
        else if(Mt2>=380.) mSigPred.at(34)+=weight;

      }else if(nJets>=6 && nBTags==1){

        if(Mt2>=200. && Mt2 < 250.) mSigPred.at(35)+=weight;
        else if(Mt2>=250. && Mt2 < 325.) mSigPred.at(36)+=weight;
        else if(Mt2>=325.) mSigPred.at(37)+=weight;

      }else if(nJets>=6 && nBTags==2){

        if(Mt2>=200. && Mt2 < 250.) mSigPred.at(38)+=weight;
        else if(Mt2>=250. && Mt2 < 300.) mSigPred.at(39)+=weight;
        else if(Mt2>=300.) mSigPred.at(40)+=weight;

      }else if(nJets>=3 && nBTags>=3){

        if(Mt2>=200. && Mt2 < 280.) mSigPred.at(41)+=weight;
        else if(Mt2>=280.) mSigPred.at(42)+=weight;

      }

    }else if(ht>=750. && ht<1200.){

      if(nJets==2 && nBTags==0){

        if(Mt2>=125. && Mt2 < 150.) mSigPred.at(43)+=weight;
        else if(Mt2>=150. && Mt2 < 180.) mSigPred.at(44)+=weight;
        else if(Mt2>=180. && Mt2 < 220.) mSigPred.at(45)+=weight;
        else if(Mt2>=220. && Mt2 < 270.) mSigPred.at(46)+=weight;
        else if(Mt2>=270. && Mt2 < 325.) mSigPred.at(47)+=weight;
        else if(Mt2>=325. && Mt2 < 425.) mSigPred.at(48)+=weight;
        else if(Mt2>=425. && Mt2 < 580.) mSigPred.at(49)+=weight;
        else if(Mt2>=580. && Mt2 < 780.) mSigPred.at(50)+=weight;
        else if(Mt2>=780.) mSigPred.at(51)+=weight;

      }else if(nJets==2 && nBTags>=1){

        if(Mt2>=100. && Mt2 < 135.) mSigPred.at(52)+=weight;
        else if(Mt2>=135. && Mt2 < 170.) mSigPred.at(53)+=weight;
        else if(Mt2>=170. && Mt2 < 260.) mSigPred.at(54)+=weight;
        else if(Mt2>=260. && Mt2 < 450.) mSigPred.at(55)+=weight;
        else if(Mt2>=450.) mSigPred.at(56)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==0){

        if(Mt2>=160. && Mt2 < 185.) mSigPred.at(57)+=weight;
        else if(Mt2>=185. && Mt2 < 215.) mSigPred.at(58)+=weight;
        else if(Mt2>=215. && Mt2 < 250.) mSigPred.at(59)+=weight;
        else if(Mt2>=250. && Mt2 < 300.) mSigPred.at(60)+=weight;
        else if(Mt2>=300. && Mt2 < 370.) mSigPred.at(61)+=weight;
        else if(Mt2>=370. && Mt2 < 480.) mSigPred.at(62)+=weight;
        else if(Mt2>=480. && Mt2 < 640.) mSigPred.at(63)+=weight;
        else if(Mt2>=640. && Mt2 < 800.) mSigPred.at(64)+=weight;
        else if(Mt2>=800.) mSigPred.at(65)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==1){

        if(Mt2>=150. && Mt2 < 175.) mSigPred.at(66)+=weight;
        else if(Mt2>=175. && Mt2 < 210.) mSigPred.at(67)+=weight;
        else if(Mt2>=210. && Mt2 < 270.) mSigPred.at(68)+=weight;
        else if(Mt2>=270. && Mt2 < 380.) mSigPred.at(69)+=weight;
        else if(Mt2>=380. && Mt2 < 600.) mSigPred.at(70)+=weight;
        else if(Mt2>=600.) mSigPred.at(71)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==2){

        if(Mt2>=130. && Mt2 < 160.) mSigPred.at(72)+=weight;
        else if(Mt2>=160. && Mt2 < 200.) mSigPred.at(73)+=weight;
        else if(Mt2>=200. && Mt2 < 270.) mSigPred.at(74)+=weight;
        else if(Mt2>=270. && Mt2 < 370.) mSigPred.at(75)+=weight;
        else if(Mt2>=370.) mSigPred.at(76)+=weight;

      }else if(nJets>=6 && nBTags==0){

        if(Mt2>=160. && Mt2 < 200.) mSigPred.at(77)+=weight;
        else if(Mt2>=200. && Mt2 < 250.) mSigPred.at(78)+=weight;
        else if(Mt2>=250. && Mt2 < 325.) mSigPred.at(79)+=weight;
        else if(Mt2>=325. && Mt2 < 425.) mSigPred.at(80)+=weight;
        else if(Mt2>=425.) mSigPred.at(81)+=weight;

      }else if(nJets>=6 && nBTags==1){

        if(Mt2>=150. && Mt2 < 190.) mSigPred.at(82)+=weight;
        else if(Mt2>=190. && Mt2 < 250.) mSigPred.at(83)+=weight;
        else if(Mt2>=250. && Mt2 < 350.) mSigPred.at(84)+=weight;
        else if(Mt2>=350.) mSigPred.at(85)+=weight;

      }else if(nJets>=6 && nBTags==2){

        if(Mt2>=130. && Mt2 < 170.) mSigPred.at(86)+=weight;
        else if(Mt2>=170. && Mt2 < 220.) mSigPred.at(87)+=weight;
        else if(Mt2>=220. && Mt2 < 300.) mSigPred.at(88)+=weight;
        else if(Mt2>=300.) mSigPred.at(89)+=weight;

      }else if(nJets>=3 && nBTags>=3){

        if(Mt2>=125. && Mt2 < 175.) mSigPred.at(90)+=weight;
        else if(Mt2>=175. && Mt2 < 275.) mSigPred.at(91)+=weight;
        else if(Mt2>=275.) mSigPred.at(92)+=weight;

      }

    }else if(ht>=1200.){

      if(nJets==2 && nBTags==0){

        if(Mt2>=120. && Mt2 < 150.) mSigPred.at(93)+=weight;
        else if(Mt2>=150. && Mt2 < 200.) mSigPred.at(94)+=weight;
        else if(Mt2>=200. && Mt2 < 260.) mSigPred.at(95)+=weight;
        else if(Mt2>=260. && Mt2 < 350.) mSigPred.at(96)+=weight;
        else if(Mt2>=350. && Mt2 < 550.) mSigPred.at(97)+=weight;
        else if(Mt2>=550.) mSigPred.at(98)+=weight;

      }else if(nJets==2 && nBTags>=1){

        if(Mt2>=100. && Mt2 < 180.) mSigPred.at(99)+=weight;
        else if(Mt2>=180.) mSigPred.at(100)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==0){

        if(Mt2>=160. && Mt2 < 185.) mSigPred.at(101)+=weight;
        else if(Mt2>=185. && Mt2 < 220.) mSigPred.at(102)+=weight;
        else if(Mt2>=220. && Mt2 < 270.) mSigPred.at(103)+=weight;
        else if(Mt2>=270. && Mt2 < 350.) mSigPred.at(104)+=weight;
        else if(Mt2>=350. && Mt2 < 450.) mSigPred.at(105)+=weight;
        else if(Mt2>=450. && Mt2 < 650.) mSigPred.at(106)+=weight;
        else if(Mt2>=650.) mSigPred.at(107)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==1){

        if(Mt2>=150. && Mt2 < 180.) mSigPred.at(108)+=weight;
        else if(Mt2>=180. && Mt2 < 230.) mSigPred.at(109)+=weight;
        else if(Mt2>=230. && Mt2 < 350.) mSigPred.at(110)+=weight;
        else if(Mt2>=350.) mSigPred.at(111)+=weight;

      }else if(nJets>=3 && nJets<=5 && nBTags==2){

        if(Mt2>=130. && Mt2 < 200.) mSigPred.at(112)+=weight;
        else if(Mt2>=200.) mSigPred.at(113)+=weight;

      }else if(nJets>=6 && nBTags==0){

        if(Mt2>=160. && Mt2 < 200.) mSigPred.at(114)+=weight;
        else if(Mt2>=200. && Mt2 < 300.) mSigPred.at(115)+=weight;
        else if(Mt2>=300.) mSigPred.at(116)+=weight;

      }else if(nJets>=6 && nBTags==1){

        if(Mt2>=150. && Mt2 < 200.) mSigPred.at(117)+=weight;
        else if(Mt2>=200. && Mt2 < 300.) mSigPred.at(118)+=weight;
        else if(Mt2>=300.) mSigPred.at(119)+=weight;

      }else if(nJets>=6 && nBTags==2){

        if(Mt2>=120. && Mt2 < 200.) mSigPred.at(120)+=weight;
        else if(Mt2>=200.) mSigPred.at(121)+=weight;

      }else if(nJets>=3 && nBTags>=3){

        if(Mt2>=125.) mSigPred.at(122)+=weight;

      }

    }


  } 
  return;
}
