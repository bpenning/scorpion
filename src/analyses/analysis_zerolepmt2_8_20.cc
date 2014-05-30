//Implementation of the 0 lepton search SUS13011
//arXiv:1305.2390 v2

#include "analysis_zerolepmt2_8_20.hh"

double lundDistance(jjet testJet, jjet refJet){

  double eTot=testJet.E()+refJet.E();
  double dTheta = fabs(testJet.Theta()-refJet.Theta());
  return (testJet.E()-testJet.P()*TMath::Cos(dTheta) )*testJet.E()/(eTot*eTot);

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
    leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
    hthist = new TH1D("hthist", ";H_{T} [GeV];Entries",250,-5.,2495.);
    mhthist = new TH1D("mhthist", ";Missing H_{T} [GeV];Entries",200,-5.,1995.);
    calomethist = new TH1D("calomethist", ";CALO Missing E_{T} [GeV];Entries",200,-5.,1995.);
    athist = new TH1D("athist", ";#alpha_{T};Normalised",200,-0.005,1.995);
    njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
    mht_over_ht = new TH1D("mht_over_ht",";MH_{T}/H_{T};Entries",200,-0.005, 1.995);
    btagrate = new TH1D("btagrate",";#b-tags;Entries",10,-0.5,9.5);
    calomet_vs_mht = new TH2D("calomet_vs_mht",";caloMET;MHT",200,-5.,1995., 200,-5.,1995.);
  }

void ZeroLepMt2::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //Get the required variables and fill the histograms

  std::vector<jjet> etmis = treereader->GetMet(); //Missing transverse energy array
  double MET = (etmis.size() == 1) ? etmis[0].E(): -1;
  if(MET < -0.5) std::cout << "MET error in 0 lepton search" << std::endl; 
  calomethist->Fill(MET, weight);

  std::vector<jlepton> elecs = treereader->GetElec();
  std::vector<jlepton> muons = treereader->GetMuon();

  std::vector<jlepton> leptons = leptonSkim(elecs, muons, 10.0, 2.4, 10.0, 2.4, 1.4442, 1.566);
  std::vector<jlepton> leptons20 = leptonSkim(elecs, muons, 20.0, 2.4, 20.0, 2.4, 1.4442, 1.566);

  std::vector<jjet> taus = goodjetsSkim(treereader->GetTauJet(), 20.0, 2.3);

  std::vector<jjet> goodjets = goodjetsSkim(treereader->GetJet(), 40.0, 2.4);
  std::vector<jjet> jets20 = goodjetsSkim(treereader->GetJet(), 20.0, 2.4);
  double HT = getht(goodjets);

  double nBTags = getnbtags(goodjets);

  hthist->Fill(HT, weight);
  btagrate->Fill(nBTags, weight);


  //If the general conditions are passed
  if(goodjets.size() >= 2 && leptons.size() == 0 &&
      taus.size() == 0 && 
      goodjets[0].Pt() > 100.0 && goodjets[1].Pt() > 100.0
    ) {

    //Check the min difference in dPhi between 4 leading jets and MET
    bool dPhiRequirement = true;
    int iMax = (goodjets.size()<4) ? goodjets.size() : 4;

    for(unsigned i=0; i<iMax; ++i){
      dPhiRequirement = (etmis[0].DeltaPhi(goodjets[i]) >= 0.3);
      if(!dPhiRequirement) break;
    }

    if(!dPhiRequirement) return;

    //Check the difference between MET and vector sum of PT of all leptons and jets

    TLorentzVector vectorSum = TLorentzVector();

    for(std::vector<jjet>::const_iterator j=goodjets.begin(); j!=goodjets.end(); j++) vectorSum+=(TLorentzVector) *j;
    for(std::vector<jjet>::const_iterator t=taus.begin(); t!=taus.end(); t++) vectorSum+=(TLorentzVector) *t;
    for(std::vector<jlepton>::const_iterator l=leptons.begin(); l!=leptons.end(); l++) vectorSum += (TLorentzVector) *l;

    if(fabs(MET-vectorSum.Pt()) > 70.0) return;

    //Calculate MT2, in this case we use massless pseudojets and a massless MET
    //
    //First we must find 2 pseudojets using the hemisphere association method, 
    //seeded with the direction of the two jets with largest dijet invariant mass
    //
    //Method described in http://iopscience.iop.org/0954-3899/34/6/S01/pdf/0954-3899_34_6_S01.pdf

    //Find the two jets that make the largest invariant mass

    double maxInvMass = 0.0;
    int jet1Index, jet2Index;

    std::vector<jjet> zeroMgoodjets = goodjets;

    for(std::vector<jjet>::iterator j=zeroMgoodjets.begin(); j!=zeroMgoodjets.end(); j++) j->setZeroMass();

    for(unsigned j1=0; j1<zeroMgoodjets.size(); j1++){
      for(unsigned j2=j1+1; j2<zeroMgoodjets.size(); j2++){

        double invMass = (zeroMgoodjets[j1]+zeroMgoodjets[j2]).M();
        if(invMass > maxInvMass){
          maxInvMass = invMass;
          jet1Index = j1;
          jet2Index = j2;
        }

      }
    }

    //pseudoJet1 and jet2 are the seed jets
    //find which hemisphere each of the other jets is associated with
    
    //Vectors to store the jets to be associated in a hemisphere
    std::vector<jjet> hemisphere1, hemisphere2;
    hemisphere1.push_back(goodjets[jet1Index]);
    hemisphere2.push_back(goodjets[jet2Index]);

    for(unsigned i=0; i<goodjets.size(); i++){
      if(i==jet1Index || i==jet2Index) continue;
      if(lundDistance(goodjets[i],goodjets[jet1Index]) <= lundDistance(goodjets[i],goodjets[jet2Index])) hemisphere1.push_back(goodjets[i]);
      else hemisphere2.push_back(goodjets[i]);
    }


    //Keep track of which hemisphere a jet was in
    std::vector<int> inHemisphere(goodjets.size(),0);

    //Iterate until no jets change hemisphere
    while(1){

      bool changeHemisphere=false;

      //Recalculate the axes and iterate
      jjet axis1,axis2;

      for(std::vector<jjet>::const_iterator it=hemisphere1.begin(); it!=hemisphere1.end();it++) axis1+=*it;
      for(std::vector<jjet>::const_iterator it=hemisphere2.begin(); it!=hemisphere2.end();it++) axis2+=*it;

      //Set the masses to 0
      axis1.setZeroMass();
      axis2.setZeroMass();

      //Reset the hemispheres
      hemisphere1.clear();
      hemisphere2.clear();

      //Associate the jets to the new hemispheres
      for(unsigned i=0; i<goodjets.size(); i++){
        if(lundDistance(goodjets[i],axis1) <= lundDistance(goodjets[i],axis2)){
          hemisphere1.push_back(goodjets[i]);
          if(inHemisphere[i]!=1) changeHemisphere=true;
          inHemisphere[i]=1;
        }else{
          hemisphere2.push_back(goodjets[i]);
          if(inHemisphere[i]!=2) changeHemisphere=true;
          inHemisphere[i]=2;
        }
      }

      if(!changeHemisphere) break;
    }

    jjet jet1,jet2;
    for(std::vector<jjet>::const_iterator it=hemisphere1.begin(); it!=hemisphere1.end();it++) jet1+=*it;
    for(std::vector<jjet>::const_iterator it=hemisphere2.begin(); it!=hemisphere2.end();it++) jet2+=*it;

    //Set the masses to 0
    jet1.setZeroMass();
    jet2.setZeroMass();

    //Calculate mt2
    mt2_bisect::mt2 mt2_event;

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
    mt2_event.set_mn( 0. );

    double Mt2 = mt2_event.get_mt2();       

    int nJets = goodjets.size();
    //Low Ht region first
    if(HT>450. && HT<750. && MET>=200.){

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

    }else if(HT>=750. && HT<1200.){

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

    }else if(HT>=1200.){

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
