#include "analysis_ss8high.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
SS8high::SS8high(const std::string & name, 
	   const std::string & experiment,
	   const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
SS8high::SS8high(const std::string & name, 
	   const std::string & experiment, 
	   const unsigned int & numBins,
	   const double & intlumi, 
	   //const std::vector<int> & datayields,
	   const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

SS8high::SS8high(const std::string & name, 
	   const std::string & experiment, 
	   const unsigned int & numBins,
	   const double & intlumi, 
	   const std::vector<double> & bgpred,
	   const std::vector<double> & bgpreduncert,
	   const std::vector<int> & datayields,
	   const std::string & fitmode,
	   const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


SS8high::~SS8high() {}

void SS8high::initHistos() {
  andir->cd();
  ptleadinglep = new TH1D("ptleadinglep",";leading p_{T} [GeV];",100, -5.0, 995.0);
  ptchargeonevschargetwo = new TH2D("ptchargeonevschargetwo",";charge leading lep;charge sub-leading lep",3, -1.5, 1.5, 3, -1.5, 1.5);
  njethist = new TH1D("njethist",";N_{Jets};",20, -0.5, 19.5);
  hthist = new TH1D("hthist",";H_{T} [GeV];",250, -5.0, 2495.0);
  methist = new TH1D("methist",";E_{T}^{miss} [GeV];", 100, -5.0, 995.0);
  leadingjetpt = new TH1D("leadingjetpt",";leading-jet p_{T} [GeV];",100, -5.0, 995.0);
  leadingbjetpt = new TH1D("leadingbjetpt",";leading-bjet p_{T} [GeV];",100, -5.0, 995.0);
  numbjets = new TH1D("numbjets",";#b-tags;", 10, -0.5, 9.5);
  sfinvmass = new TH1D("sfinvmass",";M_{SF} [GeV];", 1000, -0.5, 999.5);
  drjetcombo = new TH1D("drjetcombo",";#DeltaR;",500, -0.005, 4.995);
  double leptonEfficiencyBinEdges[]={10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 65.0, 80.0, 100.0};
  genelectronshist = new TH1D("genelectronshist",";gen electron $p_T$ (GeV);Efficiency",9,leptonEfficiencyBinEdges);
  recelectronshist = new TH1D("recelectronshist",";gen electron $p_T$ (GeV);Efficiency",9,leptonEfficiencyBinEdges);
  genmuonshist = new TH1D("genmuonshist",";gen muon $p_T$ (GeV);Efficiency",9,leptonEfficiencyBinEdges);
  recmuonshist = new TH1D("recmuonshist",";gen muon $p_T$ (GeV);Efficiency",9,leptonEfficiencyBinEdges);
  gentaushist = new TH1D("gentaushist",";gen tau $p_T$ (GeV);Efficiency", 10, 0.0, 100.0);
  rectaushist = new TH1D("rectaushist",";gen tau $p_T$ (GeV);Efficiency", 10, 0.0, 100.0);
  genbjetshist = new TH1D("genbjetshist",";gen bjet $p_T$ (GeV);Efficiency",56,40.0,600.0);
  recbjetshist = new TH1D("recbjetshist",";gen bjet $p_T$ (GeV);Efficiency",56,40.0,600.0);
}

void SS8high::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> elecs = treereader->GetElec();
  std::vector<jlepton> muons = treereader->GetMuon();
  
  std::vector<jlepton> leptons = leptonSkim(elecs, muons, 5.0, 2.4, 5.0, 2.4, 1.4442, 1.566);
  std::vector<jjet> goodjets = goodjetsSkim(treereader->GetJet(), 40.0, 2.4);
  std::vector<jjet> goodbjets = goodbjetsSkim(treereader->GetJet(), 40.0, 2.4);

  //Cross check for b-tag efficiency and muon/electron selection efficiencies
  std::vector<jlepton> recElectrons = goodleptons(elecs,10.0, 2.4, 1.4442, 1.566);
  std::vector<jparticle> genElectrons = getGenParticles(
          gentreereader->GetGenParticle(), 11, 3, 10.0, 2.4);
  std::vector<jlepton> recMuons = goodleptons(muons,10.0, 2.4);
  std::vector<jparticle> genMuons = getGenParticles(
          gentreereader->GetGenParticle(), 13, 3, 10.0, 2.4);
  std::vector<jjet> recTaus = goodjetsSkim(treereader->GetTauJet(),10.0, 2.3);
  std::vector<jparticle> genTaus = getGenParticles(
          gentreereader->GetGenParticle(), 15, 3, 10.0, 2.3);
  std::vector<jjet> recBjets = goodbjetsSkim(treereader->GetJet(), 40.0, 2.4);
  std::vector<jparticle> genBjets = getGenParticles(
          gentreereader->GetGenParticle(), 5, 3, 40.0, 2.4);

  std::vector<jparticle>::const_iterator genElectron;
  std::vector<jlepton>::const_iterator recElectron;
  for (genElectron=genElectrons.begin();genElectron!=genElectrons.end();genElectron++){
    double minDeltaR=1e9;
    for (recElectron=recElectrons.begin();recElectron!=recElectrons.end();recElectron++){
      double deltaR=recElectron->DeltaR(*genElectron);
      if (deltaR<minDeltaR)
        minDeltaR=deltaR;
    }
    genelectronshist->Fill(genElectron->Pt());
    if (minDeltaR<0.5)
      recelectronshist->Fill(genElectron->Pt());
  }

  std::vector<jparticle>::const_iterator genMuon;
  std::vector<jlepton>::const_iterator recMuon;
  for (genMuon=genMuons.begin();genMuon!=genMuons.end();genMuon++){
    double minDeltaR=1e9;
    for (recMuon=recMuons.begin();recMuon!=recMuons.end();recMuon++){
      double deltaR=recMuon->DeltaR(*genMuon);
      if (deltaR<minDeltaR)
        minDeltaR=deltaR;
    }
    genmuonshist->Fill(genMuon->Pt());
    if (minDeltaR<0.5)
      recmuonshist->Fill(genMuon->Pt());
  }

  std::vector<jparticle>::const_iterator genTau;
  std::vector<jjet>::const_iterator recTau;
  for (genTau=genTaus.begin();genTau!=genTaus.end();genTau++){
    double minDeltaR=1e9;
    for (recTau=recTaus.begin();recTau!=recTaus.end();recTau++){
      double deltaR=recTau->DeltaR(*genTau);
      if (deltaR<minDeltaR)
        minDeltaR=deltaR;
    }
    gentaushist->Fill(genTau->Pt());
    if (minDeltaR<0.5)
      rectaushist->Fill(genTau->Pt());
  }

  std::vector<jparticle>::const_iterator genBjet;
  std::vector<jjet>::const_iterator recBjet;
  for (genBjet=genBjets.begin();genBjet!=genBjets.end();genBjet++){
    double minDeltaR=1e9;
    for (recBjet=recBjets.begin();recBjet!=recBjets.end();recBjet++){
      double deltaR=recBjet->DeltaR(*genBjet);
      if (deltaR<minDeltaR)
        minDeltaR=deltaR;
    }
    genbjetshist->Fill(genBjet->Pt());
    if (minDeltaR<0.5)
      recbjetshist->Fill(genBjet->Pt());
  }

  //plot dr of each jet with every other jet
  for(unsigned int i=0; i<goodjets.size(); i++) {
    for(unsigned int j=0; j<goodjets.size(); j++) {
      if(j > i) { //prevent jet matching with itself
        drjetcombo->Fill(goodjets[i].DeltaR(goodjets[j]), weight);
      }
    }
  }

  unsigned int numjets = goodjets.size();


  if(numjets >= 2) {
    double HT = 1.05 * getht(goodjets); //5% scale correction

    std::vector<jjet> etmis = treereader->GetMet(); //Missing transverse energy array
    double MET = (etmis.size() == 1) ? etmis[0].E(): -1;

    njethist->Fill(numjets, weight);
    leadingjetpt->Fill(goodjets[0].Pt(), weight);

    //check for two or three leptons that pass selection - otherwise we ignore these events
    if(leptons.size() == 2 || leptons.size() == 3) { 
      ptleadinglep->Fill(leptons[0].Pt(), weight);
      ptchargeonevschargetwo->Fill(leptons[0].Charge(), leptons[1].Charge(), weight);
      if(leptons[0].Pt() > 20.0 && leptons[1].Pt() > 20.0) { //check both leading leptons have at least 20 GeV
        if(leptons[0].Charge() == leptons[1].Charge()) { //check first two leptons have the same charge
          double invmass = (leptons[0] + leptons[1]).M();
          sfinvmass->Fill(invmass, weight);
          if(invmass>8.0) { //if invmass > 8, we are ok
            bool passselection = true;
            if(leptons.size() == 3) { //if there are 3 leptons, need to check for Z-hypothesis if opposite charge
              if(leptons[2].Charge() != leptons[0].Charge()) { //check for opposite charge i.e. maybe consistent with Z
                //check for Z-veto
                double invmass1 = (leptons[0] + leptons[2]).M();
                bool oppflav1 = false;
                if(leptons[0].Flavour() != leptons[2].Flavour()) {
                  oppflav1 = true;
                }
                double invmass2 = (leptons[1] + leptons[2]).M();
                bool oppflav2 = false;
                if(leptons[1].Flavour() != leptons[2].Flavour()) {
                  oppflav2 = true;
                }
                if((invmass1 < 106.0 && invmass1 > 76.0 && !oppflav1 && leptons[2].Pt()>10.0) || 
                    (invmass2 < 106.0 && invmass2 > 76.0 && !oppflav2 && leptons[2].Pt()>10.0) ||
                    (invmass1 < 12.0  &&                    !oppflav1 && leptons[2].Pt()>5.0 ) ||             
                    (invmass2 < 12.0  &&                    !oppflav2 && leptons[2].Pt()>5.0 )              
                  ) {
                  passselection = false; //i.e. we have a Z candidate most likely, veto event...
                }
              }
            }
            if(passselection) {
              //we have a same-sign pair
              hthist->Fill(HT, weight);
              methist->Fill(MET, weight);
              unsigned int numbtags = getnbtags(goodjets); 
              numbjets->Fill(numbtags, weight);
              if (numbtags>=1)
                leadingbjetpt->Fill(goodbjets[0].Pt(), weight);

              //the simplified use a selection of the signal models
              enum SignalRegionsSelection {SRall,SRbTag2,SRbTag12};
              SignalRegionsSelection signalRegions=SRall;
              if (signalRegions==SRall){

                // Add events to signal regions

                //B tag 0 signal regions
                if(numbtags == 0 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT >= 200.0 && HT  <= 400.0) {
                  mSigPred.at(0)+=weight;
                }
                else if(numbtags == 0 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(1)+=weight;
                }
                else if(numbtags == 0 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(2)+=weight;
                }
                else if(numbtags == 0 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(3)+=weight;
                }
                else if(numbtags == 0 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(4)+=weight;
                }
                else if(numbtags == 0 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(5)+=weight;
                }
                else if(numbtags == 0 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(6)+=weight;
                }
                else if(numbtags == 0 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(7)+=weight;
                }
                //Btag 1 signal regions
                else if(numbtags == 1 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT >= 200.0 && HT  <= 400.0) {
                  mSigPred.at(8)+=weight;
                }
                else if(numbtags == 1 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(9)+=weight;
                }
                else if(numbtags == 1 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(10)+=weight;
                }
                else if(numbtags == 1 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(11)+=weight;
                }
                else if(numbtags == 1 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(12)+=weight;
                }
                else if(numbtags == 1 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(13)+=weight;
                }
                else if(numbtags == 1 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(14)+=weight;
                }
                else if(numbtags == 1 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(15)+=weight;
                }
                //B tag >= 2 signal regions
                else if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT >= 200.0 && HT  <= 400.0) {
                  mSigPred.at(16)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(17)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(18)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(19)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(20)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(21)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(22)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(23)+=weight;
                }

                // 	      //use this snippet to perform studies that optimise for a given bin
                // 	      if(MET > 0.0 && HT > 80.0 && numbtags >= 2) {
                // 		mSigPred.at(0)+=weight;
                // 	      }
                // 	      if(MET > 30.0 && HT > 80.0 && numbtags >= 2) {
                // 		mSigPred.at(1)+=weight;
                // 	      }
                // 	      if(MET > 120.0 && HT > 200.0 && numbtags >= 2 && numjets >= 4) {
                // 		mSigPred.at(2)+=weight;
                // 	      }
                // 	      if(MET > 50.0 && HT > 200.0 && numbtags >= 2 && numjets >= 4) {
                // 		mSigPred.at(3)+=weight;
                // 	      }
                // 	      if(MET > 50.0 && HT > 320.0 && numbtags >= 2 && numjets >= 4) {
                // 		mSigPred.at(4)+=weight;		
                // 	      }
                // 	      if(MET > 120.0 && HT > 320.0 && numbtags >= 2 && numjets >= 4) {
                // 		mSigPred.at(5)+=weight;
                // 	      }
                // 	      if(MET > 50.0 && HT > 200.0 && numbtags >= 3 && numjets >= 3) {
                // 		mSigPred.at(6)+=weight;
                // 	      }
                // 	      if(MET > 0.0 && HT > 320.0 && numbtags >=2 && numjets >= 4) {
                // 		mSigPred.at(7)+=weight;
                // 	      }
              }
              else if(signalRegions==SRbTag2){
                //B tag >= 2 signal regions
                if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT >= 200.0 && HT  <= 400.0) {
                  mSigPred.at(0)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(1)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(2)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET >= 50.0 && MET <= 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(3)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(4)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=2 && numjets <=3  &&
                    HT  > 400.0) {
                  mSigPred.at(5)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  >= 200.0 && HT <=400.0) {
                  mSigPred.at(6)+=weight;
                }
                else if(numbtags >= 2 &&  
                    MET > 120.0 && 
                    numjets >=4  &&
                    HT  > 400.0) {
                  mSigPred.at(7)+=weight;
                }
              }
            }
          }
        }
      }
    }
  }

  return;

}
