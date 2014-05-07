//Implementation of the 0 lepton search SUS13011
//arXiv:1305.2390 v2

#include "analysis_zerolep8.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
ZeroLep8::ZeroLep8(const std::string & name, 
    const std::string & experiment,
    const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

  //constructor for object where you also want to run a limit and so other information
  //is necessary
  ZeroLep8::ZeroLep8(const std::string & name, 
      const std::string & experiment, 
      const unsigned int & numBins,
		       const double & intlumi, 
		       //const std::vector<int> & datayields,
		       const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

ZeroLep8::ZeroLep8(const std::string & name, 
		       const std::string & experiment, 
           const unsigned int & numBins,
           const double & intlumi, 
           const std::vector<double> & bgpred,
           const std::vector<double> & bgpreduncert,
           const std::vector<int> & datayields,
           const std::string & fitmode,
           const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


  ZeroLep8::~ZeroLep8() {}


  void ZeroLep8::initHistos() {
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

void ZeroLep8::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //Get the required variables and fill the histograms

  std::vector<jjet> etmis = treereader->GetMet(); //Missing transverse energy array
  double MET = (etmis.size() == 1) ? etmis[0].E(): -1;
  if(MET < -0.5) std::cout << "MET error in 0 lepton search" << std::endl; 
  calomethist->Fill(MET, weight);


  std::vector<jlepton> electrons = goodleptons(treereader->GetElec(), 10.0, 2.5, 99.,99.);
  std::vector<jlepton> muons = goodleptons(treereader->GetMuon(), 10.0, 2.4, 99.,99.);
  //********************************Careful***********************************
  //The taus are used to emulate the cut on the charged track
  std::vector<jjet> taus = goodjetsSkim(treereader->GetTauJet(), 15.0, 2.4);
  //**************************************************************************

  std::vector<jjet> goodjets = goodjetsSkim(treereader->GetJet(), 50.0, 2.4);
  std::vector<jjet> loosejets = goodjetsSkim(treereader->GetJet(), 30.0, 2.4);
  double HT = getht(goodjets);

  double nBTags = getnbtags(goodjets);

  hthist->Fill(HT, weight);
  btagrate->Fill(nBTags, weight);


  //If the general conditions are passed
  if(goodjets.size() >= 3 && electrons.size() == 0 && muons.size() == 0 &&
      taus.size() == 0 && //This is an emulation of the charged track cut, may not work!
      goodjets[0].Pt() > 70.0 && goodjets[1].Pt() > 70.0 &&
      MET > 125.0 && HT > 400.0 && nBTags >= 1
    ) {

    //Make the deltaPhiHatMin variable, as defined in arXiv:1208.4859

    double deltaPhiHatMin = 999.0;

    //Look at angles for the lead 3 jets of an event
    for(unsigned j1 = 0; j1 < 3; j1++){

      //Loops over Pt > 30 jets, loose jet collection (don't say eta cut, assume 2.4)
      double sigmaT2 = 0.0;
      for(std::vector<jjet>::const_iterator j2 = loosejets.begin(); 
          j2 != loosejets.end(); j2++){

        double sigmaPt = 0.1*j2->Pt(); //Jet energy resolution

        double alpha = j2->DeltaPhi(goodjets[j1]); //Angle between j1 and j2

        double x = (sigmaPt * TMath::Sin(alpha));
        sigmaT2 += x*x; 

      }

      double sigmaDeltaPhi = TMath::ATan(TMath::Sqrt(sigmaT2)/MET);

      double deltaPhi = etmis[0].DeltaPhi(goodjets[j1]); //Azimuthal angle between MET and jet

      double deltaPhiHat = deltaPhi/sigmaDeltaPhi;

      if(deltaPhiHat < deltaPhiHatMin) deltaPhiHatMin=deltaPhiHat;
    }

    //Check the cut on the deltaPhiMin
    if( deltaPhiHatMin > 4.0 ){

      //fill the different signal regions
      if(nBTags == 1) {

        if(HT > 400.0 && HT <= 500.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(0)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(1)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(2)+=weight; 
          else if(MET > 350) mSigPred.at(3)+=weight; 
        }
        else if(HT > 500.0 && HT <= 800.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(4)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(5)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(6)+=weight; 
          else if(MET > 350) mSigPred.at(7)+=weight; 
        }
        else if(HT > 800.0 && HT <= 1000.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(8)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(9)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(10)+=weight; 
          else if(MET > 350) mSigPred.at(11)+=weight; 
        }
        else if(HT > 1000.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(12)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(13)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(14)+=weight; 
          else if(MET > 350) mSigPred.at(15)+=weight; 
        }

      } 
      else if(nBTags == 2) {

        if(HT > 400.0 && HT <= 500.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(16)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(17)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(18)+=weight; 
          else if(MET > 350) mSigPred.at(19)+=weight; 
        }
        else if(HT > 500.0 && HT <= 800.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(20)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(21)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(22)+=weight; 
          else if(MET > 350) mSigPred.at(23)+=weight; 
        }
        else if(HT > 800.0 && HT <= 1000.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(24)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(25)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(26)+=weight; 
          else if(MET > 350) mSigPred.at(27)+=weight; 
        }
        else if(HT > 1000.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(28)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(29)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(30)+=weight; 
          else if(MET > 350) mSigPred.at(31)+=weight; 
        }

      } 
      else if(nBTags >= 3) {

        if(HT > 400.0 && HT <= 500.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(32)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(33)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(34)+=weight; 
          else if(MET > 350) mSigPred.at(35)+=weight; 
        }
        else if(HT > 500.0 && HT <= 800.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(36)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(37)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(38)+=weight; 
          else if(MET > 350) mSigPred.at(39)+=weight; 
        }
        else if(HT > 800.0 && HT <= 1000.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(40)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(41)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(42)+=weight; 
          else if(MET > 350) mSigPred.at(43)+=weight; 
        }
        else if(HT > 1000.0){ 
          if(MET > 125.0 && MET <= 150.0) mSigPred.at(44)+=weight; 
          else if(MET > 150.0 && MET <= 250.0) mSigPred.at(45)+=weight; 
          else if(MET > 250.0 && MET <= 350.0) mSigPred.at(46)+=weight; 
          else if(MET > 350) mSigPred.at(47)+=weight; 
        }

      } 
    }
  }
  return;
}
