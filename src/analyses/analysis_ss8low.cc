#include "analysis_ss8low.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
SS8low::SS8low(const std::string & name, 
	   const std::string & experiment,
	   const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
SS8low::SS8low(const std::string & name, 
	   const std::string & experiment, 
	   const unsigned int & numBins,
	   const double & intlumi, 
	   //const std::vector<int> & datayields,
	   const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

SS8low::SS8low(const std::string & name, 
	   const std::string & experiment, 
	   const unsigned int & numBins,
	   const double & intlumi, 
	   const std::vector<double> & bgpred,
	   const std::vector<double> & bgpreduncert,
	   const std::vector<int> & datayields,
	   const std::string & fitmode,
	   const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


SS8low::~SS8low() {}

void SS8low::initHistos() {
  andir->cd();
  ptleadinglep = new TH1D("ptleadinglep",";leading p_{T} [GeV];",100, -5.0, 995.0);
  ptchargeonevschargetwo = new TH2D("ptchargeonevschargetwo",";charge leading lep;charge sub-leading lep",2, -0.5, 1.5, 2, -0.5, 1.5);
  njethist = new TH1D("njethist",";N_{Jets};",20, -0.5, 19.5);
  hthist = new TH1D("hthist",";H_{T} [GeV];",250, -5.0, 2495.0);
  methist = new TH1D("methist",";E_{T}^{miss} [GeV];", 100, -5.0, 995.0);
  leadingjetpt = new TH1D("leadingjetpt",";leading-jet p_{T} [GeV];",100, -5.0, 995.0);
  numbjets = new TH1D("numbjets",";#b-tags;", 10, -0.5, 9.5);
  sfinvmass = new TH1D("sfinvmass",";M_{SF} [GeV];", 1000, -0.5, 999.5);
  drjetcombo = new TH1D("drjetcombo",";#DeltaR;",500, -0.005, 4.995);
}

void SS8low::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  std::vector<jlepton> elecs = treereader->GetElec();
  std::vector<jlepton> muons = treereader->GetMuon();
  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> leptons = leptonSkim(elecs, muons, 5.0, 2.4, 5.0, 2.4, 1.4442, 1.566);
  std::vector<jjet> goodjets = goodjetsSkim(treereader->GetJet(), 40.0, 2.4);

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
      if(leptons[0].Pt() > 10.0 && leptons[1].Pt() > 10.0) { //check both leading leptons have at least 20 GeV
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

              //B tag 0 signal regions
              if(numbtags == 0 &&  
                  MET >= 50.0 && MET <= 120.0 && 
                  numjets >=2 && numjets <=3  &&
                  HT >= 250.0 && HT  <= 400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT >= 250.0 && HT  <= 400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT >= 250.0 && HT  <= 400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
                  HT  >= 250.0 && HT <=400.0) {
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
          }
        }
      }
    }
  }

  return;

}
