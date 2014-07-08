include "analysis_monojet8.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
MonoJet8::MonoJet8(const std::string & name, 
		 const std::string & experiment,
		 const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
MonoJet8::MonoJet8(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 //const std::vector<int> & datayields,
		 const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

MonoJet8::MonoJet8(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 const std::vector<double> & bgpred,
		 const std::vector<double> & bgpreduncert,
		 const std::vector<int> & datayields,
		 const std::string & fitmode,
		 const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


MonoJet8::~MonoJet8() {}


bool MonoJet8::checkforsecondjetphi(const std::vector<jjet> & goodjets, const double & dphisep) {

  if(goodjets.size() == 2) {
    
    //check if deltaphi j1, j2 < 2.5
    double phijet1 = TMath::ATan2(goodjets[0].Py(), goodjets[0].Px()); //calculate jet phi
    double phijet2 = TMath::ATan2(goodjets[1].Py(), goodjets[1].Px()); //calculate jet phi
    double dphij12 = fabs(phijet1 - phijet2); //calculate the distance in phi between the jet total momentum
    
    //for two vectors in a 2D plane, they cannot be more than pi away from each other
    if(dphij12 > TMath::Pi()) {
      dphij12 = fabs(dphij12 - (2.0 * TMath::Pi()));
    }
    
    if (dphij12 < dphisep) {
      return true;
    } else {
      return false;
    }
  }
  
  return false;
   
}

void MonoJet8::initHistos() {
  andir->cd();
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  calomet = new TH1D("calomet",";E_{T}^{miss} [GeV];Entries",200,-5.,1995.);
	histo = new TH1D("histo",";E_{Tamiss} [GeV];Entries",4,-0.5,3.5);
}

void MonoJet8::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.size() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> ele=goodleptons(treereader->GetElec(), 10.0, 2.4,1.44,1.56); //the central isolated electrons, pt > PT_ELEC GeV
  std::vector<jlepton> mu=goodleptons(treereader->GetMuon(), 10.0, 2.1,2.5,2.6); //the central isolated muons, pt > PT_MUON GeV
  //TSimpleArray<TRootTau> mu=SubArrayMu(treereader.Tau(), 20.0, 2.3); //the central isolated taus, pt > PT_MUON GeV
  std::vector<jjet> goodjets=goodjetsSkim(treereader->GetJet(), 30.0, 4.5); //check for jets which we should analyse
  std::vector<jjet> etmis = treereader->GetMet(); //Missing transverse energy array
  double calo_met = (etmis.size() == 1) ? etmis[0].E(): -1;

  if(goodjets.size() <= 2 && goodjets.size() > 0 && mu.size() == 0 && ele.size() == 0) {

      if(goodjets[0].Pt() > 110.0 && fabs(goodjets[0].Eta()) < 2.4) {a
	histo->Fill(0);

	  if(checkforsecondjetphi(goodjets, 2.5) || goodjets.size()==1) { 
	histo->Fill(2);

	      leadingjetpt->Fill(goodjets[0].Pt(), weight);
	      calomet->Fill(calo_met, weight);
	      if(calo_met > 250.0) { mSigPred.at(0)+=weight; }
	      if(calo_met > 300.0) { mSigPred.at(1)+=weight; }
	      if(calo_met > 350.0) { mSigPred.at(2)+=weight; }
	      if(calo_met > 400.0) { mSigPred.at(3)+=weight; }
	      if(calo_met > 450.0) { mSigPred.at(4)+=weight; }
	      if(calo_met > 500.0) { mSigPred.at(5)+=weight; }
	      if(calo_met > 550.0) { mSigPred.at(6)+=weight; }

	  }

      }

  }


  return;
}
