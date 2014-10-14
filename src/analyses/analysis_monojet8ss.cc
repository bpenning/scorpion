#include "analysis_monojet8ss.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
MonoJet8ss::MonoJet8ss(const std::string & name, 
    const std::string & experiment,
    const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

  //constructor for object where you also want to run a limit and so other information
  //is necessary
  MonoJet8ss::MonoJet8ss(const std::string & name, 
      const std::string & experiment, 
      const unsigned int & numBins,
      const double & intlumi, 
      //const std::vector<int> & datayields,
      const std::vector<double> & bgpred) :
    AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

    MonoJet8ss::MonoJet8ss(const std::string & name, 
	const std::string & experiment, 
	const unsigned int & numBins,
	const double & intlumi, 
	const std::vector<double> & bgpred,
	const std::vector<double> & bgpreduncert,
	const std::vector<int> & datayields,
	const std::string & fitmode,
	const bool & calculateR) :
      AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


      MonoJet8ss::~MonoJet8ss() {}


      bool MonoJet8ss::checkforsecondjetphi(const std::vector<jjet> & goodjets, const double & dphisep) {

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

void MonoJet8ss::initHistos() {
  andir->cd();
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  calomet = new TH1D("calomet",";E_{T}^{miss} [GeV];Entries",200,-5.,1995.);
  event_weight = new TH1D("event_weight",";weight;entries",1000,1940e9,1960e9);
  cutflow = new TH1D("cutflow",";cut;entries",5,-0.5,4.5);
  event_weight = new TH1D("event_weight",";weight_num;entries",1000,1940e9,1960e9);
  njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
  calomet = new TH1D("calomet",";E_{T}^{miss} [GeV];Entries",30,250,1000); 
  event_weight->SetBit(TH1::kCanRebin);
}

void MonoJet8ss::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.size() << std::endl;
  andir->cd();
  //std::cout << weight << std::endl;
  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> ele=goodleptons(treereader->GetElec(), 10.0, 2.5,1.44,1.56); //the central isolated electrons, pt > PT_ELEC GeV
  std::vector<jlepton> mu=goodleptons(treereader->GetMuon(), 10.0, 2.4,99,99); //the central isolated muons, pt > PT_MUON GeV
  //TSimpleArray<TRootTau> mu=SubArrayMu(treereader.Tau(), 20.0, 2.3); //the central isolated taus, pt > PT_MUON GeV
  std::vector<jjet> goodjets=goodjetsSkim(treereader->GetJet(), 60.0, 4.5); //check for jets which we should analyse
  std::vector<jjet> taujets=goodjetsSkim(treereader->GetTauJet(), 20.0, 2.3); //check for jets which we should analyse
  std::vector<jjet> etmis = treereader->GetMet(); //Missing transverse energy array
  double calo_met = (etmis.size() == 1) ? etmis[0].E(): -1;


     
 
  njets->Fill(goodjets.size());


	cutflow->Fill(0.,weight);
        event_weight->Fill(weight); 
  if (calo_met > 250) {
   if(goodjets.size() <= 2 && goodjets.size() > 0 && mu.size() == 0 && ele.size() == 0 && taujets.size() == 0) {
		cutflow->Fill(1.,weight);
      if(goodjets[0].Pt() > 110.0 && fabs(goodjets[0].Eta()) < 2.4) {
			cutflow->Fill(2.,weight);
	  if(checkforsecondjetphi(goodjets, 2.5) || goodjets.size()==1) { 
                        cutflow->Fill(3.,weight);

	      calomet->Fill(calo_met, weight);
	      leadingjetpt->Fill(goodjets[0].Pt(), weight);
	      if(goodjets[0].Pt() > 250.0) { mSigPred.at(0)+=weight; }
	      if(goodjets[0].Pt() > 300.0) { mSigPred.at(1)+=weight; }
	      if(goodjets[0].Pt() > 350.0) { mSigPred.at(2)+=weight; }
	      if(goodjets[0].Pt() > 400.0) { mSigPred.at(3)+=weight;

                                     cutflow->Fill(4., weight); }

	      if(goodjets[0].Pt() > 450.0) { mSigPred.at(4)+=weight; }
	      if(goodjets[0].Pt() > 500.0) { mSigPred.at(5)+=weight; }
	      if(goodjets[0].Pt() > 550.0) { mSigPred.at(6)+=weight; }






      }

    }
   }
  }


  return;
}

