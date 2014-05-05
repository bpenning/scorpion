#include "analysis_atlas5.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
ATLAS5::ATLAS5(const std::string & name, 
	       const std::string & experiment,
	       const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
ATLAS5::ATLAS5(const std::string & name, 
	       const std::string & experiment, 
	       const unsigned int & numBins,
	       const double & intlumi, 
	       //const std::vector<int> & datayields,
	       const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

ATLAS5::ATLAS5(const std::string & name, 
	       const std::string & experiment, 
	       const unsigned int & numBins,
	       const double & intlumi, 
	       const std::vector<double> & bgpred,
	       const std::vector<double> & bgpreduncert,
	       const std::vector<int> & datayields,
	       const std::string & fitmode,
	       const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


ATLAS5::~ATLAS5() {}



/*bool ATLAS5::jetdrsep(const TSimpleArray<TRootJet> & jets, const double & DR) {

  bool failDR = false; //variable to tell us if any jet pair separation is <DR
  double deltarjj=100.0; //initialise it to some large number for sanity

  for(unsigned int i=0; i < jets.size(); i++) {
    for(unsigned int j=0; j < jets.size(); j++) {
      if(j>i) { //do this to avoid repeating the same pair

	//Find distance (deltaR) between jets
	deltarjj = TMath::Sqrt( (((jets[i]->Phi)-(jets[j]->Phi))*(((jets[i]->Phi)-(jets[j]->Phi)))) + (((jets[i]->Eta)-(jets[j]->Eta))*(((jets[i]->Eta)-(jets[j]->Eta)))) );
	
	if (deltarjj <= DR) { return false; } //if we find a jet pair with sep<DR, no point carrying on!

      }
    }
  }

  return true; //if we haven't returned false, we must not have any jet pairs with sep<DR!

}
*/



void ATLAS5::initHistos() {
  andir->cd();
  njets = new TH1D("njets", ";N_{jets};Entries",20,-0.5,19.5);
  PHIPT = new TH1D( "Phi", ";#Phi; Entries",140,-7.05, 6.95);
  DELTAPHI = new TH1D( "Delta Phi", ";#Delta #Phi; Entries",140,-7.05, 6.95);
  PHI = new TH1D( "Phijets", ";#Phi; Entries",140,-7.05, 6.95);
  PHI40 = new TH1D( "Phijets40", ";#Phi; Entries",140,-7.05, 6.95);
  hthist = new TH1D("hthist", ";H_{T} [GeV];Entries",300,-5,2995);
  etmisshist = new TH1D("etmisshist", ";E_{T}^{miss} [GeV];Entries",200,-5.,1995.);
  ptdiff = new TH1D("ptdiff", ";P_{T}^{miss} difference [GeV];Entries",200,-.5,199.5);
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T}^{miss} difference [GeV];Entries",200,-.5,199.5);
  etdiff = new TH1D("etdiff", ";E_{T}^{miss} difference [GeV];Entries",200,-.5,199.5);
  ATLASA = new TH1D("meffincA", ";m_{eff} [GeV];Entries",300,-5,2995);
  ATLASAd = new TH1D("meffincAd", ";m_{eff} [GeV];Entries",300,-5.,2995.);
  ATLASB = new TH1D("meffincB", ";m_{eff} [GeV];Entries",300,-5.,2995.);
  ATLASC = new TH1D("meffincC", ";m_{eff} [GeV];Entries",300,-5.,2995.);
  ATLASD = new TH1D("meffincD", ";m_{eff} [GeV];Entries",300,-5.,2995.);
  ATLASE = new TH1D("meffincE", ";m_{eff} [GeV];Entries",300,-5.,2995.);
}

void ATLAS5::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  std::vector<jlepton> eletemp=goodleptons(treereader->Elec(), 20.0, 2.47); //all electrons with pt > 20GeV and eta < 2.47 
  std::vector<jlepton>     muontemp=goodleptons(treereader->Muon(), 10.0, 2.4); //all muons  with pt > 10GeV and eta < 2.47 

  std::vector<jjet>      goodjets=goodjetsDR(treereader->Jet(),20.0,2.8, 0.2, eletemp); //check for jets which we should analyse with pt(>20) and eta cuts and electron proximity veto
  std::vector<jjet>      goodjets40=goodjetsDR(treereader->Jet(),40.0,2.8, 0.2, eletemp); //check for jets to make deltaphi with pt(>40) and eta cuts and electron jet proximity veto
  std::vector<jjet>      cleanjets=goodjetsDR(treereader->Jet(),20.0,2.8, 0.2, eletemp); //Jets used to veto leptons not themselves vetoed by electron jet proximity with pt(>20) and eta cuts
  std::vector<jlepton> ele=goodleptonsDR(treereader->Elec(), 20.0, 2.47, 0.4, cleanjets); //Electron Jets with pt and eta cuts and vetoed by cleanjets proximity 
  std::vector<jlepton>     mu=goodleptonsDR(treereader->Muon(), 10.0, 2.4, 0.4, cleanjets);  //Muon Jets with pt and eta cuts and vetoed by cleanjets proximity 
  std::vector<jjet>    etmis=treereader->ETMis();//Missing transverse energy array 
  
  double total_met = etmis[0].Et();//Get missing transverse energy from array
  double diffe = eletemp.size()-ele.size(); //Number of electron jets vetoed by proximity to cleanjets

  if (ele.size() == 0 && mu.size() == 0 && goodjets.size() > 0) {
    njets->Fill(goodjets.size(), weight);//number of jets surviving so far

    double meff = 0.0; //m effective variable (sum of Etmiss and relevant good jets PT)
    double meffinc = total_met;
    double totalht = 0.0;
    for (int nn =0; nn < goodjets40.size(); nn++) {
      meffinc += goodjets40[nn].Pt();//Effective mass defined as sum of PT of all jets with PT>40
      totalht += goodjets40[nn].Pt();
    }

    //Create Et miss from all calorimeter deposits
    double met_x1 = etmis[0].Px(); ///Make x-comp of calorimeter ET
    double met_y1 = etmis[0].Py();//Make y-comp of calorimeter ET

    etmisshist->Fill(total_met, weight); // Histogram of calorimeter Missing Transverse Energy 
    hthist->Fill(totalht, weight);
    leadingjetpt->Fill(goodjets[0].Pt(), weight);//Leading Jet PT histogram

    if(total_met > 160.0 && goodjets.size() > 0 && goodjets[0].Pt() > 130.0 ) { //Apply Etmiss >130 cut and pt >130GeV for leading jet cut with check that jets survive  
      double dphiret = 10.0;
      double phiret = 10.0;

      //CHANNEL A and A'
      if(goodjets.size() >= 2 && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 2, phiret, dphiret)) {//Apply No of jets cut 

	PHIPT->Fill(phiret, weight);
	PHI->Fill(goodjets[0].Phi(), weight); 

	DELTAPHI->Fill(dphiret, weight);
	if(goodjets[1].Pt() > 60.0) { // Second jet pt > 60 cut
	  meff = total_met + goodjets[0].Pt() + goodjets[1].Pt(); //Make meff variable from PT of leading 2 goodjets
	  if(total_met/meff > 0.3) { //SRA cut on transverse momentum/effective mass
	    ATLASA->Fill(meffinc, weight);
	    if(meffinc > 1400.0) {  
	      mSigPred.at(1)+=weight;  // SRA medium
	    }
	    if(meffinc > 1900.0) {  
	      mSigPred.at(0)+=weight;
	    } //SRA tight
	  }
	  if(total_met/meff > 0.4) { //SRA' cut on transverse momentum/effective mass
	    ATLASAd->Fill(meffinc, weight);
	    if(meffinc > 1200.0) {  
	      mSigPred.at(2)+=weight; //SRA'medium 	      
	    }
	  }
	}
      }

      //SIGNAL REGION B:
      if(goodjets.size() >= 3 && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1].Pt() > 60.0 && goodjets[2].Pt() > 60.0 )  {// second and third jet pt > 60 cut

	  meff = total_met + goodjets[0].Pt() + goodjets[1].Pt() + goodjets[2].Pt(); // Make effective mass variable from PT of leading 3 jets

	  if(total_met/meff > 0.25) {//SRB cut on transverse momentum/effective mass

	    ATLASB->Fill(meffinc, weight);

	    if(meffinc > 1900.0) {  
	      mSigPred.at(3)+=weight; //SRB tight
	    }
	  }
	}
      }

      //SIGNAL REGION C:
      if(goodjets.size() >= 4 && deltaphiptjet(goodjets40, met_x1, met_y1, 0.2, 1000, phiret, dphiret) && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1].Pt() > 60.0 && goodjets[2].Pt() > 60.0 && goodjets[3].Pt() > 60.0) {// second, third and fourth jet pt > 60 cut

	  meff = total_met + goodjets[0].Pt() + goodjets[1].Pt() + goodjets[2].Pt() + goodjets[3].Pt(); //Make effective mass variable from PT of leading 4 jets
	  
	  if(total_met/meff > 0.25) {//SRC cut on transverse momentum/effective mass

	    ATLASC->Fill(meffinc, weight);

	    if(meffinc > 900.0) {
	      mSigPred.at(6)+=weight;//SRC Loose
	      if(meffinc > 1200.0) {
		mSigPred.at(5)+=weight;//SRC Medium
		if(meffinc > 1500.0) {
		  mSigPred.at(4)+=weight;//SRC Tight
		}
	      }
	    }
	  }
	}
      }

      //SIGNAL REGION D:
      if(goodjets.size() >= 5 && deltaphiptjet(goodjets40, met_x1, met_y1, 0.2, 100, phiret, dphiret) && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1].Pt() > 60.0 && goodjets[2].Pt() > 60.0 && goodjets[3].Pt() > 60.0 && goodjets[4].Pt() > 40.0) {// second, third, fourth and fifth jet pt > 60 cut
	  
	  meff = total_met + goodjets[0].Pt() + goodjets[1].Pt() + goodjets[2].Pt() + goodjets[3].Pt() + goodjets[4].Pt();//Make effective mass variable from PT of leading 5 jets

	  if(total_met/meff > 0.2) {//SRD cut on transverse momentum/effective mass

	    ATLASD->Fill(meffinc, weight);

	    if(meffinc > 1500.0) {  
	      mSigPred.at(7)+=weight;//SRD Tight
	    }
	  }
	}	   
      }

      //SIGNAL REGION E:
      if(goodjets.size() >= 6 && deltaphiptjet(goodjets40, met_x1, met_y1, 0.2, 1000, phiret, dphiret) && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1].Pt() > 60.0 && goodjets[2].Pt() > 60.0 && goodjets[3].Pt() > 60.0 && goodjets[4].Pt() > 40.0 && goodjets[5].Pt() > 40.0) {// second, third, fourth, fifth and sixth jet pt > 60 cut

	  meff = total_met + goodjets[0].Pt() + goodjets[1].Pt() + goodjets[2].Pt() + goodjets[3].Pt()+ goodjets[4].Pt()+ goodjets[5].Pt(); //Make effective mass variable from PT of leading 6 jets

	  if(total_met/meff > 0.15) {//SRE cut on transverse momentum/effective mass

	    ATLASE->Fill(meffinc, weight);

	    if(meffinc > 900.0) {
	      mSigPred.at(10)+=weight;//SRE Loose
	      if(meffinc > 1200.0) {
		mSigPred.at(9)+=weight;//SRE Medium
		if(meffinc > 1400.0) {
		  mSigPred.at(8)+=weight;//SRE Tight
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
