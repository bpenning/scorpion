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


bool ATLAS5::jetdrsep(const TSimpleArray<TRootJet> & jets, const double & DR) {

  bool failDR = false; //variable to tell us if any jet pair separation is <DR
  double deltarjj=100.0; //initialise it to some large number for sanity

  for(unsigned int i=0; i < jets.GetEntries(); i++) {
    for(unsigned int j=0; j < jets.GetEntries(); j++) {
      if(j>i) { //do this to avoid repeating the same pair

	//Find distance (deltaR) between jets
	deltarjj = TMath::Sqrt( (((jets[i]->Phi)-(jets[j]->Phi))*(((jets[i]->Phi)-(jets[j]->Phi)))) + (((jets[i]->Eta)-(jets[j]->Eta))*(((jets[i]->Eta)-(jets[j]->Eta)))) );
	
	if (deltarjj <= DR) { return false; } //if we find a jet pair with sep<DR, no point carrying on!

      }
    }
  }

  return true; //if we haven't returned false, we must not have any jet pairs with sep<DR!

}

bool ATLAS5::deltaphiptjet(const TSimpleArray<TRootJet> & jets, const double & px, const double & py, const double & deltaphi, const int & jetcuts, double & dphiret, double & phiret) {

  double phietmiss = TMath::ATan2(py,px); //calculate etmiss phi - do it here as it doesn't change per jet!

  double phijet = 100.0; //initialise to large number for sanity
  double dphitemp = 100.0; //initialise to large number for sanity
  double dphi = 100.0; //initialise to large number for sanity

  for(unsigned int ii=0; ii < jets.GetEntries(); ii++) { //loop over all jet entries

    phijet = TMath::ATan2(jets[ii]->Py, jets[ii]->Px); //calculate jet phi
    dphitemp = fabs(phijet - phietmiss); //calculate the distance in phi between the jet total momentum

    //for two vectors in a 2D plane, they cannot be more than pi away from each other
    if(dphitemp > TMath::Pi()) {
      dphitemp = fabs(dphitemp - (2.0 * TMath::Pi()));
    }
    
    //check for the smallest dphi between jet and etmiss
    if (dphitemp < dphi) {   
      dphi = dphitemp;
    }

    //we only consider the first jetcuts jet - but only if they exist so keep the for loop over the GetEntries
    if (ii == (jetcuts-1)) {
      break;
    }
  }

  //populate some variables for plotting
  phiret = phietmiss;
  dphiret = dphi;
  
  if (dphi <= deltaphi) { return false; } //if we find a jet with sep<DR relative to etmiss, return false
  return true; //otherwise, we haven't found any, so return true

}

TSimpleArray<TRootElectron> ATLAS5::SubArrayEl(const TClonesArray *ELEC, float pt, float eta) {
    TIter itElec((TCollection*)ELEC);
    TRootElectron *elec;
    itElec.Reset();
    TSimpleArray<TRootElectron> array;
      while( (elec = (TRootElectron*) itElec.Next()) )
        {
	  if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
	  array.Add(elec);
        }
    return array;
}

TSimpleArray<TRootElectron> ATLAS5::SubArrayEl2(const TClonesArray *ELEC, float pt, float eta, double DRLmin, const TSimpleArray<TRootJet> &jet) {

  TIter itElec((TCollection*)ELEC);
  TRootElectron *elec;
  itElec.Reset();
  double DRele, DRelemin;
  TSimpleArray<TRootElectron> array;
  while( (elec = (TRootElectron*) itElec.Next()) ) {
    DRelemin = 10.0;
    for(unsigned int i=0; i < jet.GetEntries(); i++) {

      DRele = TMath::Sqrt( (((elec->Phi)-(jet[i]->Phi))*((elec->Phi)-(jet[i]->Phi))) + (((elec->Eta)-(jet[i]->Eta))*((elec->Eta)-(jet[i]->Eta))) );
      if (DRelemin > DRele) {
	DRelemin = DRele;
      }
    }
	 		
    if (DRelemin > DRLmin) { //if the electron is away from all the jets
      //if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
      if(elec->PT>pt && elec->IsolFlag && fabs(elec->Eta) < eta) {
	array.Add(elec);
      }
    }
  }
  
  return array;
}


TSimpleArray<TRootMuon> ATLAS5::SubArrayMu(const TClonesArray *MUON, float pt, float eta) {
   TIter itMuon((TCollection*)MUON);
   TRootMuon *muon;
   itMuon.Reset();
   TSimpleArray<TRootMuon> array;
   while( (muon = (TRootMuon*) itMuon.Next()) )
     {
       if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
       array.Add(muon);
     }
   return array;
}

TSimpleArray<TRootMuon> ATLAS5::SubArrayMu2(const TClonesArray *MUON, float pt, float eta, double DRLmin, const TSimpleArray<TRootJet> &jet) {

  TIter itMuon((TCollection*)MUON);
  TRootMuon *muon;
  itMuon.Reset();
  TSimpleArray<TRootMuon> array;
  double DRelemin, DRele;
  while( (muon = (TRootMuon*) itMuon.Next()) ) {
    DRelemin = 10.0;
    for(unsigned int i=0; i <jet.GetEntries();i++) {

      DRele = TMath::Sqrt( (((muon->Phi)-(jet[i]->Phi))*((muon->Phi)-(jet[i]->Phi))) + (((muon->Eta)-(jet[i]->Eta))*((muon->Eta)-(jet[i]->Eta))) );
      if(DRelemin > DRele) {
	DRelemin = DRele;
      }
      
    }
			
    if (DRelemin > DRLmin) { //if the muon is away from all the jets
      //if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
      if(muon->PT > pt && muon->IsolFlag && fabs(muon->Eta) < eta) {
	array.Add(muon);
      }
    }
  }

  return array;
}

TSimpleArray<TRootJet> ATLAS5::SubArrayGoodJets2(const TClonesArray *JET, float pt, float eta, double DRlim, const TSimpleArray<TRootElectron> &ele ) { 
  
  TIter itJet((TCollection*)JET); //define JET iterator
  TRootJet *jet; //pointer to TRootJet object
  itJet.Reset(); //set it back to the start
  TSimpleArray<TRootJet> array;
  double DRele, DRelemin;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    
    DRele = 0.0;
    DRelemin = 10.0;
    
    for (unsigned int ee=0; ee < ele.GetEntries(); ee++) {
      DRele = TMath::Sqrt( (((ele[ee]->Phi)-(jet->Phi))*((ele[ee]->Phi))-(jet->Phi)) + (((ele[ee]->Eta)-(jet->Eta))*((ele[ee]->Eta)-(jet->Eta))) ); 
      if(DRelemin > DRele) {
	DRelemin = DRele;
      }
    }

    //check if any jet has Pt>50 and |eta|<3.0 and if the min is above threshold between electron and jet (otherwise it is a jet)
    if(jet->PT > pt && fabs(jet->Eta) < eta && DRelemin > DRlim) {
      array.Add(jet); 
    }
  } 
  return array; 
}

TSimpleArray<TRootETmis> ATLAS5::makeETM(const TClonesArray *ETMISS) {
  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();
  TSimpleArray<TRootETmis> array;
  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
      array.Add(etm);
  }
  return array;
}


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

void ATLAS5::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  TSimpleArray<TRootElectron> eletemp=SubArrayEl(treereader.Elec(), 20.0, 2.47); //all electrons with pt > 20GeV and eta < 2.47 
  TSimpleArray<TRootMuon>     muontemp=SubArrayMu(treereader.Muon(), 10.0, 2.4); //all muons  with pt > 10GeV and eta < 2.47 
  TSimpleArray<TRootJet>      goodjets=SubArrayGoodJets2(treereader.Jet(),20.0,2.8, 0.2, eletemp); //check for jets which we should analyse with pt(>20) and eta cuts and electron proximity veto
  TSimpleArray<TRootJet>      goodjets40=SubArrayGoodJets2(treereader.Jet(),40.0,2.8, 0.2, eletemp); //check for jets to make deltaphi with pt(>40) and eta cuts and electron jet proximity veto
  TSimpleArray<TRootJet>      cleanjets=SubArrayGoodJets2(treereader.Jet(),20.0,2.8, 0.2, eletemp); //Jets used to veto leptons not themselves vetoed by electron jet proximity with pt(>20) and eta cuts
  TSimpleArray<TRootElectron> ele=SubArrayEl2(treereader.Elec(), 20.0, 2.47, 0.4, cleanjets); //Electron Jets with pt and eta cuts and vetoed by cleanjets proximity 
  TSimpleArray<TRootMuon>     mu=SubArrayMu2(treereader.Muon(), 10.0, 2.4, 0.4, cleanjets);  //Muon Jets with pt and eta cuts and vetoed by cleanjets proximity 
  TSimpleArray<TRootETmis>    etmis=makeETM(treereader.ETMis());//Missing transverse energy array 
  
  double total_met = etmis[0]->ET;//Get missing transverse energy from array
  double diffe = eletemp.GetEntries()-ele.GetEntries(); //Number of electron jets vetoed by proximity to cleanjets

  if (ele.GetEntries() == 0 && mu.GetEntries() == 0 && goodjets.GetEntries() > 0) {
    njets->Fill(goodjets.GetEntries(), weight);//number of jets surviving so far

    double meff = 0.0; //m effective variable (sum of Etmiss and relevant good jets PT)
    double meffinc = total_met;
    double totalht = 0.0;
    for (int nn =0; nn < goodjets40.GetEntries(); nn++) {
      meffinc += goodjets40[nn]->PT;//Effective mass defined as sum of PT of all jets with PT>40
      totalht += goodjets40[nn]->PT;
    }

    //Create Et miss from all calorimeter deposits
    double met_x1 = etmis[0]->Px; ///Make x-comp of calorimeter ET
    double met_y1 = etmis[0]->Py;//Make y-comp of calorimeter ET

    etmisshist->Fill(total_met, weight); // Histogram of calorimeter Missing Transverse Energy 
    hthist->Fill(totalht, weight);
    leadingjetpt->Fill(goodjets[0]->PT, weight);//Leading Jet PT histogram

    if(total_met > 160.0 && goodjets.GetEntries() > 0 && goodjets[0]->PT > 130.0 ) { //Apply Etmiss >130 cut and pt >130GeV for leading jet cut with check that jets survive  
      double dphiret = 10.0;
      double phiret = 10.0;

      //CHANNEL A and A'
      if(goodjets.GetEntries() >= 2 && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 2, phiret, dphiret)) {//Apply No of jets cut 

	PHIPT->Fill(phiret, weight);
	PHI->Fill(goodjets[0]->Phi, weight); 

	DELTAPHI->Fill(dphiret, weight);
	if(goodjets[1]->PT > 60.0) { // Second jet pt > 60 cut
	  meff = total_met + goodjets[0]->PT + goodjets[1]->PT; //Make meff variable from PT of leading 2 goodjets
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
      if(goodjets.GetEntries() >= 3 && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1]->PT > 60.0 && goodjets[2]->PT > 60.0 )  {// second and third jet pt > 60 cut

	  meff = total_met + goodjets[0]->PT + goodjets[1]->PT + goodjets[2]->PT; // Make effective mass variable from PT of leading 3 jets

	  if(total_met/meff > 0.25) {//SRB cut on transverse momentum/effective mass

	    ATLASB->Fill(meffinc, weight);

	    if(meffinc > 1900.0) {  
	      mSigPred.at(3)+=weight; //SRB tight
	    }
	  }
	}
      }

      //SIGNAL REGION C:
      if(goodjets.GetEntries() >= 4 && deltaphiptjet(goodjets40, met_x1, met_y1, 0.2, 1000, phiret, dphiret) && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1]->PT > 60.0 && goodjets[2]->PT > 60.0 && goodjets[3]->PT > 60.0) {// second, third and fourth jet pt > 60 cut

	  meff = total_met + goodjets[0]->PT + goodjets[1]->PT + goodjets[2]->PT + goodjets[3]->PT; //Make effective mass variable from PT of leading 4 jets
	  
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
      if(goodjets.GetEntries() >= 5 && deltaphiptjet(goodjets40, met_x1, met_y1, 0.2, 100, phiret, dphiret) && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1]->PT > 60.0 && goodjets[2]->PT > 60.0 && goodjets[3]->PT > 60.0 && goodjets[4]->PT > 40.0) {// second, third, fourth and fifth jet pt > 60 cut
	  
	  meff = total_met + goodjets[0]->PT + goodjets[1]->PT + goodjets[2]->PT + goodjets[3]->PT + goodjets[4]->PT;//Make effective mass variable from PT of leading 5 jets

	  if(total_met/meff > 0.2) {//SRD cut on transverse momentum/effective mass

	    ATLASD->Fill(meffinc, weight);

	    if(meffinc > 1500.0) {  
	      mSigPred.at(7)+=weight;//SRD Tight
	    }
	  }
	}	   
      }

      //SIGNAL REGION E:
      if(goodjets.GetEntries() >= 6 && deltaphiptjet(goodjets40, met_x1, met_y1, 0.2, 1000, phiret, dphiret) && deltaphiptjet(goodjets, met_x1, met_y1, 0.4, 3, phiret, dphiret)) { //Apply No of jets cut

	if(goodjets[1]->PT > 60.0 && goodjets[2]->PT > 60.0 && goodjets[3]->PT > 60.0 && goodjets[4]->PT > 40.0 && goodjets[5]->PT > 40.0) {// second, third, fourth, fifth and sixth jet pt > 60 cut

	  meff = total_met + goodjets[0]->PT + goodjets[1]->PT + goodjets[2]->PT + goodjets[3]->PT+ goodjets[4]->PT+ goodjets[5]->PT; //Make effective mass variable from PT of leading 6 jets

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
