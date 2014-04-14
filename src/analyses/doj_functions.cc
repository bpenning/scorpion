#include "doj_functions.hh"

double getht(const std::vector<jjet> & jets) {
  double HT=0.0;
  for(unsigned int i=0;i<jets.size();i++) {
    HT += jets[i].Pt();
  }
  return HT;
}

double getmht(const std::vector<jjet> & jets) {
  double MHTx=0.0;
  double MHTy=0.0;
  for(unsigned int i=0;i<jets.size();i++) {
    MHTx += jets[i].Px();
    MHTy += jets[i].Py();
  }
  return TMath::Sqrt((MHTx * MHTx) + (MHTy * MHTy));
}


double getmet(const TClonesArray *ETMISS) {
  
  std::vector<double> etmvec;

  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();

  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
    etmvec.push_back(etm->ET);
  }

  if(etmvec.size() == 1) {
    return etmvec[0];
  } else {
    std::cout << "more than one MET vector defined in this event!" << std::endl;

  }
  return -1.0;

}

unsigned int getnbtags(const std::vector<jjet> & jets) {

  unsigned int numbtags = 0;

  for(unsigned int i=0; i<jets.size(); i++) {
    if(jets[i].Btag()) {
      numbtags++;
    }
  }
  
  return numbtags;

}

std::vector<jjet> getjets(const TClonesArray *JET, const float & pt, const float & eta, const bool & odd) {

  std::vector<jjet> jets;

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    if(!odd) {
      //check if any jet has Pt>50 and |eta|<3.0
      if(jet->PT > pt && fabs(jet->Eta) < eta) {
	jjet myjet(jet->Px, jet->Py, jet->Pz, jet->E, jet->Btag);
	jets.push_back(myjet);
      }
    } else {
      //check if any jet has Pt>50 and |eta|>3.0
      if(jet->PT > pt && fabs(jet->Eta) > eta) {
	jjet myjet(jet->Px, jet->Py, jet->Pz, jet->E, jet->Btag);
	jets.push_back(myjet);
      }
    }
  }

  std::sort(jets.begin(), jets.end(), std::greater<jjet>()); //operators defined in the jjets class

  return jets;
}

bool checkDRjetlep(const std::vector<jjet> & jets, const std::vector<jlepton> & leptons, const double & DRmin) {

  for(unsigned int i=0; i<jets.size(); i++) {
    for(unsigned int j=0; j<leptons.size(); j++) {
      double deltar = jets[i].DeltaR(leptons[j]);
      //deltarjetlep->Fill(deltar);
      if(deltar < DRmin) {
	return false;
      }

    }
  }

  return true; //i.e. no lepton found within DRmin of a jet

}

std::vector<jjet> getjetsdrsep(const TClonesArray *JET, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim) {

  //this method drops jets with DR<0.4 with any leptons that pass the selection

  std::vector<jjet> jets;

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    //check if any jet has Pt>50 and |eta|<3.0
    if(jet->PT > pt && fabs(jet->Eta) < eta) {
      jjet myjet(jet->Px, jet->Py, jet->Pz, jet->E, jet->Btag);
      bool veto = false;
      for(unsigned int j=0; j<lep.size(); j++) {
	if(myjet.DeltaR(lep[j]) < drlim) {
	  veto = true;
	}
      }
      if(!veto) {
	jets.push_back(myjet);
      }
    }
  }

  std::sort(jets.begin(), jets.end(), std::greater<jjet>()); //operators defined in the jjets class

  return jets;

}

std::vector<jlepton> getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam, const jcharge & charge, const jflavour & flavour) {

  std::vector<jlepton> leptons;

  if(flavour == FLAVOURELECTRON || flavour == FLAVOURLEPTON) {
    TIter itElec((TCollection*)ELEC);
    TRootElectron *elec;
    itElec.Reset();
    
    while( elec = ((TRootElectron*) itElec.Next()) ) {
      //if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
      //
      //double reliso = (elec->SumEt + elec->SumPt)/elec->PT;
      //std::cout << "reliso: " << reliso << std::endl;
      //reliso > riso
      if((elec->Charge > 0 && charge == CHARGENEGATIVE) || (elec->Charge < 0 && charge == CHARGEPOSITIVE)) continue;
      if(elec->PT < pte || !elec->IsolFlag || fabs(elec->Eta) > etae || (fabs(elec->Eta) > 1.4 && fabs(elec->Eta) < 1.6)) continue;
      bool poschargeele = true;
      if(elec->Charge < 0) {
	poschargeele = false;
      }
      jlepton lepton(elec->Px, elec->Py, elec->Pz, elec->E, true, elec->Charge,elec->IsolFlag);
      leptons.push_back(lepton);
    }
  }
  
  if(flavour == FLAVOURMUON || flavour == FLAVOURLEPTON) {
    TIter itMuon((TCollection*)MUON);
    TRootMuon *muon;
    itMuon.Reset();
  
    while( muon = ((TRootMuon*) itMuon.Next()) ) {
      //double reliso = (muon->SumEt + muon->SumPt)/muon->PT;
      //if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
      //reliso > riso
      if((muon->Charge > 0 && charge == CHARGENEGATIVE) || (muon->Charge < 0 && charge == CHARGEPOSITIVE)) continue;
      if(muon->PT<ptm || !muon->IsolFlag || fabs(muon->Eta) > etam) continue;
      bool poschargemu = true;
      if(muon->Charge < 0) {
	poschargemu = false;
      }
      jlepton lepton(muon->Px, muon->Py, muon->Pz, muon->E, false, muon->Charge,muon->IsolFlag);
      leptons.push_back(lepton);
    }
  }
  std::sort(leptons.begin(), leptons.end(), std::greater<jlepton>()); //operators defined in the jlepton class
 
  return leptons;
}



double getAlphaT(const std::vector<jjet> & jets) {
 
  std::vector<double> etvec;
  std::vector<double> pxvec;
  std::vector<double> pyvec;
  std::vector<bool> pseudo;

  for(unsigned int i=0;i<jets.size();i++) {
    etvec.push_back(jets[i].Pt());
    pxvec.push_back(jets[i].Px());
    pyvec.push_back(jets[i].Py());
  }

  return alphat(etvec, pxvec, pyvec, pseudo, true);
}
