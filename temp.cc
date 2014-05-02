
std::vector<jlepton> SS::getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam) {

  std::vector<jlepton> leptons;

  TIter itElec((TCollection*)ELEC);
  TRootElectron *elec;
  itElec.Reset();

  while( elec = ((TRootElectron*) itElec.Next()) ) {
    //if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
    //
    //double reliso = (elec->SumEt + elec->SumPt)/elec->PT;
    //std::cout << "reliso: " << reliso << std::endl;
    //reliso > riso
    if(elec->PT < pte || !elec->IsolFlag || fabs(elec->Eta) > etae || (fabs(elec->Eta) > 1.4 && fabs(elec->Eta) < 1.6)) continue;
    bool poscharge = true;
    if(elec->Charge < 0) { poscharge = false; }
    jlepton lepton(elec->Px, elec->Py, elec->Pz, elec->E, true, poscharge);
    leptons.push_back(lepton);
  }
  
  TIter itMuon((TCollection*)MUON);
  TRootMuon *muon;
  itMuon.Reset();
  
  while( muon = ((TRootMuon*) itMuon.Next()) ) {
    //double reliso = (muon->SumEt + muon->SumPt)/muon->PT;
    //if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
    //reliso > riso
    if(muon->PT<ptm || !muon->IsolFlag || fabs(muon->Eta) > etam) continue;
    bool poscharge = true;
    if(muon->Charge < 0) { 
      poscharge = false; 
    }
    jlepton lepton(muon->Px, muon->Py, muon->Pz, muon->E, false, poscharge);
    leptons.push_back(lepton);
  }

  std::sort(leptons.begin(), leptons.end(), std::greater<jlepton>()); //operators defined in the jlepton class
 
  return leptons;
}
