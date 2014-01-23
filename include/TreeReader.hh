#ifndef __TREEREADER__
#define __TREEREADER__

#include "TROOT.h"
#include "TTree.h"
#include "TString.h"

#include <map>

class TreeReader {
public :

  TreeReader(TTree *tree);
  ~TreeReader();

  Long64_t GetEntries() const { return fChain ? static_cast<Long64_t>(fChain->GetEntries()) : 0; }
  Bool_t ReadEntry(Long64_t entry);

  TClonesArray *UseBranch(const TString &branchName);

  const TClonesArray * Jet() const;
  const TClonesArray * Elec() const;
  const TClonesArray * Muon() const;
  const TClonesArray * ETMis() const;
  const TClonesArray * GenParticles() const;
  const TClonesArray * GenEvent() const;

private:

  void Init();
  void Notify();

  TTree *fChain;  // pointer to the analyzed TTree or TChain
  Int_t fCurrentTree; // current Tree number in a TChain

  Bool_t fIsInitDone;

  typedef std::map<TString, std::pair<TBranch*, TClonesArray*> > TBranchMap;

  TBranchMap fBranchMap;
  
  const TClonesArray * JET;   
  //const TClonesArray * TAUJET;
  //const TClonesArray * PHOTON;
  const TClonesArray * ELEC; 
  const TClonesArray * MUON; 
  const TClonesArray * ETMIS;

  const TClonesArray * GENPARTICLE;
  const TClonesArray * GENEVENT;

};

#endif
