#ifndef __D2READER__
#define __D2READER__

#include "TROOT.h"
#include "TTree.h"
#include "TString.h"
#include "jad_jet_class.hh"
#include "jad_particle_class.hh"
#include "jad_lepton_class.hh"
#include "jad_track_class.hh"
#include "jad_photon_class.hh"
#include "jad_object_class.hh"
#include "jBlockClasses.hh"
#include "Reader.hh"
#include <map>

class D2Reader : public Reader {

public :

  D2Reader(TTree *tree);
  ~D2Reader();

  Long64_t GetEntries() const;
  Bool_t ReadEntry(Long64_t entry);

  TClonesArray *UseBranch(const TString &branchName);

  std::vector<jjet> GetJet() const;
  std::vector<jjet> GetTauJet() const;
  std::vector<jlepton> GetElec() const;
  std::vector<jlepton> GetMuon() const;
  std::vector<jtrack> GetIsoChargedTrack() const;
  std::vector<jphoton> GetPhoton() const;
  std::vector<jjet> GetMet() const;
  std::vector<jparticle> GetGenParticle() const;
  double GetWeight() const;

private:

  void Init();
  void Notify();

  TTree *fChain;  // pointer to the analyzed TTree or TChain
  Int_t fCurrentTree; // current Tree number in a TChain

  Bool_t fIsInitDone;

  typedef std::map<TString, std::pair<TBranch*, TClonesArray*> > TBranchMap;

  TBranchMap fBranchMap;
  
  const TClonesArray * JET;   
  const TClonesArray * TAUJET;
  const TClonesArray * PHOTON;
  const TClonesArray * CTRACK;
  const TClonesArray * ELEC; 
  const TClonesArray * MUON; 
  const TClonesArray * ETMIS;

  const TClonesArray * GENPARTICLE;
  const TClonesArray * GENEVENT;

};

#endif
