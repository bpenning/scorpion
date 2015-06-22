#ifndef D3Reader_h
#define D3Reader_h

/** \class D3Reader
 *
 *  Class simplifying access to ROOT tree branches
 *
 *  $Date: 2008-06-04 13:57:27 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "TROOT.h"
#include "TNamed.h"
#include "TChain.h"
#include "TFile.h"
#include "Reader.hh"
#include "jad_jet_class.hh"
#include "jad_particle_class.hh"
#include "jad_lepton_class.hh"
#include "jad_photon_class.hh"
#include "jad_object_class.hh"
#include "jad_track_class.hh"

#include <map>

class D3Reader : public Reader
{
public :

  D3Reader(TTree *tree);
  D3Reader(TTree *tree,bool gen);
  D3Reader();
  ~D3Reader();

  void SetTree(TTree *tree) { fChain = tree; }

  Long64_t GetEntries() const { return fChain ? static_cast<Long64_t>(fChain->GetEntries()) : 0; }
  Bool_t ReadEntry(Long64_t entry);

  TClonesArray *UseBranch(const char *branchName);
  std::vector<jjet> GetJet() const;
  std::vector<jparticle> GetGenParticle() const;
  std::vector<jtrack> GetIsoChargedTrack() const;
  std::vector<jjet> GetTauJet() const;
  std::vector<jjet> GetMet() const;
  std::vector<jlepton> GetElec() const;
  std::vector<jlepton> GetMuon() const;
  std::vector<jphoton> GetPhoton() const;
  std::vector<double> GetScalarHT() const;
  double GetWeight() const;

private:

  Bool_t Notify();

  TTree *fChain; //! pointer to the analyzed TTree or TChain
  Int_t fCurrentTree; //! current Tree number in a TChain

  typedef std::map<TString, std::pair<TBranch*, TClonesArray*> > TBranchMap;

  TBranchMap fBranchMap; //!

  const TClonesArray * JET;   
  const TClonesArray * TAUJET;
  const TClonesArray * CTRACK;
  const TClonesArray * PHOTON;
  const TClonesArray * ELEC; 
  const TClonesArray * MUON; 
  const TClonesArray * ETMIS;
  const TClonesArray * SCALARHT;

  const TClonesArray * GENPARTICLE;
  const TClonesArray * GENEVENT;
//  const TClonesArray * GENEVENT;
};

#endif // D3Reader_h
