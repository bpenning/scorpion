#ifndef __READER__
#define __READER__

#include "TROOT.h"
#include "TTree.h"
#include "TString.h"
#include "jad_jet_class.hh"
#include "jad_particle_class.hh"
#include "jad_lepton_class.hh"
#include "jad_photon_class.hh"
#include "jad_object_class.hh"
#include "jBlockClasses.hh"

#include <map>

class Reader {
public :
  virtual Long64_t GetEntries() const =0;
  virtual Bool_t ReadEntry(Long64_t entry)=0;

  virtual std::vector<jjet> GetJet() const = 0;
  virtual std::vector<jjet> GetTauJet() const = 0;
  virtual std::vector<jlepton> GetElec() const = 0;
  virtual std::vector<jlepton> GetMuon() const = 0;
  virtual std::vector<jphoton> GetPhoton() const = 0;
  virtual std::vector<jjet> GetMet() const = 0;
  virtual std::vector<jparticle> GetGenParticle() const = 0;
  virtual double GetWeight() const = 0;

};

#endif
