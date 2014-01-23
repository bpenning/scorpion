#ifndef __JBLOCKCLASSES__
#define __JBLOCKCLASSES__

#include "TLorentzVector.h"
#include "TObject.h"

const float UNDEFINED=-9999.;

class TRootParticle : public TObject {

public:

  TRootParticle() {};
  float E;  // particle energy in GeV
  float Px; // particle momentum vector (x component) in GeV
  float Py; // particle momentum vector (y component) in GeV
  float Pz; // particle momentum vector (z component) in GeV

  float Eta; // particle pseudorapidity  
  float Phi; // particle azimuthal angle in rad 

  void Set(const TLorentzVector& momentum);
  void Set(const float & px, const float & py, const float & pz, const float & e);
  void SetEtaPhi(const float & eta, const float & phi);
  void SetEtaPhiEET(const float & eta, const float & phi, const float & e, const float & et);
  float PT; // particle transverse momentum in GeV
  ClassDef(TRootParticle, 2)
};
		     
//--------------------------------------------------------------------------

//class TRootGenParticle;

namespace TRootC {
class GenParticle: public TRootParticle {

public:
  GenParticle() {};
  //GenParticle(const TRootGenParticle& p); //JM commented out don't need
  int PID; // particle HEP ID number [RawHepEventParticle::pid()]
  int Status; // particle status [RawHepEventParticle::status()]
  int M1; // particle 1st mother [RawHepEventParticle::mother1() - 1]
  int M2; // particle 2nd mother [RawHepEventParticle::mother2() - 1]
  int D1; // particle 1st daughter [RawHepEventParticle::daughter1() - 1]
  int D2; // particle 2nd daughter [RawHepEventParticle::daughter2() - 1]

  float Charge;

  float T; // particle vertex position (t component) [RawHepEventParticle::t()]
  float X; // particle vertex position (x component) [RawHepEventParticle::x()]
  float Y; // particle vertex position (y component) [RawHepEventParticle::y()]
  float Z; // particle vertex position (z component) [RawHepEventParticle::z()]
  float M;
 

  //static TCompare *fgCompare; //! JM commented out as don't need
  
  ClassDef(GenParticle, 2)
};
}

//---------------------------------------------------------------------------
/*
class TRootGenParticle: public TRootParticle {

public:
  TRootGenParticle() {_initialised=false; M=-9999.; }
  TRootGenParticle(const int pid): PID(pid) {_initialised=false;}
  TRootGenParticle(TRootC::GenParticle* part);

  int PID; // particle HEP ID number [RawHepEventParticle::pid()]
  int Status; // particle status [RawHepEventParticle::status()]
  int M1; // particle 1st mother [RawHepEventParticle::mother1() - 1]
  int M2; // particle 2nd mother [RawHepEventParticle::mother2() - 1]
  int D1; // particle 1st daughter [RawHepEventParticle::daughter1() - 1]
  int D2; // particle 2nd daughter [RawHepEventParticle::daughter2() - 1]

  float T; // particle vertex position (t component) [RawHepEventParticle::t()]
  float X; // particle vertex position (x component) [RawHepEventParticle::x()]
  float Y; // particle vertex position (y component) [RawHepEventParticle::y()]
  float Z; // particle vertex position (z component) [RawHepEventParticle::z()]
  float M;
  void setFractions();
  const float getFem()  {if(!_initialised) setFractions(); return _Fem;}
  const float getFhad() {if(!_initialised) setFractions(); return _Fhad;}

  float EtaCalo; // particle pseudorapidity when entering the calo, 
  float PhiCalo; // particle azimuthal angle in rad when entering the calo
  void SetEtaPhiCalo(const float eta, const float phi) {EtaCalo=eta; PhiCalo=phi;};

  void print();//

  static TCompare *fgCompare; //!

  float Charge;      // electrical charge
 protected:
  float _Fem, _Fhad; // fractions of energy deposit
  bool _initialised;
  ClassDef(TRootGenParticle, 1)
};
*/

//------------------------------------------------------------------------------

class TRootElectron : public TRootParticle {
public:
  TRootElectron():Charge(-999), IsolFlag(false), EtaCalo(UNDEFINED), PhiCalo(UNDEFINED), EHoverEE(UNDEFINED){};

  int Charge; // particle Charge [RawHepEventParticle::pid()]
  bool IsolFlag; // stores the result of the isolation test
  float EtaCalo; // particle pseudorapidity when entering the calo, 
  float PhiCalo; // particle azimuthal angle in rad when entering the calo
  float EHoverEE;
  float EtRatio;
  float SumEt;
  float SumPt;

  void SetEtaPhiCalo(const float & eta, const float & phi);
  ClassDef(TRootElectron, 2)
};

//------------------------------------------------------------------------------

class TRootPhoton : public TRootParticle {
public:
  TRootPhoton() : EHoverEE(UNDEFINED) {};

  float EHoverEE;
  ClassDef(TRootPhoton, 2)
};

//------------------------------------------------------------------------------

class TRootMuon : public TRootParticle {
public:
  TRootMuon():Charge(-999), IsolFlag(false), EtaCalo(UNDEFINED), PhiCalo(UNDEFINED), EHoverEE(UNDEFINED), EtRatio(UNDEFINED) {};
  int Charge; // particle Charge [RawHepEventParticle::pid()]
  bool IsolFlag;
  float EtaCalo; // particle pseudorapidity when entering the calo, 
  float PhiCalo; // particle azimuthal angle in rad when entering the calo
  float EHoverEE; // hadronic energy over electromagnetic energy
  float EtRatio;  // calo Et in NxN-tower grid around the muon over the muon Et
  float SumEt;
  float SumPt;

  void SetEtaPhiCalo(const float & eta, const float & phi);
  ClassDef(TRootMuon, 2)
};

//---------------------------------------------------------------------------

class TRootTauJet : public TRootParticle {
public:
  TRootTauJet() {};
  float Charge; // normally, using the charge of the track ; here using gen-level tau charge
  int NTracks;
  int NCalo;

  float EHoverEE;
  //void Set(const TLorentzVector& momentum);// { return TRootParticle::Set(momentum); }
  ClassDef(TRootTauJet, 2)
};

//---------------------------------------------------------------------------

class TRootJet : public TRootParticle {
public:
  TRootJet() {};

  bool Btag;
  int NTracks;
  int NCalo;
  
  float EHoverEE;
  ClassDef(TRootJet, 2)
};

class TRootETmis : public TObject {
public:
  TRootETmis() {};
  float ET;
  float Phi;
  float Px;
  float Py;
  ClassDef(TRootETmis, 2)
};

#endif
