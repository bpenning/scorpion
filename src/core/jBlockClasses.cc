#include "jBlockClasses.hh"

void TRootParticle::Set(const TLorentzVector& momentum) {
  
  E   = momentum.E();
  Px  = momentum.Px();
  Py  = momentum.Py();
  Pz  = momentum.Pz();
  PT  = momentum.Pt();
  Eta = momentum.Eta();
  Phi = momentum.Phi();
  return;

}

void TRootParticle::Set(const float & px, const float & py, const float & pz, const float & e) {
  
  //fill TLorentzVector for convenience functions Pt() Eta() Phi()
  TLorentzVector toFill;
  toFill.SetPxPyPzE(px,py,pz,e);
  E   = e;
  Px  = px;
  Py  = py;
  Pz  = pz;
  PT  = toFill.Pt();
  Eta = toFill.Eta();
  Phi = toFill.Phi();
  return;
  
}

void TRootParticle::SetEtaPhi(const float & eta, const float & phi) {
  
  Eta = eta;
  Phi = phi;
  return;
  
}


void TRootParticle::SetEtaPhiEET(const float & eta, const float & phi, const float & e, const float & et) {
  
  Eta = eta; 
  Phi = phi; 
  E = e; 
  PT = et;
  Px = PT*cos(Phi);
  Py = PT*sin(Phi);
  return;
  
}

void TRootElectron::SetEtaPhiCalo(const float & eta, const float & phi) {
  
  EtaCalo = eta; 
  PhiCalo = phi;
  return;
  
}

void TRootMuon::SetEtaPhiCalo(const float & eta, const float & phi) {
  
  EtaCalo = eta; 
  PhiCalo = phi;
  return;

};

/*
namespace TRootC {
  GenParticle::GenParticle(const TRootGenParticle& p) :
    PID(p.PID), Status(p.Status), M1(p.M1), M2(p.M2), D1(p.D1), D2(p.D2), Charge(p.Charge),
    T(p.T), X(p.X), Y(p.Y), Z(p.Z), M(p.M) {}
}

TRootGenParticle::TRootGenParticle(TRootC::GenParticle* part) :
PID(part->PID),Status(part->Status),M1(part->M1),M2(part->M2),D1(part->D1),D2(part->D2),
T(part->T),X(part->X),Y(part->Y),Z(part->Z),M(part->M){
    E=part->E;
    Px=part->Px;
    Py=part->Py;
    Pz=part->Pz;
    Eta=part->Eta;
    Phi=part->Phi;
    PT=part->PT;
    _initialised=false;
}


void TRootGenParticle::setFractions() {
  switch(abs(PID)) {
    default:  _Fem = 0; _Fhad=1; break;
    case(pE):
    case(pGAMMA):
    case(pPI0):
        _Fem = 1; _Fhad=0; break;
    case(pNU1):
    case(pNU2):
    case(pNU3):
    case(pMU):
        _Fem =0; _Fhad=0; break;
    case(pK0S):
    case(pLAMBDA):
	_Fem=0.3; _Fhad=0.7; break;
  }
  _initialised=true;
}


void TRootGenParticle:: print(){
  cout << "pid = " << PID << "\tM=" << M << "\tQ= " << Charge << "\tStatus = " << Status << "\t E= " << E << endl;
//  cout << "T   = " << T << "\t X = " << X << "\t Y = " << Y << "\t Z = " << Z << endl;
}
*/

/*
void TRootTauJet::Set(const TLorentzVector& momentum) {

        E   = momentum.E();
        Px  = momentum.Px();
        Py  = momentum.Py();
        Pz  = momentum.Pz();
        PT  = momentum.Pt();
        Eta = momentum.Eta();
        Phi = momentum.Phi();

}
*/
