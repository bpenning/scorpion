#include "jad_jet_class.hh"

jjet::jjet() {}
jjet::~jjet() {}
jjet::jjet(double Px, double Py, double Pz, double E, bool isbtag) : TLorentzVector(Px, Py, Pz, E), mIsBtagged(isbtag) {}

bool jjet::Btag() const {
  return mIsBtagged;
}

bool jjet::operator<(const jjet & rhs) const {
  return this->Pt() < rhs.Pt();
}

bool jjet::operator>(const jjet & rhs) const {
  return this->Pt() > rhs.Pt();
}
