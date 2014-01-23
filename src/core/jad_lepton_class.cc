#include "jad_lepton_class.hh"

jlepton::jlepton() {}
jlepton::~jlepton() {}
jlepton::jlepton(double Px, double Py, double Pz, double E, bool iselectron, bool poscharge) : TLorentzVector(Px, Py, Pz, E), mIsElectron(iselectron), mPosCharge(poscharge) {}

bool jlepton::Charge() const {
  return mPosCharge;
}

std::string jlepton::Flavour() const {

  std::string flavour = "electron";
  if(!mIsElectron) {
    flavour = "muon";
  }
 
  //if(mPosCharge) {
  //  flavour += "+";
  //} else {
  //  flavour += "-";
  //}

  return flavour;

}

bool jlepton::operator<(const jlepton & rhs) const {
  return this->Pt() < rhs.Pt();
}

bool jlepton::operator>(const jlepton & rhs) const {
  return this->Pt() > rhs.Pt();
}
