#include "jad_particle_class.hh"

jparticle::jparticle() {}
jparticle::~jparticle() {}
jparticle::jparticle(double Px, double Py, double Pz, double E, int PID, int Status, double Charge) : TLorentzVector(Px, Py, Pz, E), mPID(PID), mStatus(Status), mCharge(Charge) {}

double jparticle::Charge() const {
  return mCharge;
}

int jparticle::PID() const {
  return mPID;
}

int jparticle::Status() const {
  return mStatus;
}

bool jparticle::operator<(const jparticle & rhs) const {
  return this->Pt() < rhs.Pt();
}

bool jparticle::operator>(const jparticle & rhs) const {
  return this->Pt() > rhs.Pt();
}
