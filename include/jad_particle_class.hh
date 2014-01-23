#ifndef __JPARTICLECLASS__
#define __JPARTICLECLASS__

#include "TLorentzVector.h"

class jparticle : public TLorentzVector {
  
public:
  jparticle();
  ~jparticle();
  jparticle(double Px, double Py, double Pz, double E, int PID, int Status, double Charge);

  bool operator<(const jparticle & rhs) const;
  bool operator>(const jparticle & rhs) const;
  
  double Charge() const;
  int PID() const;
  int Status() const;

private:

  double mCharge;
  int mPID;
  int mStatus;

};

#endif
