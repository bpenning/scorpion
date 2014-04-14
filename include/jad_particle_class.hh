#ifndef __JPARTICLECLASS__
#define __JPARTICLECLASS__

#include "jad_object_class.hh"

class jparticle : public jobject {
  
public:
  jparticle();
  ~jparticle();
  jparticle(double Px, double Py, double Pz, double E, int PID, int Status, double Charge);

  double Charge() const;
  int PID() const;
  int Status() const;

private:

  double mCharge;
  int mPID;
  int mStatus;

};

#endif
