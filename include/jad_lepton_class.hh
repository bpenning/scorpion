#ifndef __JLEPTONCLASS__
#define __JLEPTONCLASS__

#include "TLorentzVector.h"

class jlepton : public TLorentzVector {
  
public:
  jlepton();
  ~jlepton();
  jlepton(double Px, double Py, double Pz, double E, bool iselectron, bool poscharge);

  bool operator<(const jlepton & rhs) const;
  bool operator>(const jlepton & rhs) const;
  
  std::string Flavour() const;
  bool Charge() const;

private:

  bool mIsElectron;
  bool mPosCharge;

};

#endif
