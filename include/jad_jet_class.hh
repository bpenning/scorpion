#ifndef __JJETCLASS__
#define __JJETCLASS__

#include "TLorentzVector.h"

class jjet : public TLorentzVector {
  
public:
  jjet();
  ~jjet();
  jjet(double Px, double Py, double Pz, double E, bool isbtag);

  bool operator<(const jjet & rhs) const;
  bool operator>(const jjet & rhs) const;
  bool Btag() const;

private:

  bool mIsBtagged;

};

#endif
