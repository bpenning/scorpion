#ifndef __JJETCLASS__
#define __JJETCLASS__

#include "jad_object_class.hh"

class jjet : public jobject {

public:
  jjet();
  jjet(double Px, double Py, double Pz, double E, bool isbtag, bool istautag);
  ~jjet();

  bool Btag() const;
  bool TauTag() const;
  void setZeroMass(); 

private:

  bool mIsBtagged;
  bool mIsTautagged;

};

#endif
