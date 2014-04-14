#ifndef __JJETCLASS__
#define __JJETCLASS__

#include "jad_object_class.hh"

class jjet : public jobject {

public:
  jjet();
  jjet(double Px, double Py, double Pz, double E, bool isbtag);
  ~jjet();

  bool Btag() const;

private:

  bool mIsBtagged;

};

#endif
