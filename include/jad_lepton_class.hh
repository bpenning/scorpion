#ifndef __JLEPTONCLASS__
#define __JLEPTONCLASS__

#include "jad_object_class.hh"

class jlepton : public jobject {
  
public:
  jlepton();
  ~jlepton();
  jlepton(double Px, double Py, double Pz, double E, bool iselectron, int charge, bool IsolFlag);

  std::string Flavour() const;
  int Charge() const;
  bool IsolFlag() const; 

private:

  bool mIsElectron;
  bool mIsolFlag;
  int mCharge;

};

#endif
