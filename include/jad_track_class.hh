#ifndef __JTRACKCLASS__
#define __JTRACKCLASS__

#include "jad_object_class.hh"

class jtrack : public jobject {
  
public:
  jtrack();
  ~jtrack();
  jtrack(double Px, double Py, double Pz, double E, int charge, bool IsolFlag);

  int Charge() const;
  bool IsolFlag() const; 

private:

  bool mIsolFlag;
  int mCharge;

};

#endif
