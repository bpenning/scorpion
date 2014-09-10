#ifndef __JPHOTONCLASS__
#define __JPHOTONCLASS__

#include "jad_object_class.hh"

class jphoton : public jobject {
  
public:
  jphoton();
  ~jphoton();
  jphoton(double Px, double Py, double Pz, double E);

};

#endif
