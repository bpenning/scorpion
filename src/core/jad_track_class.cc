#include "jad_track_class.hh"

jtrack::jtrack(){}
jtrack::~jtrack(){}
jtrack::jtrack(double Px, double Py, double Pz, double E, int charge, bool isolFlag) : jobject(Px,Py,Pz,E,mtrack) {
    mCharge = charge;
    mIsolFlag = isolFlag;
}

int jtrack::Charge() const {
    return mCharge;
}
bool jtrack::IsolFlag() const {
    return mIsolFlag;
}
