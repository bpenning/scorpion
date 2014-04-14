#include "jad_jet_class.hh"

jjet::jjet(){}
jjet::~jjet(){}
jjet::jjet(double Px, double Py, double Pz, double E, bool isbtag) : jobject(Px,Py,Pz,E,mjet) {

mIsBtagged=isbtag;

}

bool jjet::Btag() const {
    return mIsBtagged;
}


