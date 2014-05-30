#include "jad_jet_class.hh"

jjet::jjet(){}
jjet::~jjet(){}
jjet::jjet(double Px, double Py, double Pz, double E, bool isbtag, bool istautag) : jobject(Px,Py,Pz,E,mjet) {

mIsBtagged=isbtag;
mIsTautagged=istautag;

}

bool jjet::Btag() const {
    return mIsBtagged;
}
bool jjet::TauTag() const {
    return mIsTautagged;
}

void jjet::setZeroMass(){
   SetXYZM(this->Px(),this->Py(),this->Pz(),0); 
}
