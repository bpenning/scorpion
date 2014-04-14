#include "jad_lepton_class.hh"

jlepton::jlepton() {}
jlepton::~jlepton() {}
jlepton::jlepton(double Px, double Py, double Pz, double E, bool iselectron, int charge, bool IsolFlag): jobject(Px,Py,Pz,E,mlepton)
{
    mIsElectron=iselectron;
    mCharge=charge;
    mIsolFlag=IsolFlag;
}

int jlepton::Charge() const {
    return mCharge;
}
bool jlepton::IsolFlag() const {
    return mIsolFlag;
}

std::string jlepton::Flavour() const {

    std::string flavour = "electron";
    if(!mIsElectron) {
	flavour = "muon";
    }

    //if(mPosCharge) {
    //  flavour += "+";
    //} else {
    //  flavour += "-";
    //}

    return flavour;

}

