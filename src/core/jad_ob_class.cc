#include "jad_object_class.hh"


jobject::jobject(){}

jobject::~jobject(){}

jobject::jobject(double Px, double Py, double Pz, double E, jobjectType type) : TLorentzVector(Px,Py,Pz,E)
{
    //TLorentzVector(Px,Py,Pz,E);
    jtype=type;
}

jobject::jobjectType jobject::gettype()
{
    return jtype;
}

bool jobject::operator<(const jobject & rhs) const {
    return this->Pt() < rhs.Pt();
}

bool jobject::operator>(const jobject & rhs) const {
    return this->Pt() > rhs.Pt();
}

