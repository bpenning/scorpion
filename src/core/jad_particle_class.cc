#include "jad_particle_class.hh"

jparticle::jparticle() {}
jparticle::~jparticle() {}
jparticle::jparticle(double Px, double Py, double Pz, double E, int PID, int Status, double Charge) : jobject(Px,Py,Pz,E,mparticle)
{
    mPID=PID;
    mStatus=Status;
    mCharge=Charge; 
}

double jparticle::Charge() const {
    return mCharge;
}

int jparticle::PID() const {
    return mPID;
}

int jparticle::Status() const {
    return mStatus;
}

