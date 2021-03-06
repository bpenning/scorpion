#include "jad_particle_functions.hh"
#include "TMath.h"
std::vector<jparticle> getGenParticles(const std::vector<jparticle> & particle_collection, const unsigned int & PIDfilter) {

  std::vector<jparticle> myparticles;

    for (int ii=0; ii<particle_collection.size();ii++)
{
    bool addparticle = true;
    if(PIDfilter != 0 && abs(particle_collection[ii].PID()) != PIDfilter) {
      addparticle = false;
    }
    
    if(addparticle) {
      if(particle_collection[ii].Status() == 3 && fabs(particle_collection[ii].Eta()) < 2.4 && particle_collection[ii].Pt() > 10.0) {
	    jparticle particle(particle_collection[ii].Px(), particle_collection[ii].Py(), particle_collection[ii].Pz(), particle_collection[ii].E(), particle_collection[ii].PID(), particle_collection[ii].Status(), particle_collection[ii].Charge());
	myparticles.push_back(particle);
      }
    }
  }

  return myparticles;

}

std::vector<jparticle> getGenParticles(const std::vector<jparticle> & 
        particle_collection, const unsigned int & PIDfilter, const int & status, 
        const double & pt, const double & eta) {

  std::vector<jparticle> myparticles;
  for (int ii=0; ii<particle_collection.size();ii++){
    bool addparticle = false;
    if(PIDfilter == 0 || abs(particle_collection[ii].PID()) == PIDfilter) 
      addparticle = true;
  
    if(addparticle) {
      if(particle_collection[ii].Status() == status && 
              fabs(particle_collection[ii].Eta()) < eta && 
              particle_collection[ii].Pt() > pt) {
	    jparticle particle(particle_collection[ii].Px(), 
                particle_collection[ii].Py(), particle_collection[ii].Pz(), 
                particle_collection[ii].E(), particle_collection[ii].PID(), 
                particle_collection[ii].Status(), particle_collection[ii].Charge());
	    myparticles.push_back(particle);
      }
    }
  }
  return myparticles;
}
