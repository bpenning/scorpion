#ifndef __JPARTFUNC__
#define __JPARTFUNC__
#include "jad_particle_class.hh"

std::vector<jparticle> getGenParticles(const std::vector<jparticle> & 
        particle_collection, const unsigned int & PIDfilter);
std::vector<jparticle> getGenParticles(const std::vector<jparticle> & 
        particle_collection, const unsigned int & PIDfilter, const int & status, 
        const double & pt, const double & eta) ;
#endif
