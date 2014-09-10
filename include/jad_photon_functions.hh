#ifndef __JPHOTFUNC
#define __JPHOTFUNC
#include <utility>
#include <iostream>
#include "jad_photon_class.hh"
#include "jad_jet_class.hh"

std::vector <jphoton> goodphotons(const std::vector <jphoton> & photon_collection,
            const double &ptcut, const double & etacut, const float & etaw1, const float & etaw2);

#endif
