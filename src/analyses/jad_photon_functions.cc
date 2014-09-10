#include "jad_photon_functions.hh"
#include "TMath.h"
    //Normal Jet Collection (eta and pt cuts)
    std::vector <jphoton> goodphotons(const std::vector <jphoton> & photon_collection,
            const double &ptcut, const double & etacut, const float & etaw1, const float & etaw2)
    {
	std::vector<jphoton> mphotonvec;
	for (int ii=0; ii<photon_collection.size();ii++)
	{
	    if(photon_collection[ii].Pt() < ptcut || fabs(photon_collection[ii].Eta()) > etacut || (fabs(photon_collection[ii].Eta()) > etaw1 && fabs(photon_collection[ii].Eta()) < etaw2)) continue;
	    mphotonvec.push_back(photon_collection[ii]);
	}
	return mphotonvec;
    }


