#include "TMath.h"
#include "jad_track_functions.hh"


std::vector<jtrack> trackSkim(const std::vector<jtrack> & track, const float & pte) {
	std::vector<jtrack> tracks;

	for (int ii=0; ii<track.size();ii++)
	{
	    if(track[ii].Pt() < pte || !track[ii].IsolFlag()) continue;
	    tracks.push_back(track[ii]);
	}

	std::sort(tracks.begin(), tracks.end(), std::greater<jtrack>()); //operators defined in the jtrack class
	return tracks;
}

 std::vector<jtrack> trackSkim(const std::vector<jtrack> & track, const double & pte) {

	std::vector<jtrack> tracks;

	for (int ii=0; ii<track.size();ii++)
	{
	    if(track[ii].Pt() < pte || !track[ii].IsolFlag()) continue;
	    tracks.push_back(track[ii]);
	}

	std::sort(tracks.begin(), tracks.end(), std::greater<jtrack>()); //operators defined in the jtrack class
	return tracks;
    }


