#ifndef __JTRACKFUNC
#define __JTRACKFUNC
#include <utility>
#include <iostream>
#include "jad_track_class.hh"

std::vector<jtrack> trackSkim(const std::vector<jtrack> & track, const float & pte);
std::vector<jtrack> trackSkim(const std::vector<jtrack> & track, const double & pte);
#endif
