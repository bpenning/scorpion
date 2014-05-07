#ifndef __JJETFUNC__
#define __JJETFUNC__

#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
std::vector <jjet> goodjetsSkim(const std::vector <jjet> & jet_collection,double ptcut, double etacut);
std::vector <jjet> badjetsSkim(const std::vector <jjet> & jet_collection,double ptcut, double etacut);

std::vector<jjet> goodjetsSkimDRcut(const std::vector<jjet> & jet_collection, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim);

bool deltaphiptjet(const std::vector<jjet> & jets, const double & px, const double & py, const double & deltaphi, const int & jetcuts, double & dphiret, double & phiret);
unsigned int getnbtags(const std::vector<jjet> & jets);
double getht(const std::vector<jjet> & jets);
double getmet(const std::vector<jjet> &etmvec);
#endif
