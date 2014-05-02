#ifndef __JJETFUNC__
#define __JJETFUNC__

#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
std::vector <jjet> goodjetsSkim(std::vector <jjet> jet_collection,double ptcut, double etacut);
std::vector <jjet> badjetsSkim(std::vector <jjet> jet_collection,double ptcut, double etacut);

std::vector<jjet> goodjetsSkimDRcut(std::vector<jjet> jet_collection, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim);

bool deltaphiptjet(std::vector<jjet> & jets, const double & px, const double & py, const double & deltaphi, const int & jetcuts, double & dphiret, double & phiret);
unsigned int getnbtags(std::vector<jjet> & jets);
double getht( std::vector<jjet> & jets);
double getmet(std::vector<jjet> etmvec);
#endif
