#ifndef __DOJFUNCTIONS__
#define __DOJFUNCTIONS__

#include <iostream>
#include <vector>
#include "alphat_functions.hh"
#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"


enum jcharge { CHARGEPOSITIVE, CHARGENEGATIVE, CHARGEBOTH};
enum jflavour { FLAVOURELECTRON, FLAVOURMUON, FLAVOURLEPTON};

double getht(const std::vector<jjet> & jets);
double getmht(const std::vector<jjet> & jets);
double getmet(const TClonesArray *ETMISS);
unsigned int getnbtags(const std::vector<jjet> & jets);
std::vector<jjet> getjets(const TClonesArray *JET, const float & pt, const float & eta, const bool & odd=false );
bool checkDRjetlep(const std::vector<jjet> & jets, const std::vector<jlepton> & leptons, const double & DRmin);
std::vector<jjet> getjetsdrsep(const TClonesArray *JET, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim);
std::vector<jlepton> getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam, const jcharge & charge, const jflavour & flavour);
double getAlphaT(const std::vector<jjet> & jets);
//std::vector<jparticle> getparticles(const TClonesArray * GEN, const unsigned int & pidfilter=0);

#endif
