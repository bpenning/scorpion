#ifndef __JLEPFUNC
#define __JLEPFUNC
#include "jad_lepton_class.hh"
#include "jad_jet_class.hh"
#include <utility>

std::vector <jlepton> goodleptons(const std::vector <jlepton> & lepton_collection,
            const double &ptcut, const double & etacut, const float & etaw1, const float & etaw2);
std::vector <jlepton> goodleptonsNoIso(const std::vector <jlepton> & lepton_collection,
            const double &ptcut, const double & etacut, const float & etaw1, const float & etaw2);
std::vector <jlepton> goodleptons(const std::vector <jlepton> & lepton_collection,double ptcut, double etacut);
std::vector <jlepton> badleptons(std::vector <jlepton> & lepton_collection,double ptcut, double etacut);

std::vector<jlepton> leptonSkim(std::vector<jlepton> & elec,std::vector<jlepton> & muon, const float & pte, const float & etae, const float & ptm, const float & etam, const float & etaw1,const float & etaw2);
std::vector<jlepton> leptonSkim(const std::vector<jlepton> & elec,const std::vector<jlepton> & muon, const double & pte, const double & etae, const double & ptm, const double & etam);
std::vector<std::pair<jlepton,jlepton> > get_ossf_pairs(std::vector<jlepton> leptons);
#endif
