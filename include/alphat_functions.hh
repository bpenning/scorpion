#ifndef __ALPHATFUNCTIONS__
#define __ALPHATFUNCTIONS__

#include <iostream>
#include <numeric>

#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "jad_object_collection.hh"
#include "jad_lepton_class.hh"
#include "jad_jet_class.hh"
#include "jad_object_collection.hh"


struct energy_sums{
  std::vector<double> etvec;
  std::vector<double> pxvec;
  std::vector<double> pyvec;
  double total_ht;
  double total_mht;
  unsigned int njets;
  unsigned int nbtags;
  bool pass_quality_cuts;
};

double alphat( const std::vector<double>& et,
	       const std::vector<double>& px,
	       const std::vector<double>& py,
	       std::vector<bool>& pseudo_jet1,
	       bool list);

double biasedDPhi(std::vector<jjet> inJets);
energy_sums make_energy_sums(const std::vector<jjet> & ht275, const std::vector<jjet> & ht325, const std::vector<jjet> & ht375);

energy_sums make_energy_sums(const std::vector<jjet> & ht40, const std::vector<jjet> & ht100);

energy_sums make_energy_sums_20(const std::vector<jjet> & ht200, const std::vector<jjet> & ht275, const std::vector<jjet> & ht325, const std::vector<jjet> & ht375);
#endif
