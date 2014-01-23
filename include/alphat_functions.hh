#ifndef __ALPHATFUNCTIONS__
#define __ALPHATFUNCTIONS__

#include <iostream>
#include <numeric>

#include "jBlockClasses.hh"
#include "TSimpleArray.hh"

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

energy_sums make_energy_sums(const TSimpleArray<TRootJet> & ht275, const TSimpleArray<TRootJet> & ht325, const TSimpleArray<TRootJet> & ht375);

#endif
