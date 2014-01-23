#ifndef __PAIRINFO__
#define __PAIRINFO__

#include <iostream>
#include <vector>
#include <algorithm>

#include "jad_jet_class.hh"
#include "jad_particle_class.hh"

class pair_info {

public:
  pair_info(int id1, int id2, double dr, int pidid1);
  int ID1();
  int ID2();
  int PIDID1();
  double DR();
  void Print();

private:
  int mId1;
  int mId2;
  double mDr;
  int mPIDID1;
};

bool sortDR(pair_info vec1, pair_info vec2);
std::vector<pair_info> make_pairs(const std::vector<jparticle> & gen_vec, const std::vector<jjet> & reco_vec);
std::vector<int> analyse_pairs(std::vector<pair_info> & pairs, int reco_size, double DRmax);
bool in_array(std::vector<int> array, int value);

#endif
