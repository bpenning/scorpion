#include "pair_info.hh"

pair_info::pair_info(int id1, int id2, double dr, int pidid1) {
  mId1 = id1;
  mId2 = id2;
  mDr = dr;
  mPIDID1 = pidid1;
}

double pair_info::DR() {
  return mDr;
}

int pair_info::ID1() {
  return mId1;
}

int pair_info::ID2() {
  return mId2;
}

int pair_info::PIDID1() {
  return mPIDID1;
}

void pair_info::Print() {

  std::cout << "ID1 = " << this->ID1() << ", ID2 = " << this->ID2() << ", DeltaR = " << this->DR() << ", PIDID1 = " << this->PIDID1() << std::endl;
  return;

}

bool sortDR(pair_info vec1, pair_info vec2) {

  return (vec1.DR() < vec2.DR());

}

std::vector<pair_info> make_pairs(const std::vector<jparticle> & gen_vec, const std::vector<jjet> & reco_vec) {

  std::vector<pair_info> pinfo;

  for(unsigned int i=0; i<gen_vec.size(); i++) {
    for(unsigned int j=0; j<reco_vec.size(); j++) {
      pinfo.push_back(pair_info(i,j,gen_vec[i].DeltaR(reco_vec[j]),gen_vec[i].PID()));
    }
  }

  std::sort(pinfo.begin(), pinfo.end(), sortDR);

  return pinfo;

}

bool in_array(std::vector<int> array, int value) {

  for(unsigned int i=0; i<array.size(); i++) {
    if(array[i] == value) {
      return true;
    }
  }

  return false;

}

std::vector<int> analyse_pairs(std::vector<pair_info> & pairs, int reco_size, double DRmax) {

  std::vector<int> reco_matched_index(reco_size, 0);

  std::vector<int> gen_matched;
  std::vector<int> reco_matched;

  for(unsigned int i=0; i < pairs.size(); i++) {
    if(!in_array(gen_matched, pairs[i].ID1()) && !in_array(reco_matched, pairs[i].ID2()) && pairs[i].DR() < DRmax) {
      gen_matched.push_back(pairs[i].ID1());
      reco_matched.push_back(pairs[i].ID2());
      reco_matched_index[pairs[i].ID2()] = pairs[i].PIDID1(); //true;//get_flavour(pairs[i].ID1());
    }
  }

  return reco_matched_index;

}
