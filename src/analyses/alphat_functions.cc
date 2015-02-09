#include "alphat_functions.hh"

double alphat( const std::vector<double>& et,
               const std::vector<double>& px,
               const std::vector<double>& py,
               std::vector<bool>& pseudo_jet1,
               bool list) {

  // Clear pseudo-jet container
  pseudo_jet1.clear();
  
  // Momentum sums in transverse plane
  const double sum_et = accumulate( et.begin(), et.end(), 0. );
  const double sum_px = accumulate( px.begin(), px.end(), 0. );
  const double sum_py = accumulate( py.begin(), py.end(), 0. );
  
  // Minimum Delta Et for two pseudo-jets
  double min_delta_sum_et = -1.;
  for ( unsigned i=0; i < unsigned(1<<(et.size()-1)); i++ ) { // iterate through different combinations
    double delta_sum_et = 0.;
    std::vector<bool> jet;
    for ( unsigned j=0; j < et.size(); j++ ) { // iterate through jets
      delta_sum_et += et[j] * ( 1 - 2 * (int(i>>j)&1) );
      if ( list ) { jet.push_back( (int(i>>j)&1) == 0 ); }
    }
    if ( ( fabs(delta_sum_et) < min_delta_sum_et || min_delta_sum_et < 0. ) ) {
      min_delta_sum_et = fabs(delta_sum_et);
      if ( list && jet.size() == et.size() ) {
        pseudo_jet1.resize(jet.size());
        std::copy( jet.begin(), jet.end(), pseudo_jet1.begin() );
      }
    }
  }
  if ( min_delta_sum_et < 0. ) { return 0.; }
  
  // Alpha_T
  return ( 0.5 * ( sum_et - min_delta_sum_et ) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) ) );
  
}

energy_sums make_energy_sums(const std::vector<jjet> & ht275, const std::vector<jjet> & ht325, const std::vector<jjet> & ht375) {

  energy_sums esums;
  //initialise the energy sum quantities
  double total_ht = 0.0;
  double mht_x = 0.0;
  double mht_y = 0.0;
  unsigned int njets = 0;
  unsigned int nbtags = 0;
  std::vector<double> etvec;
  std::vector<double> pxvec;
  std::vector<double> pyvec;

  //first, loop over all jets based on threshold of >50 GeV and calculate HT
  for(unsigned int numjets = 0; numjets < ht375.size(); numjets++) {
    total_ht += ht375.at(numjets).Pt();
    mht_x += ht375.at(numjets).Px();
    mht_y += ht375.at(numjets).Py();    
    etvec.push_back(ht375.at(numjets).Pt());
    pxvec.push_back(ht375.at(numjets).Px());
    pyvec.push_back(ht375.at(numjets).Py());
    if(ht375.at(numjets).Btag()) {
      nbtags++;
    }
  }

  if(total_ht > 375.0 && ht375.size() >=2) {
    //now check the jet quality criteria. if ok, fill esum struct
    if(ht375.at(0).Pt() > 100.0 && fabs(ht375.at(0).Eta()) < 2.4 && ht375.at(1).Pt() > 100.0) {
      esums.pass_quality_cuts = true;
      esums.njets = ht375.size();
      esums.nbtags = nbtags;
      esums.total_ht = total_ht;
      esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
      esums.etvec = etvec;
      esums.pxvec = pxvec;
      esums.pyvec = pyvec;
    } else {
      esums.pass_quality_cuts = false;
    }    
  } else {

    //now recalculate things with a lower jet threshold > 43 GeV
    total_ht = 0.0;
    mht_x = 0.0;
    mht_y = 0.0;
    njets = 0;
    nbtags = 0;
    etvec.clear();
    pxvec.clear();
    pyvec.clear();

    for(unsigned int numjets = 0; numjets < ht325.size(); numjets++) {
      total_ht += ht325.at(numjets).Pt();
      mht_x += ht325.at(numjets).Px();
      mht_y += ht325.at(numjets).Py();    
      etvec.push_back(ht325.at(numjets).Pt());
      pxvec.push_back(ht325.at(numjets).Px());
      pyvec.push_back(ht325.at(numjets).Py());
      if(ht325.at(numjets).Btag()) {
        nbtags++;
      }
    }

    //now check the ht
    if(total_ht > 325.0 && total_ht < 375.0 && ht325.size() >=2) {

      //now check the jet quality criteria. if ok, fill esum struct
      if(ht325.at(0).Pt() > 86.0 && fabs(ht325.at(0).Eta()) < 2.4 && ht325.at(1).Pt() > 86.0) {
        esums.pass_quality_cuts = true;
        esums.njets = ht325.size();
        esums.nbtags = nbtags;
        esums.total_ht = total_ht;
        esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
        esums.etvec = etvec;
        esums.pxvec = pxvec;
        esums.pyvec = pyvec;
      } else {
        esums.pass_quality_cuts = false;
      }
      
    } else if(total_ht > 375.0) {
      esums.pass_quality_cuts = false;
    } else {

      //now recalculate things with yet another lower jet threshold > 37 GeV
      total_ht = 0.0;
      mht_x = 0.0;
      mht_y = 0.0;
      nbtags = 0;
      njets = 0;
      etvec.clear();
      pxvec.clear();
      pyvec.clear();
      
      for(unsigned int numjets = 0; numjets < ht275.size(); numjets++) {
        total_ht += ht275.at(numjets).Pt();
        mht_x += ht275.at(numjets).Px();
        mht_y += ht275.at(numjets).Py();    
        etvec.push_back(ht275.at(numjets).Pt());
        pxvec.push_back(ht275.at(numjets).Px());
        pyvec.push_back(ht275.at(numjets).Py());
        if(ht275.at(numjets).Btag()) {
          nbtags++;
        }
      }
      
      //now check the ht
      //if(total_ht > 275.0 && total_ht < 325.0) { //to be able to see full distribution of HT in histograms
      if(total_ht < 325.0 && ht275.size() >=2) {

        //now check the jet quality criteria. if ok, fill esum struct
        if(ht275.at(0).Pt() > 74.0 && fabs(ht275.at(0).Eta()) < 2.4 && ht275.at(1).Pt() > 74.0) {
          esums.pass_quality_cuts = true;
          esums.njets = ht275.size();
          esums.nbtags = nbtags;
          esums.total_ht = total_ht;
          esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
          esums.etvec = etvec;
          esums.pxvec = pxvec;
          esums.pyvec = pyvec;
        } else {
          esums.pass_quality_cuts = false;
        }
	
      } else {
        esums.pass_quality_cuts = false;
      }

    }

  }

  return esums;

}


energy_sums make_energy_sums_20(const std::vector<jjet> & ht200, const std::vector<jjet> & ht275, const std::vector<jjet> & ht325, const std::vector<jjet> & ht375) {

  energy_sums esums;
  //initialise the energy sum quantities
  double total_ht = 0.0;
  double mht_x = 0.0;
  double mht_y = 0.0;
  unsigned int njets = 0;
  unsigned int nbtags = 0;
  std::vector<double> etvec;
  std::vector<double> pxvec;
  std::vector<double> pyvec;

  //first, loop over all jets based on threshold of >50 GeV and calculate HT
  for(unsigned int numjets = 0; numjets < ht375.size(); numjets++) {
    total_ht += ht375.at(numjets).Pt();
    mht_x += ht375.at(numjets).Px();
    mht_y += ht375.at(numjets).Py();    
    etvec.push_back(ht375.at(numjets).Pt());
    pxvec.push_back(ht375.at(numjets).Px());
    pyvec.push_back(ht375.at(numjets).Py());
    if(ht375.at(numjets).Btag()) {
      nbtags++;
    }
  }

  if(total_ht > 375.0 && ht375.size() >=2) {
    //now check the jet quality criteria. if ok, fill esum struct
    if(ht375.at(0).Pt() > 100.0 && fabs(ht375.at(0).Eta()) < 2.4 && ht375.at(1).Pt() > 100.0) {
      esums.pass_quality_cuts = true;
      esums.njets = ht375.size();
      esums.nbtags = nbtags;
      esums.total_ht = total_ht;
      esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
      esums.etvec = etvec;
      esums.pxvec = pxvec;
      esums.pyvec = pyvec;
    } else {
      esums.pass_quality_cuts = false;
    }    
  } else {

    //now recalculate things with a lower jet threshold > 43 GeV
    total_ht = 0.0;
    mht_x = 0.0;
    mht_y = 0.0;
    njets = 0;
    nbtags = 0;
    etvec.clear();
    pxvec.clear();
    pyvec.clear();

    for(unsigned int numjets = 0; numjets < ht325.size(); numjets++) {
      total_ht += ht325.at(numjets).Pt();
      mht_x += ht325.at(numjets).Px();
      mht_y += ht325.at(numjets).Py();    
      etvec.push_back(ht325.at(numjets).Pt());
      pxvec.push_back(ht325.at(numjets).Px());
      pyvec.push_back(ht325.at(numjets).Py());
      if(ht325.at(numjets).Btag()) {
        nbtags++;
      }
    }

    //now check the ht
    if(total_ht > 325.0 && total_ht < 375.0 && ht325.size() >=2) {

      //now check the jet quality criteria. if ok, fill esum struct
      if(ht325.at(0).Pt() > 100.0 && fabs(ht325.at(0).Eta()) < 2.4 && ht325.at(1).Pt() > 100.0) {
        esums.pass_quality_cuts = true;
        esums.njets = ht325.size();
        esums.nbtags = nbtags;
        esums.total_ht = total_ht;
        esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
        esums.etvec = etvec;
        esums.pxvec = pxvec;
        esums.pyvec = pyvec;
      } else {
        esums.pass_quality_cuts = false;
      }
      
    } else if(total_ht > 375.0) {
      esums.pass_quality_cuts = false;
    } else {

      //now recalculate things with yet another lower jet threshold > 37 GeV
      total_ht = 0.0;
      mht_x = 0.0;
      mht_y = 0.0;
      nbtags = 0;
      njets = 0;
      etvec.clear();
      pxvec.clear();
      pyvec.clear();
      
      for(unsigned int numjets = 0; numjets < ht275.size(); numjets++) {
        total_ht += ht275.at(numjets).Pt();
        mht_x += ht275.at(numjets).Px();
        mht_y += ht275.at(numjets).Py();    
        etvec.push_back(ht275.at(numjets).Pt());
        pxvec.push_back(ht275.at(numjets).Px());
        pyvec.push_back(ht275.at(numjets).Py());
        if(ht275.at(numjets).Btag()) {
          nbtags++;
        }
      }
      
      //now check the ht
      //if(total_ht > 275.0 && total_ht < 325.0)  //to be able to see full distribution of HT in histograms
      if(total_ht < 325.0 && ht275.size() >=2) {

        //now check the jet quality criteria. if ok, fill esum struct
        if(ht275.at(0).Pt() > 100.0 && fabs(ht275.at(0).Eta()) < 2.4 && ht275.at(1).Pt() > 100.0) {
          esums.pass_quality_cuts = true;
          esums.njets = ht275.size();
          esums.nbtags = nbtags;
          esums.total_ht = total_ht;
          esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
          esums.etvec = etvec;
          esums.pxvec = pxvec;
          esums.pyvec = pyvec;
        } else {
          esums.pass_quality_cuts = false;
        }
      }
      else if(total_ht > 325.0) {
        esums.pass_quality_cuts = false;
       
      } else {


        //now recalculate things with yet another lower jet threshold > 37 GeV
        total_ht = 0.0;
        mht_x = 0.0;
        mht_y = 0.0;
        nbtags = 0;
        njets = 0;
        etvec.clear();
        pxvec.clear();
        pyvec.clear();
      
        for(unsigned int numjets = 0; numjets < ht200.size(); numjets++) {
          total_ht += ht200.at(numjets).Pt();
          mht_x += ht200.at(numjets).Px();
          mht_y += ht200.at(numjets).Py();    
          etvec.push_back(ht200.at(numjets).Pt());
          pxvec.push_back(ht200.at(numjets).Px());
          pyvec.push_back(ht200.at(numjets).Py());
          if(ht200.at(numjets).Btag()) {
            nbtags++;
          }
        }
      
        //now check the ht
        //if(total_ht > 200.0 && total_ht < 275.0)  //to be able to see full distribution of HT in histograms
        if(total_ht < 275.0 && ht200.size() >=2) {

          //now check the jet quality criteria. if ok, fill esum struct
          if(ht200.at(0).Pt() > 100.0 && fabs(ht200.at(0).Eta()) < 2.4 && ht200.at(1).Pt() > 100.0) {
            esums.pass_quality_cuts = true;
            esums.njets = ht200.size();
            esums.nbtags = nbtags;
            esums.total_ht = total_ht;
            esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
            esums.etvec = etvec;
            esums.pxvec = pxvec;
            esums.pyvec = pyvec;
          } else {
            esums.pass_quality_cuts = false;
          }
	
        } else {
          esums.pass_quality_cuts = false;
        }


      }

    

    }


  }
  return esums;

}





energy_sums make_energy_sums(const std::vector<jjet> & ht40, const std::vector<jjet> & ht100) {
  energy_sums esums;
  //initialise the energy sum quantities
  double total_ht = 0.0;
  double mht_x = 0.0;
  double mht_y = 0.0;
  unsigned int njets = 0;
  unsigned int nbtags = 0;
  std::vector<double> etvec;
  std::vector<double> pxvec;
  std::vector<double> pyvec;

  //first, loop over all jets based on threshold of >50 GeV and calculate HT
  for(unsigned int numjets = 0; numjets < ht100.size(); numjets++) {
    total_ht += ht100.at(numjets).Pt();
    mht_x += ht100.at(numjets).Px();
    mht_y += ht100.at(numjets).Py();    
    etvec.push_back(ht100.at(numjets).Pt());
    pxvec.push_back(ht100.at(numjets).Px());
    pyvec.push_back(ht100.at(numjets).Py());
    if(ht100.at(numjets).Btag()) {
      nbtags++;
    }
  }

  if(total_ht > 200.0 && ht100.size() >=2) {
    //now check the jet quality criteria. if ok, fill esum struct
    if(ht100.at(0).Pt() > 100.0 && fabs(ht100.at(0).Eta()) < 2.4 && ht100.at(1).Pt() > 100.0) {
      esums.pass_quality_cuts = true;
      esums.njets = ht100.size();
      esums.nbtags = nbtags;
      esums.total_ht = total_ht;
      esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
      esums.etvec = etvec;
      esums.pxvec = pxvec;
      esums.pyvec = pyvec;
    } else {
      esums.pass_quality_cuts = false;
    }    
  } 
//   else {
//     //now recalculate things with a lower jet threshold > 43 GeV
//     total_ht = 0.0;
//     mht_x = 0.0;
//     mht_y = 0.0;
//     njets = 0;
//     nbtags = 0;
//     etvec.clear();
//     pxvec.clear();
//     pyvec.clear();
    
//     for(unsigned int numjets = 0; numjets < ht40.size(); numjets++) {
//       total_ht += ht40.at(numjets).Pt();
//       mht_x += ht40.at(numjets).Px();
//       mht_y += ht40.at(numjets).Py();    
//       etvec.push_back(ht40.at(numjets).Pt());
//       pxvec.push_back(ht40.at(numjets).Px());
//       pyvec.push_back(ht40.at(numjets).Py());
//       if(ht40.at(numjets).Btag()) {
//         nbtags++;
//       }
//     }
    
//     //now check the ht
//     if(total_ht > 200 && ht40.size() >=2) {
//       //now Check jet quality criteria. if ok, fill esum struct
//       if(ht40.at(0).Pt() > 100 && ht100.at(0).Pt() > 100 && ht40.at(0).Pt()<fabs(40.at(0).Eta()) < 2.4 ) {
//         esums.pass_quality_cuts = true;
//         esums.njets = ht40.size();
//         esums.nbtags = nbtags;
//         esums.total_ht = total_ht;
//         esums.total_mht = TMath::Sqrt((mht_x * mht_x) +(mht_y * mht_y));
//         esums.etvec = etvec;
//         esums.pxvec = pxvec;
//         esums.pyvec = pyvec;
//       } 
//       else {
//         esums.pass_quality_cuts = false;
//       }
//     } 
//   }



  return esums;  
 
}
