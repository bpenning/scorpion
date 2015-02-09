#include "jad_jet_functions.hh"
//Normal Jet Collection (eta and pt cuts)
std::vector <jjet> goodjetsSkim(const std::vector<jjet> & jet_collection,double ptcut, double etacut)
{
  std::vector<jjet> mjetvec;
  for (int ii=0; ii<jet_collection.size();ii++)
    {
      if(jet_collection[ii].Pt() < ptcut || fabs(jet_collection[ii].Eta()) > etacut) continue;
      mjetvec.push_back(jet_collection[ii]);
    }
  return mjetvec;
}

std::vector <jjet> goodbjetsSkim(const std::vector<jjet> & jet_collection,double ptcut, double etacut)
{
  std::vector<jjet> mjetvec;
  for (int ii=0; ii<jet_collection.size();ii++)
    {
      if(jet_collection[ii].Pt() < ptcut || fabs(jet_collection[ii].Eta()) > etacut || !jet_collection[ii].Btag()) continue;
      mjetvec.push_back(jet_collection[ii]);
    }
  return mjetvec;
}

//Bad Jet Collection (pt cut normal but eta cut reversed)
std::vector <jjet> badjetsSkim(const std::vector<jjet> & jet_collection,double ptcut, double etacut)
{
  std::vector<jjet> mjetvec;
  for (int ii=0; ii<jet_collection.size();ii++)
    {
      if(jet_collection[ii].Pt() < ptcut || fabs(jet_collection[ii].Eta()) < etacut) continue;
      mjetvec.push_back(jet_collection[ii]);
    }
  return mjetvec;
}
double getht(const std::vector<jjet> & jets) {
  double HT=0.0;
  for(unsigned int i=0;i<jets.size();i++) {
    HT += jets[i].Pt();
  }
  return HT;
}
//Jet collection with extra deltaR cut
std::vector<jjet> goodjetsSkimDRcut(const std::vector<jjet> & jet_collection, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim) {

  //this method drops jets with DR<0.4 with any leptons that pass the selection

  std::vector<jjet> jets;

  for (int ii=0; ii<jet_collection.size();ii++)
    {
      //check if any jet has Pt>50 and |eta|<3.0
      if(jet_collection[ii].Pt() > pt && fabs(jet_collection[ii].Eta()) < eta) {
        jjet myjet(jet_collection[ii].Px(), jet_collection[ii].Py(), jet_collection[ii].Pz(), jet_collection[ii].E(), jet_collection[ii].Btag(),jet_collection[ii].TauTag());
        bool veto = false;
        for(unsigned int j=0; j<lep.size(); j++) {
          if(myjet.DeltaR(lep[j]) < drlim) {
            veto = true;
          }
        }
        if(!veto) {
          jets.push_back(myjet);
        }
      }
    }

  std::sort(jets.begin(), jets.end(), std::greater<jjet>()); //operators defined in the jjets class

  return jets;

}

unsigned int getnbtags(const std::vector<jjet> & jets) {

  unsigned int numbtags = 0;

  for(unsigned int i=0; i<jets.size(); i++) {
    if(jets[i].Btag()) {
	    numbtags++;
    }
  }

  return numbtags;

}
/*
//Jets for ATLAS5 (strip using DR with leptons)
std::vector<jjet> goodjetsDR(std::vector<jjet> & jet_collection, float pt, float eta, double DRlim, std::vector<jlepton> &ele ) { 

std::vector<jjet> array;

double DRele, DRelemin;
for (int ii=0; ii<jet_collection.size();ii++)
{
DRele = 0.0;
DRelemin = 10.0;

for (unsigned int ee=0; ee < ele.size(); ee++) {
DRele = TMath::Sqrt( (((ele[ee].Phi())-(jet_collection[ii].Phi()))*((ele[ee].Phi()))-(jet_collection[ii].Phi())) + (((ele[ee].Eta())-(jet_collection[ii].Eta()))*((ele[ee].Eta())-(jet_collection[ii].Eta()))) ); 
if(DRelemin > DRele) {
DRelemin = DRele;
}
}

//check if any jet has Pt>50 and |eta|<3.0 and if the min is above threshold between electron and jet (otherwise it is a jet)
if(jet_collection[ii].Pt() > pt && fabs(jet_collection[ii].Eta()) < eta && DRelemin > DRlim) {
array.push_back(jet_collection[ii]); 
}
} 
return array; 
}
//Get deltaphi etmis vector and jet 
bool deltaphiptjet(std::vector<jjet> & jets, const double & px, const double & py, const double & deltaphi, const int & jetcuts, double & dphiret, double & phiret) {

double phietmiss = TMath::ATan2(py,px); //calculate etmiss phi - do it here as it doesn't change per jet!

double phijet = 100.0; //initialise to large number for sanity
double dphitemp = 100.0; //initialise to large number for sanity
double dphi = 100.0; //initialise to large number for sanity

for(unsigned int ii=0; ii < jets.size(); ii++) { //loop over all jet entries

phijet = TMath::ATan2(jets[ii].Py(), jets[ii].Px()); //calculate jet phi
dphitemp = fabs(phijet - phietmiss); //calculate the distance in phi between the jet total momentum

//for two vectors in a 2D plane, they cannot be more than pi away from each other
if(dphitemp > TMath::Pi()) {
dphitemp = fabs(dphitemp - (2.0 * TMath::Pi()));
}

//check for the smallest dphi between jet and etmiss
if (dphitemp < dphi) {   
dphi = dphitemp;
}

//we only consider the first jetcuts jet - but only if they exist so keep the for loop over the size
if (ii == (jetcuts-1)) {
break;
}
}

//populate some variables for plotting
phiret = phietmiss;
dphiret = dphi;

if (dphi <= deltaphi) { return false; } //if we find a jet with sep<DR relative to etmiss, return false
return true; //otherwise, we haven't found any, so return true

}*/
