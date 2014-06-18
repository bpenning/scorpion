#include "jad_lepton_functions.hh"
#include "TMath.h"
    //Normal Jet Collection (eta and pt cuts)
    std::vector <jlepton> goodleptons(const std::vector <jlepton> & lepton_collection,
            const double &ptcut, const double & etacut, const float & etaw1, const float & etaw2)
    {
	std::vector<jlepton> mleptonvec;
	for (int ii=0; ii<lepton_collection.size();ii++)
	{
	    if(!(lepton_collection[ii].IsolFlag()) || lepton_collection[ii].Pt() < ptcut || fabs(lepton_collection[ii].Eta()) > etacut || (fabs(lepton_collection[ii].Eta()) > etaw1 && fabs(lepton_collection[ii].Eta()) < etaw2)) continue;
	    mleptonvec.push_back(lepton_collection[ii]);
	}
	return mleptonvec;
    }

    std::vector <jlepton> goodleptonsNoIso(const std::vector <jlepton> & lepton_collection,
            const double &ptcut, const double & etacut, const float & etaw1, const float & etaw2)
    {
	std::vector<jlepton> mleptonvec;
	for (int ii=0; ii<lepton_collection.size();ii++)
	{
	    if( lepton_collection[ii].Pt() < ptcut || fabs(lepton_collection[ii].Eta()) > etacut || (fabs(lepton_collection[ii].Eta()) > etaw1 && fabs(lepton_collection[ii].Eta()) < etaw2)) continue;
	    mleptonvec.push_back(lepton_collection[ii]);
	}
	return mleptonvec;
    }
    //Normal Jet Collection (eta and pt cuts)
    std::vector <jlepton> goodleptons(const std::vector <jlepton> & lepton_collection,double ptcut, double etacut)
    {
	std::vector<jlepton> mleptonvec;
	for (int ii=0; ii<lepton_collection.size();ii++)
	{
	    if(!(lepton_collection[ii].IsolFlag()) || lepton_collection[ii].Pt() < ptcut || fabs(lepton_collection[ii].Eta()) > etacut ) continue;
	    mleptonvec.push_back(lepton_collection[ii]);
	}
	return mleptonvec;
    }

    //Bad Jet Collection (pt cut normal but eta cut reversed)
    std::vector <jlepton> badleptons(std::vector <jlepton> & lepton_collection,double ptcut, double etacut)
    {
	std::vector<jlepton> mleptonvec;
	for (int ii=0; ii<lepton_collection.size();ii++)
	{
	    if(!(lepton_collection[ii].IsolFlag()) || lepton_collection[ii].Pt() < ptcut || fabs(lepton_collection[ii].Eta()) < etacut) continue;
	    mleptonvec.push_back(lepton_collection[ii]);
	}
	return mleptonvec;
    }
    //Generic leptons (with veto window for electrons)
std::vector<jlepton> leptonSkim(const std::vector<jlepton> & elec, 
        const std::vector<jlepton> & muon, const float & pte, 
        const float & etae, const float & ptm, const float & etam, 
        const float & etaw1,const float & etaw2) {
	std::vector<jlepton> leptons;

	for (int ii=0; ii<elec.size();ii++)
	{
	    if(elec[ii].Pt() < pte || !elec[ii].IsolFlag() || fabs(elec[ii].Eta()) > etae || (fabs(elec[ii].Eta()) > etaw1 && fabs(elec[ii].Eta()) < etaw2)) continue;
	    leptons.push_back(elec[ii]);
	}


	for (int ii=0; ii<muon.size();ii++)
	{
	    if(muon[ii].Pt()<ptm || !muon[ii].IsolFlag() || fabs(muon[ii].Eta()) > etam) continue;
	    leptons.push_back(muon[ii]);
	}

	std::sort(leptons.begin(), leptons.end(), std::greater<jlepton>()); //operators defined in the jlepton class
	return leptons;
}

    std::vector<jlepton> leptonSkim(const std::vector<jlepton> & elec, const std::vector<jlepton> & muon, 
            const double & pte, const double & etae, const double & ptm, const double & etam) {

	std::vector<jlepton> leptons;

	for (int ii=0; ii<elec.size();ii++)
	{
	    if(elec[ii].Pt() < pte || !elec[ii].IsolFlag() || fabs(elec[ii].Eta()) > etae ) continue;
	    leptons.push_back(elec[ii]);
	}


	for (int ii=0; ii<muon.size();ii++)
	{
	    if(muon[ii].Pt()<ptm || !muon[ii].IsolFlag() || fabs(muon[ii].Eta()) > etam) continue;
	    leptons.push_back(muon[ii]);
	}

	std::sort(leptons.begin(), leptons.end(), std::greater<jlepton>()); //operators defined in the jlepton class
	return leptons;
    }
/*
//Leptons for ATLAS5 (strip using DR with jets)
std::vector<jlepton> goodleptonsDR(std::vector<jlepton> & lepton_collection, float pt, float eta, double DRLmin, std::vector<jjet> &jet) {
    std::vector<jlepton> array;
    for (int ii=0; ii<lepton_collection.size();ii++)
    {
	double DRelemin = 10.0;
	for(unsigned int i=0; i <jet.size();i++) {

	    double DRele = TMath::Sqrt( (((lepton_collection[ii].Phi())-(jet[i].Phi()))*((lepton_collection[ii].Phi())-(jet[i].Phi()))) + (((lepton_collection[ii].Eta())-(jet[i].Eta()))*((lepton_collection[ii].Eta())-(jet[i].Eta()))) );
	    if(DRelemin > DRele) {
		DRelemin = DRele;
	    }

	}

	if (DRelemin > DRLmin) { //if the lepton_collection is away from all the jets
	    //if(lepton_collection[ii].PT<pt || !lepton_collection[ii].IsolFlag || fabs(lepton_collection[ii].Eta) > eta) continue;
	    if(lepton_collection[ii].Pt() > pt && lepton_collection[ii].IsolFlag() && fabs(lepton_collection[ii].Eta()) < eta) {
		array.push_back(lepton_collection[ii]);
	    }
	}
    }

    return array;
}
*/

std::vector<std::pair<jlepton,jlepton> > get_ossf_pairs(std::vector<jlepton> leptons){
     std::vector<std::pair<jlepton,jlepton> > ossf_pairs;
     std::vector<jlepton>::iterator l1_it;
     std::vector<jlepton>::iterator l2_it;
     std::vector<std::vector<jlepton>::const_iterator> paired_leptons;
     for (l1_it=leptons.begin();l1_it!=leptons.end();l1_it++){
        for (l2_it=l1_it+1;l2_it!=leptons.end();l2_it++){
            //if opposite sign same flavour
            if (l1_it->Charge()==-l2_it->Charge() && l1_it->Flavour()==l2_it->Flavour()){
                bool isPaired=false;
                std::vector<std::vector<jlepton>::const_iterator>::iterator lp_it;
                for (lp_it=paired_leptons.begin();lp_it!=paired_leptons.end();lp_it++){
                    if (*lp_it==l1_it || *lp_it==l2_it){
                        isPaired=true;
                        break;
                    }
                }
                if (!isPaired){
                    std::pair<jlepton,jlepton> ossf_pair((*l1_it),(*l2_it));
                    ossf_pairs.push_back(ossf_pair);
                    paired_leptons.push_back(l1_it);
                    paired_leptons.push_back(l2_it);
                }
            }
        }
     }
     return ossf_pairs;
}

std::vector<jlepton> leptonChargeSkim(const std::vector<jlepton> & elec,
        const std::vector<jlepton> & muon, const float & pte, 
        const float & etae, const float & ptm, const float & etam, 
        const float & etaw1, const float & etaw2, int charge){
    if (!(charge==1 || charge == -1))
        std::cout << "ERROR: LEPTON IS ASKED FOR NON-UNIT CHARGE" << std::endl;
	std::vector<jlepton> leptons;
	for (int ii=0; ii<elec.size();ii++)	{
	    if(elec[ii].Pt() < pte || !elec[ii].IsolFlag() || 
                fabs(elec[ii].Eta()) > etae || (fabs(elec[ii].Eta()) > etaw1 &&
                fabs(elec[ii].Eta()) < etaw2 ) || elec[ii].Charge() != charge) 
            continue;
	    leptons.push_back(elec[ii]);
	}

	for (int ii=0; ii<muon.size();ii++)	{
	    if(muon[ii].Pt()<ptm || !muon[ii].IsolFlag() || fabs(muon[ii].Eta()) > etam
                || muon[ii].Charge() != charge) 
            continue;
	    leptons.push_back(muon[ii]);
	}

	std::sort(leptons.begin(), leptons.end(), std::greater<jlepton>()); //operators defined in the jlepton class
	return leptons;
}
