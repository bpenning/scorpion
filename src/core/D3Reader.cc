
/** \class D3Reader
 *
 *  Class simplifying access to ROOT tree branches
 *
 *  $Date: 2008-06-04 13:57:57 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "D3Reader.hh"
#include "DelphesClasses.h"

//#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TBranchElement.h"

#include <iostream>

using namespace std;

//------------------------------------------------------------------------------
D3Reader::D3Reader(){} 

D3Reader::D3Reader(TTree *tree) :
  fChain(tree), fCurrentTree(-1)
{
  std::string anstring = "Analysis";
    JET    = this->UseBranch("Jet");
//    TAUJET = this->UseBranch("TauJet");
    PHOTON = this->UseBranch("Photon");
    ELEC   = this->UseBranch("Electron");
    MUON   = this->UseBranch("Muon");
    ETMIS  = this->UseBranch("MissingET"); 
    CTRACK  = this->UseBranch("ChargedTracks"); 
    GENEVENT  = this->UseBranch("Event"); 
    SCALARHT  = this->UseBranch("ScalarHT");
//Add this info later as takes time

//    GENEVENT = this->UseBranch("Event");
    GENPARTICLE = this->UseBranch("Particle");

}
//TEST
D3Reader::D3Reader(TTree *tree, bool gen) :
    fChain(tree), fCurrentTree(-1)
{
    if (!gen)
    {
	JET    = this->UseBranch("Jet");
	//    TAUJET = this->UseBranch("TauJet");
	PHOTON = this->UseBranch("Photon");
	ELEC   = this->UseBranch("Electron");
	MUON   = this->UseBranch("Muon");
	ETMIS  = this->UseBranch("MissingET"); 
	CTRACK  = this->UseBranch("ChargedTracks"); 
	SCALARHT = this->UseBranch("ScalarHT");
    }
    else 
    {
	GENEVENT  = this->UseBranch("Event"); 
	//Add this info later as takes time

	//    GENEVENT = this->UseBranch("Event");
	GENPARTICLE = this->UseBranch("Particle");
    }
}

//------------------------------------------------------------------------------

D3Reader::~D3Reader()
{
    TBranchMap::iterator itBranchMap;

    for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap)
    {
	delete itBranchMap->second.second;
    }
}

//------------------------------------------------------------------------------

Bool_t D3Reader::ReadEntry(Long64_t entry)
{
    // Read contents of entry.
    if(!fChain) return kFALSE;

    Int_t treeEntry = fChain->LoadTree(entry);
    if(treeEntry < 0) return kFALSE;

    if(fChain->IsA() == TChain::Class())
    {
	TChain *chain = static_cast<TChain*>(fChain);
	if(chain->GetTreeNumber() != fCurrentTree)
	{
	    fCurrentTree = chain->GetTreeNumber();
	    Notify();
	}
    }

    TBranchMap::iterator itBranchMap;
    TBranch *branch;

    for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap)
    {
	branch = itBranchMap->second.first;
	if(branch)
	{
	    branch->GetEntry(treeEntry);
	}
    }

    return kTRUE;
}

//------------------------------------------------------------------------------

TClonesArray *D3Reader::UseBranch(const char *branchName)
{
    TClonesArray *array = 0;

    TBranchMap::iterator itBranchMap = fBranchMap.find(branchName);

    if(itBranchMap != fBranchMap.end())
    {
	cout << "** WARNING: branch '" << branchName << "' is already in use" << endl;
	array = itBranchMap->second.second;
    }
    else
    {
	TBranch *branch = fChain->GetBranch(branchName);
	if(branch)
	{
	    if(branch->IsA() == TBranchElement::Class())
	    {
		TBranchElement *element = static_cast<TBranchElement*>(branch);
		const char *className = element->GetClonesName();
		Int_t size = element->GetMaximum();
		TClass *cl = gROOT->GetClass(className);
		if(cl)
		{
		    array = new TClonesArray(cl, size);
		    array->SetName(branchName);
		    fBranchMap.insert(make_pair(branchName, make_pair(branch, array)));
		    branch->SetAddress(&array);
		}
	    }
	}
    }

    if(!array)
    {
	cout << "** WARNING: cannot access branch '" << branchName << "', return NULL pointer" << endl;
    }

    return array;
}

//------------------------------------------------------------------------------

Bool_t D3Reader::Notify()
{
    // Called when loading a new file.
    // Get branch pointers.
    if(!fChain) return kFALSE;

    TBranchMap::iterator itBranchMap;
    TBranch *branch;

    for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap)
    {
	branch = fChain->GetBranch(itBranchMap->first);
	if(branch)
	{
	    itBranchMap->second.first = branch;
	    branch->SetAddress(&(itBranchMap->second.second));
	}
	else
	{
	    cout << "** WARNING: cannot get branch '" << itBranchMap->first << "'" << endl;
	}
    }
    return kTRUE;
}
std::vector<jjet> D3Reader::GetJet() const {
    std::vector<jjet> jet_collection;

    for(int i = 0; i < JET->GetEntries(); ++i)
    {
	Jet *jet = (Jet*) JET->At(i);
	if (!jet->TauTag) jet_collection.push_back(jjet(jet->P4().Px(),jet->P4().Py(),jet->P4().Pz(),jet->P4().E(),jet->BTag,false));
    }
    std::sort(jet_collection.begin(), jet_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return jet_collection;
}
std::vector<jtrack> D3Reader::GetIsoChargedTrack() const {
    std::vector<jtrack> track_collection;
    if(CTRACK)
    {
	for(int i = 0; i < CTRACK->GetEntries(); ++i)
	{
	    Track *track = (Track*) CTRACK->At(i);
	    track_collection.push_back(jtrack(track->P4().Px(),track->P4().Py(),track->P4().Pz(),track->P4().E(),track->Charge,true));
	}
	std::sort(track_collection.begin(), track_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    }
    /*else
    {
	std::cerr<<"No charged tracks info - returning empty vector" << std::endl;
    }*/
    return track_collection;
}

std::vector<jlepton> D3Reader::GetElec() const {
    std::vector<jlepton> electron_collection; 

    for(int i = 0; i < ELEC->GetEntries(); ++i)
    {
	Electron *elec = (Electron*) ELEC->At(i);
	electron_collection.push_back(jlepton(elec->P4().Px(), elec->P4().Py(), elec->P4().Pz(), elec->P4().E(), true, elec->Charge,true)); //elec->IsolFlag));
    }

    std::sort(electron_collection.begin(), electron_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return electron_collection;
}

std::vector<jlepton> D3Reader::GetMuon() const {
    std::vector<jlepton> muon_collection; 

    for(int i = 0; i < MUON->GetEntries(); ++i)
    {
	Muon *muon = (Muon*) MUON->At(i);
	muon_collection.push_back(jlepton(muon->P4().Px(), muon->P4().Py(), muon->P4().Pz(), muon->P4().E(), false, muon->Charge,true));
    }

    std::sort(muon_collection.begin(), muon_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return muon_collection;
}
std::vector<jphoton> D3Reader::GetPhoton() const {
    std::vector<jphoton> photon_collection; 

    for(int i = 0; i < PHOTON->GetEntries(); ++i)
    {
	Photon *photon = (Photon*) PHOTON->At(i);
	photon_collection.push_back(jphoton(photon->P4().Px(), photon->P4().Py(), photon->P4().Pz(), photon->P4().E()));
    }

    std::sort(photon_collection.begin(), photon_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return photon_collection;
}

std::vector<jjet> D3Reader::GetMet() const {
    std::vector<jjet> etmiss_collection;
    for(int i = 0; i < ETMIS->GetEntries(); ++i)
    {
	MissingET *etm = (MissingET*) ETMIS->At(i);
	etmiss_collection.push_back(jjet(etm->P4().Px(),etm->P4().Py(),0.,etm->MET,false,false));
    }

    return etmiss_collection;
}

std::vector<jjet> D3Reader::GetTauJet() const {
    std::vector<jjet> jet_collection;

    for(int i = 0; i < JET->GetEntries(); ++i)
    {
	Jet *jet = (Jet*) JET->At(i);
	if (jet->TauTag) jet_collection.push_back(jjet(jet->P4().Px(),jet->P4().Py(),jet->P4().Pz(),jet->P4().E(),jet->BTag,true));
    }
    std::sort(jet_collection.begin(), jet_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return jet_collection;
}

std::vector<jparticle> D3Reader::GetGenParticle() const {

    std::vector<jparticle> particle_collection;

    for(int i = 0; i < GENPARTICLE->GetEntries(); ++i)
    {
	GenParticle *gen = (GenParticle*) GENPARTICLE->At(i);
	particle_collection.push_back(jparticle(gen->P4().Px(), gen->P4().Py(), gen->P4().Pz(), gen->P4().E(), gen->PID, gen->Status, gen->Charge));
    }
    std::sort(particle_collection.begin(), particle_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return particle_collection;
}

std::vector<double> D3Reader::GetScalarHT() const {
    std::vector<double> scalarht_collection;
    for(int i = 0; i < SCALARHT->GetEntries(); ++i)
    {
	ScalarHT *scalarht = (ScalarHT*) SCALARHT->At(i);
	scalarht_collection.push_back(scalarht->HT);
    }

    return scalarht_collection;
}

double D3Reader::GetWeight() const {
    if (GENEVENT->GetEntries()>0)
    {
	HepMCEvent *event = (HepMCEvent *)GENEVENT->At(0); 
	return event->Weight;
    } 
    return 1.0;
}

//   const TClonesArray * D2Reader::GenEvent() const {
//   return GENEVENT;
//   }

//------------------------------------------------------------------------------

