#include <iostream>
#include <iomanip>

#include "TBranchElement.h"
#include "TClonesArray.h"
#include "TChain.h"

#include "D2Reader.hh"
#include "Reader.hh"

D2Reader::D2Reader(TTree *tree) :
  fChain(tree), fCurrentTree(-1), fIsInitDone(kFALSE)
{

  //std::cout << "tree = " << tree->GetName() << std::endl;
  std::string anstring = "Analysis";
  std::string genstring = "GEN";
  //use the compare function as ->GetName() returns const char*

  if(!anstring.compare(tree->GetName())) {
    JET    = this->UseBranch("Jet");
    //TAUJET = this->UseBranch("TauJet");
    //PHOTON = this->UseBranch("Photon");
    ELEC   = this->UseBranch("Electron");
    MUON   = this->UseBranch("Muon");
    ETMIS  = this->UseBranch("ETmis"); 
  } else if(!genstring.compare(tree->GetName())) {
    GENEVENT = this->UseBranch("Event");
    GENPARTICLE = this->UseBranch("Particle");
  } else {
    std::cout << "Invalid tree name specified!" << std::endl;
  }

  
}

D2Reader::~D2Reader()
{
  TBranchMap::iterator it_map;

  for(it_map = fBranchMap.begin(); it_map != fBranchMap.end(); ++it_map)
  {
    delete it_map->second.second;
  }

  //std::cout << "destructor called" << std::endl;

}
Long64_t D2Reader::GetEntries() const 
{
return fChain ? static_cast<Long64_t>(fChain->GetEntries()) : 0; 
}

Bool_t D2Reader::ReadEntry(Long64_t entry)
{
  // Read contents of entry.
  if(!fChain) return kFALSE;

  if(!fIsInitDone) Init();

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

  TBranchMap::iterator it_map;
  TBranch *branch;

  for(it_map = fBranchMap.begin(); it_map != fBranchMap.end(); ++it_map)
  {
    branch = it_map->second.first;
    if(branch)
    {
      branch->GetEntry(treeEntry);
    }
  }

  return kTRUE;
}

TClonesArray *D2Reader::UseBranch(const TString &branchName)
{
  fIsInitDone = kFALSE;

  TClonesArray *array = 0;

  TBranchMap::iterator it_map = fBranchMap.find(branchName);

  if(it_map != fBranchMap.end())
  {
    array = it_map->second.second;
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
          fBranchMap.insert(std::make_pair(branchName, std::make_pair(branch, array)));
        }
      }
    }
  }

  if(!array)
  {
    std::cerr << std::left  << std::setw(35) <<"**   ERROR: cannot access branch"<<""
	      << std::left  << std::setw(15) << branchName                       <<""
	      << std::right << std::setw(19) <<" **"<<""<< std::endl;
  }

  return array;
}

void D2Reader::Init()
{
  TBranchMap::iterator it_map;

  for(it_map = fBranchMap.begin(); it_map != fBranchMap.end(); ++it_map)
  {
    fChain->SetBranchAddress(it_map->first, &(it_map->second.second));
  }

  Notify();

  fIsInitDone = kTRUE;
}

void D2Reader::Notify()
{
  // Called when loading a new file.
  // Get branch pointers.
  if(!fChain) return;

  TBranchMap::iterator it_map;

  for(it_map = fBranchMap.begin(); it_map != fBranchMap.end(); ++it_map)
  {
    it_map->second.first = fChain->GetBranch(it_map->first);
    if(!it_map->second.first)
    {
      std::cerr << std::left  << std::setw(35) <<"**   ERROR: cannot get branch"<<""
		<< std::left  << std::setw(15) << it_map->first                  <<""
		<< std::right << std::setw(19) <<" **"<<""<< std::endl;
    }
  }
}

std::vector<jjet> D2Reader::Jet() const {
    std::vector<jjet> jet_collection;

    TIter itJet(JET);
    TRootJet *jet;
    itJet.Reset();
    while( (jet = (TRootJet*) itJet.Next()) ) {
	jet_collection.push_back(jjet(jet->Px,jet->Py,jet->Pz,jet->E,jet->Btag));
    }
    std::sort(jet_collection.begin(), jet_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return jet_collection;
}

std::vector<jlepton> D2Reader::Elec() const {
    std::vector<jlepton> electron_collection; 

    TIter itElec((TCollection*)ELEC);
    TRootElectron *elec;
    itElec.Reset();

    while( elec = ((TRootElectron*) itElec.Next()) ) {

	electron_collection.push_back(jlepton(elec->Px, elec->Py, elec->Pz, elec->E, true, elec->Charge,elec->IsolFlag));
    }

    std::sort(electron_collection.begin(), electron_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return electron_collection;
}

std::vector<jlepton> D2Reader::Muon() const {
    std::vector<jlepton> muon_collection; 

    TIter itMuon((TCollection*)MUON);
    TRootMuon *muon;
    itMuon.Reset();

    while( muon = ((TRootMuon*) itMuon.Next()) ) {
	muon_collection.push_back(jlepton(muon->Px, muon->Py, muon->Pz, muon->E, false, muon->Charge,muon->IsolFlag));
    }

    std::sort(muon_collection.begin(), muon_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return muon_collection;
}

std::vector<jjet> D2Reader::ETMis() const {

    std::vector<jjet> etmiss_collection;
    TIter itEtMiss((TCollection*)ETMIS);
    TRootETmis *etm;
    itEtMiss.Reset();

    while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
	etmiss_collection.push_back(jjet(etm->Px,etm->Py,0.,etm->ET,false));
    }

    return etmiss_collection;
}
std::vector<jparticle> D2Reader::GenParticles() const {

    std::vector<jparticle> particle_collection;

    TIter itGen((TCollection*)GENPARTICLE);
    TRootC::GenParticle *gen;
    itGen.Reset();
    while( (gen = (TRootC::GenParticle*) itGen.Next()) ) {

	particle_collection.push_back(jparticle(gen->Px, gen->Py, gen->Pz, gen->E, gen->PID, gen->Status, gen->Charge));
    }
    std::sort(particle_collection.begin(), particle_collection.end(), std::greater<jobject>()); //operators defined in the jobject class
    return particle_collection;
}
/*
   const TClonesArray * D2Reader::GenEvent() const {
   return GENEVENT;
   }*/
