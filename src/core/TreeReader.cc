#include <iostream>
#include <iomanip>

#include "TBranchElement.h"
#include "TClonesArray.h"
#include "TChain.h"

#include "TreeReader.hh"

TreeReader::TreeReader(TTree *tree) :
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

TreeReader::~TreeReader()
{
  TBranchMap::iterator it_map;

  for(it_map = fBranchMap.begin(); it_map != fBranchMap.end(); ++it_map)
  {
    delete it_map->second.second;
  }

  //std::cout << "destructor called" << std::endl;

}

Bool_t TreeReader::ReadEntry(Long64_t entry)
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

TClonesArray *TreeReader::UseBranch(const TString &branchName)
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

void TreeReader::Init()
{
  TBranchMap::iterator it_map;

  for(it_map = fBranchMap.begin(); it_map != fBranchMap.end(); ++it_map)
  {
    fChain->SetBranchAddress(it_map->first, &(it_map->second.second));
  }

  Notify();

  fIsInitDone = kTRUE;
}

void TreeReader::Notify()
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

const TClonesArray * TreeReader::Jet() const {
  return JET;
}

const TClonesArray * TreeReader::Elec() const {
  return ELEC;
}

const TClonesArray * TreeReader::Muon() const {
  return MUON;
}

const TClonesArray * TreeReader::ETMis() const {
  return ETMIS;
}

const TClonesArray * TreeReader::GenParticles() const {
  return GENPARTICLE;
}

const TClonesArray * TreeReader::GenEvent() const {
  return GENEVENT;
}
