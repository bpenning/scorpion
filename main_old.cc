#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
#include "jad_particle_class.hh"
#include <iostream>
#include <string>
#include <TRandom1.h>
#include "TFile.h"
#include "TChain.h"
#include <TApplication.h>
#include <TH1.h>

#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "D3Reader.hh"
#include "DelphesClasses.h"
std::vector <jjet> goodjetsSkim(std::vector<jjet> jet_collection,double ptcut, double etacut)
{
    std::vector<jjet> mjetvec;
    for (int ii=0; ii<jet_collection.size();ii++)
    {
	if(jet_collection[ii].Pt() < ptcut || fabs(jet_collection[ii].Eta()) > etacut) continue;
	//std::cout << jet_collection[ii].Pt() << std::endl;
	mjetvec.push_back(jet_collection[ii]);
    }
    return mjetvec;
}

int main(int argc,char* argv[])
{

    //gSystem->Load("libDelphes");
    TApplication* rootapp = new TApplication("pt test",&argc, argv);
    TH1D * pthist = new TH1D("jettest", ";P_{T} [GeV];Entries",500,-2505.,2495.);
    std::string delphes_out = "DMoutput.root";
    TChain chainRec("Delphes");
    chainRec.Add((delphes_out).c_str());


    D3Reader treeReader(&chainRec);

    Long64_t numberOfEntries = treeReader.GetEntries();

    // Get pointers to branches used in this analysis
    //TClonesArray *branchJet = treeReader->UseBranch("Jet");

    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
	// Load selected branches with data from specified event
	treeReader.ReadEntry(entry);
    std::vector <jjet> jets = treeReader.GetJet();

	// If event contains at least 1 jet
	if(jets.size() > 0)
	{
	    // Take first jet
	    //Jet *jet = (Jet*) branchJet->At(0);

	    // Plot jet transverse momentum
	    pthist->Fill(jets[0].Pt());

	    // Print jet transverse momentum
	    //std::cout << jet->PT << std::endl;
	}
    }
    pthist->Draw();
    rootapp->Run();
    return 0;
}

