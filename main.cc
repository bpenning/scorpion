#include "jad_jet_class.hh"
#include "jad_lepton_class.hh"
#include "jad_particle_class.hh"
#include "jBlockClasses.hh"
#include "TSimpleArray.hh"
#include "TreeReader.hh"
#include <iostream>
#include <string>
#include <TRandom1.h>
#include "TFile.h"
#include "TChain.h"
#include <TApplication.h>
#include <TH2.h>
#include <TH1.h>

TSimpleArray<TRootJet> SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta) {

    TIter itJet((TCollection*)JET);
    TRootJet *jet;
    itJet.Reset();
    TSimpleArray<TRootJet> array;
    while( (jet = (TRootJet*) itJet.Next()) ) {
	//check if any jet has Pt>50 and |eta|<3.0
	if(jet->PT > pt && fabs(jet->Eta) < eta) {
	    array.Add(jet);
	}
    }
    return array;
}

int main(int argc,char* argv[]){

    TApplication* rootapp = new TApplication("pt test",&argc, argv);
    TH1D * ptold = new TH1D("jettestold", ";P_{T} [GeV];Entries",500,-0.,2495.);
    std::string delphes_out = "delphes-output.root";
    TChain chainRec("Analysis");
    TChain chainGen("GEN");
    chainRec.Add((delphes_out).c_str());
    chainGen.Add((delphes_out).c_str());
    TreeReader oldtreereader(&chainRec);

    int numevents = oldtreereader.GetEntries();

    for(unsigned int event=0; event<numevents; event++) {
	oldtreereader.ReadEntry(event);	  
	TSimpleArray<TRootJet>     goodjets=SubArrayGoodJets(oldtreereader.Jet(),37.0,3.0); //check for jets which we should analyse
	for(unsigned int ii = 0; ii < goodjets.GetEntries();++ii)
	{ 
	    ptold->Fill(goodjets[ii]->PT);
	} 
    }
    ptold->Draw();
    rootapp->Run();
    return 0;

}
