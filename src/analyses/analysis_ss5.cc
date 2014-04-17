#include "analysis_ss5.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
SS::SS(const std::string & name, 
       const std::string & experiment,
       const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
SS::SS(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       //const std::vector<int> & datayields,
       const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

SS::SS(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       const std::vector<double> & bgpred,
       const std::vector<double> & bgpreduncert,
       const std::vector<int> & datayields,
       const std::string & fitmode,
       const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


SS::~SS() {}

std::vector<jparticle> SS::getparticles(const TClonesArray * GEN, const unsigned int & PIDfilter) {

  std::vector<jparticle> myparticles;

  TIter itGen((TCollection*)GEN);
  TRootC::GenParticle *gen;
  itGen.Reset();
  while( (gen = (TRootC::GenParticle*) itGen.Next()) ) {
    bool addparticle = true;
    if(PIDfilter != 0 && abs(gen->PID) != PIDfilter) {
      addparticle = false;
    }
    
    if(addparticle) {
      if(gen->Status == 2 && fabs(gen->Eta) < 2.5 && gen->PT > 2.0) {
	jparticle particle(gen->Px, gen->Py, gen->Pz, gen->E, gen->PID, gen->Status, gen->Charge);
	myparticles.push_back(particle);
      }
    }
  }

  return myparticles;

}

std::vector<jlepton> SS::getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam) {

  std::vector<jlepton> leptons;

  TIter itElec((TCollection*)ELEC);
  TRootElectron *elec;
  itElec.Reset();

  while( elec = ((TRootElectron*) itElec.Next()) ) {
    //if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
    //
    //double reliso = (elec->SumEt + elec->SumPt)/elec->PT;
    //std::cout << "reliso: " << reliso << std::endl;
    //reliso > riso
    if(elec->PT < pte || !elec->IsolFlag || fabs(elec->Eta) > etae || (fabs(elec->Eta) > 1.4 && fabs(elec->Eta) < 1.6)) continue;
    bool poscharge = true;
    if(elec->Charge < 0) { poscharge = false; }
    jlepton lepton(elec->Px, elec->Py, elec->Pz, elec->E, true, poscharge);
    leptons.push_back(lepton);
  }
  
  TIter itMuon((TCollection*)MUON);
  TRootMuon *muon;
  itMuon.Reset();
  
  while( muon = ((TRootMuon*) itMuon.Next()) ) {
    //double reliso = (muon->SumEt + muon->SumPt)/muon->PT;
    //if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
    //reliso > riso
    if(muon->PT<ptm || !muon->IsolFlag || fabs(muon->Eta) > etam) continue;
    bool poscharge = true;
    if(muon->Charge < 0) { 
      poscharge = false; 
    }
    jlepton lepton(muon->Px, muon->Py, muon->Pz, muon->E, false, poscharge);
    leptons.push_back(lepton);
  }

  std::sort(leptons.begin(), leptons.end(), std::greater<jlepton>()); //operators defined in the jlepton class
 
  return leptons;
}

std::vector<jjet> SS::getjets(const TClonesArray *JET, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim) {

  //this method drops jets with DR<0.4 with any leptons that pass the selection

  std::vector<jjet> jets;

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    //check if any jet has Pt>50 and |eta|<3.0
    if(jet->PT > pt && fabs(jet->Eta) < eta) {
      jjet myjet(jet->Px, jet->Py, jet->Pz, jet->E, jet->Btag);
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

double SS::getht(const std::vector<jjet> & jets) {
  double HT=0.0;
  for(unsigned int i=0;i<jets.size();i++) {
    HT += jets[i].Pt();
  }
  return HT;
}

double SS::getmet(const TClonesArray *ETMISS) {
  
  std::vector<double> etmvec;

  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();

  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
    etmvec.push_back(etm->ET);
  }

  if(etmvec.size() == 1) {
    return etmvec[0];
  } else {
    std::cout << "more than one MET vector defined in this event!" << std::endl;

  }
  return -1.0;

}

unsigned int SS::getnbtags(const std::vector<jjet> & jets) {

  unsigned int numbtags = 0;

  for(unsigned int i=0; i<jets.size(); i++) {
    if(jets[i].Btag()) {
      numbtags++;
    }
  }
  
  return numbtags;

}

void SS::initHistos() {
  andir->cd();
  ptleadinglep = new TH1D("ptleadinglep",";leading p_{T} [GeV];",100, -5.0, 995.0);
  ptchargeonevschargetwo = new TH2D("ptchargeonevschargetwo",";charge leading lep;charge sub-leading lep",2, -0.5, 1.5, 2, -0.5, 1.5);
  njethist = new TH1D("njethist",";N_{Jets};",20, -0.5, 19.5);
  hthist = new TH1D("hthist",";H_{T} [GeV];",250, -5.0, 2495.0);
  methist = new TH1D("methist",";E_{T}^{miss} [GeV];", 100, -5.0, 995.0);
  leadingjetpt = new TH1D("leadingjetpt",";leading-jet p_{T} [GeV];",100, -5.0, 995.0);
  numbjets = new TH1D("numbjets",";#b-tags;", 10, -0.5, 9.5);
  jetptc = new TH1D("jetptc",";jet p_{T} [GeV];", 100, -5.0, 995.0);
  bjetptc = new TH1D("bjetptc",";jet p_{T} [GeV];",100, -5.0, 995.0);
  jetptb = new TH1D("jetptb",";jet p_{T} [GeV];", 100, -5.0, 995.0);
  bjetptb = new TH1D("bjetptb",";jet p_{T} [GeV];",100, -5.0, 995.0);
  jetptm = new TH1D("jetptm",";jet p_{T} [GeV];", 100, -5.0, 995.0);
  bjetptm = new TH1D("bjetptm",";jet p_{T} [GeV];",100, -5.0, 995.0);
  sfinvmass = new TH1D("sfinvmass",";M_{SF} [GeV];", 1000, -0.5, 999.5);
  drjetcombo = new TH1D("drjetcombo",";#DeltaR;",500, -0.005, 4.995);
}

void SS::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //get generator particle information
  std::vector<jparticle> cquarks = getparticles(gentreereader.GenParticles(), 4);
  std::vector<jparticle> bquarks = getparticles(gentreereader.GenParticles(), 5);

  std::vector<jparticle> quarks(bquarks);
  quarks.insert(quarks.end(), cquarks.begin(), cquarks.end());

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> leptons = getleptons(treereader.Elec(), treereader.Muon(), 10.0, 2.4, 10.0, 2.4);
  std::vector<jjet> goodjets = getjets(treereader.Jet(), 40.0, 2.5, leptons, 0.4);

  //do some b quark matching to reco-jets
  std::vector<pair_info> pairs = make_pairs(quarks, goodjets);

  //std::cout << "new event: " << std::endl;
  //for(unsigned int i=0; i<pairs.size(); i++) {
  //  pairs[i].Print();
  //}

  std::vector<int> bmatched = analyse_pairs(pairs, goodjets.size(), 0.5);

  //for(unsigned int i=0; i<bmatched.size(); i++) {
  //  std::cout << "reco jet with index " << i << " is matched to flavour: " << bmatched[i] << std::endl;
  //}
  
  //plot dr of each jet with every other jet
  for(unsigned int i=0; i<goodjets.size(); i++) {
    for(unsigned int j=0; j<goodjets.size(); j++) {
      if(j > i) { //prevent jet matching with itself
	drjetcombo->Fill(goodjets[i].DeltaR(goodjets[j]), weight);
      }
    }
  }


  for(unsigned int i=0;i<goodjets.size();i++) {    
    if(abs(bmatched[i]) == 4) {
      jetptc->Fill(goodjets[i].Pt(), weight);
      if(goodjets[i].Btag()) {
	bjetptc->Fill(goodjets[i].Pt(), weight);
      }
    } else if(abs(bmatched[i]) == 5) {
      jetptb->Fill(goodjets[i].Pt(), weight);
      if(goodjets[i].Btag()) {
	bjetptb->Fill(goodjets[i].Pt(), weight);
      }
    } else {
      jetptm->Fill(goodjets[i].Pt(), weight);
      if(goodjets[i].Btag()) {
	bjetptm->Fill(goodjets[i].Pt(), weight);
      }
    }
  }
  

  if(goodjets.size() >= 2 ) {
    double HT = 1.05 * getht(goodjets); //5% scale correction
    double MET = getmet(treereader.ETMis());
    njethist->Fill(goodjets.size(), weight);
    leadingjetpt->Fill(goodjets[0].Pt(), weight);
    
    //check for two or three leptons that pass selection - otherwise we ignore these events
    if(leptons.size() == 2 || leptons.size() == 3) { 
      ptleadinglep->Fill(leptons[0].Pt(), weight);
      ptchargeonevschargetwo->Fill(leptons[0].Charge(), leptons[1].Charge(), weight);
      if(leptons[0].Pt() > 20.0) { //check leading lepton has at least 20 GeV
	if(leptons[0].Charge() == leptons[1].Charge()) { //check first two leptons have the same charge
	  double invmass = (leptons[0] + leptons[1]).M();
	  sfinvmass->Fill(invmass, weight);
	  if(invmass>8.0) { //if invmass > 8, we are ok
	    bool passselection = true;
	    if(leptons.size() == 3) { //if there are 3 leptons, need to check for Z-hypothesis if opposite charge
	      if(leptons[2].Charge() != leptons[0].Charge()) { //check for opposite charge i.e. maybe consistent with Z
		//check for Z-veto
		double invmass1 = (leptons[0] + leptons[2]).M();
		bool oppflav1 = false;
		if(leptons[0].Flavour() != leptons[2].Flavour()) {
		  oppflav1 = true;
		}
		double invmass2 = (leptons[1] + leptons[2]).M();
		bool oppflav2 = false;
		if(leptons[1].Flavour() != leptons[2].Flavour()) {
		  oppflav2 = true;
		}
		if((invmass1 < 106.0 && invmass1 > 76.0 && oppflav1) || (invmass2 < 106.0 && invmass2 > 76.0 && oppflav2)) {
		  passselection = false; //i.e. we have a Z candidate most likely, veto event...
		}
	      }
	    }
	    if(passselection) {
	      //we have a same-sign pair
	      hthist->Fill(HT, weight);
	      methist->Fill(MET, weight);
	      numbjets->Fill(getnbtags(goodjets), weight);
	      //for(unsigned int i=0;i<goodjets.size();i++) {
	      //  jetpt->Fill(goodjets[i].Pt(), weight);
	      //  if(goodjets[i].Btag()) {
	      //    bjetpt->Fill(goodjets[i].Pt(), weight);
	      //  }
	      //}
	      if(MET > 120.0 && HT> 450.0) {
		//our SR for the CMSSM
		if(leptons[0].Flavour() == leptons[1].Flavour()) {
		  if(leptons[0].Flavour() == "electron") {
		    //SR1: ee
		    mSigPred.at(0)+=weight;
		  } else {
		    //SR2: mumu
		    mSigPred.at(1)+=weight;
		  }
		} else {
		  //SR3: either emu or mue:
		  mSigPred.at(2)+=weight;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return;

}
