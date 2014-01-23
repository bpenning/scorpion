#include "analysis_os5.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
OS::OS(const std::string & name, 
       const std::string & experiment,
       const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
OS::OS(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       //const std::vector<int> & datayields,
       const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

OS::OS(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       const std::vector<double> & bgpred,
       const std::vector<double> & bgpreduncert,
       const std::vector<int> & datayields,
       const std::string & fitmode,
       const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


OS::~OS() {}

std::vector<jlepton> OS::getleptons(const TClonesArray *ELEC, const TClonesArray *MUON, const float & pte, const float & etae, const float & ptm, const float & etam, bool poscharge) {

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
    if((elec->Charge > 0 && !poscharge) || (elec->Charge < 0 && poscharge)) continue;
    if(elec->PT < pte || !elec->IsolFlag || fabs(elec->Eta) > etae || (fabs(elec->Eta) > 1.4 && fabs(elec->Eta) < 1.6)) continue;
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
    if((muon->Charge > 0 && !poscharge) || (muon->Charge < 0 && poscharge)) continue;
    if(muon->PT<ptm || !muon->IsolFlag || fabs(muon->Eta) > etam) continue;
    jlepton lepton(muon->Px, muon->Py, muon->Pz, muon->E, false, poscharge);
    leptons.push_back(lepton);
  }

  std::sort(leptons.begin(), leptons.end(), std::greater<jlepton>()); //operators defined in the jlepton class
 
  return leptons;
}

std::vector<jjet> OS::getjets(const TClonesArray *JET, const float & pt, const float & eta) {

  std::vector<jjet> jets;

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    //check if any jet has Pt>50 and |eta|<3.0
    if(jet->PT > pt && fabs(jet->Eta) < eta) {
      jjet myjet(jet->Px, jet->Py, jet->Pz, jet->E, jet->Btag);
      jets.push_back(myjet);
    }
  }

  std::sort(jets.begin(), jets.end(), std::greater<jjet>()); //operators defined in the jjets class

  return jets;
}

std::vector<jjet> OS::getjets2(const TClonesArray *JET, const float & pt, const float & eta, const std::vector<jlepton> & lep, const double & drlim) {

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

double OS::getht(const std::vector<jjet> & jets) {
  double HT=0.0;
  for(unsigned int i=0;i<jets.size();i++) {
    HT += jets[i].Pt();
  }
  return HT;
}

double OS::getmet(const TClonesArray *ETMISS) {
  
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

bool OS::checkDRjetlep(const std::vector<jjet> & jets, const std::vector<jlepton> & leptons, const double & DRmin) {

  for(unsigned int i=0; i<jets.size(); i++) {
    for(unsigned int j=0; j<leptons.size(); j++) {
      double deltar = jets[i].DeltaR(leptons[j]);
      //deltarjetlep->Fill(deltar);
      if(deltar < DRmin) {
	return false;
      }

    }
  }

  return true; //i.e. no lepton found within DRmin of a jet

}

void OS::initHistos() {
  andir->cd();
  ptleadingleppos = new TH1D("ptleadingleppos",";leading p_{T}+ [GeV];",100, -5.0, 995.0);
  ptleadinglepneg = new TH1D("ptleadinglepneg",";leading p_{T}- [GeV];",100, -5.0, 995.0);
  ptleadinglepposvsneg = new TH2D("ptleadinglepposvsneg",";leading p_{T}+;leading p_{T}-",100, -5.0, 995.0, 100, -5.0, 995.0);
  njethist = new TH1D("njethist",";N_{Jets};",20, -0.5, 19.5);
  deltarjetlep = new TH1D("deltarjetlep",";#DeltaR(jet, lep);",100, -0.005, 0.995);
  hthist = new TH1D("hthist",";H_{T} [GeV];",250, -5.0, 2495.0);
  methist = new TH1D("methist",";E_{T}^{miss} [GeV];", 100, -5.0, 995.0);
  leadingjetpt = new TH1D("leadingjetpt",";leading-jet p_{T} [GeV];",100, -5.0, 995.0);
  htvsmetsr1sf = new TH2D("htvsmetsr1sf",";H_{T} [GeV];E_{T}^{miss} [GeV]",100, -5.0, 995.0, 100, -5.0, 995.0);
  htvsmetsr2sf = new TH2D("htvsmetsr2sf",";H_{T} [GeV];E_{T}^{miss} [GeV]",100, -5.0, 995.0, 100, -5.0, 995.0);
  htvsmetsr3sf = new TH2D("htvsmetsr3sf",";H_{T} [GeV];E_{T}^{miss} [GeV]",100, -5.0, 995.0, 100, -5.0, 995.0);
  htvsmetsr1of = new TH2D("htvsmetsr1of",";H_{T} [GeV];E_{T}^{miss} [GeV]",100, -5.0, 995.0, 100, -5.0, 995.0);
  htvsmetsr2of = new TH2D("htvsmetsr2of",";H_{T} [GeV];E_{T}^{miss} [GeV]",100, -5.0, 995.0, 100, -5.0, 995.0);
  htvsmetsr3of = new TH2D("htvsmetsr3of",";H_{T} [GeV];E_{T}^{miss} [GeV]",100, -5.0, 995.0, 100, -5.0, 995.0);
  sfinvmass = new TH1D("sfinvmass",";M_{SF} [GeV];", 1000, -0.5, 999.5);
}

void OS::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  //std::vector<jjet> goodjets = getjets(treereader.Jet(), 30.0, 3.0);
  std::vector<jlepton> poslep = getleptons(treereader.Elec(), treereader.Muon(), 10.0, 2.5, 10.0, 2.4, true);
  std::vector<jlepton> neglep = getleptons(treereader.Elec(), treereader.Muon(), 10.0, 2.5, 10.0, 2.4, false);

  std::vector<jlepton> leptons = poslep;
  leptons.insert(leptons.end(), neglep.begin(), neglep.end());


  std::vector<jjet> goodjets = getjets2(treereader.Jet(), 30.0, 3.0, leptons, 0.4);

  if(goodjets.size() >= 2 ) {
    double HT = 1.05 * getht(goodjets);
    double MET = getmet(treereader.ETMis());
    njethist->Fill(goodjets.size(), weight);
    leadingjetpt->Fill(goodjets[0].Pt(), weight);

    if(poslep.size() && neglep.size()) { //there has to be at least one lepton of each sign
      ptleadingleppos->Fill(poslep[0].Pt(), weight);
      ptleadinglepneg->Fill(neglep[0].Pt(), weight);
      ptleadinglepposvsneg->Fill(poslep[0].Pt(), neglep[0].Pt(), weight);
      if(poslep[0].Pt() > 10.0 && neglep[0].Pt() > 10.0) { //both of the leading leptons have to have at least 10 GeV
	if(poslep[0].Pt() > 20.0 || neglep[0].Pt() > 20.0) { //one of the leading leptons must have 20 GeV
	  //now we should check that each jet is separated by DR>0.4 from each lepton which passes the event selection
	  //if(checkDRjetlep(goodjets, poslep, 0.4) && checkDRjetlep(goodjets, neglep, 0.4)) { //hold off on this - not clear
	    hthist->Fill(HT, weight);
	    methist->Fill(MET, weight);
	    if(poslep[0].Flavour() == neglep[0].Flavour()) {
	      //SF bin
	      double invmass = (poslep[0] + neglep[0]).M(); //check its outside Z-window
	      sfinvmass->Fill(invmass, weight);
	      if((invmass > 12.0 && invmass < 76.0) || invmass > 106.0) {
		if(MET > 275.0) {
		  if(HT > 300.0 && HT < 600.0) {
		    //SR1SF
		    htvsmetsr1sf->Fill(HT, MET, weight);
		    mSigPred.at(0)+=weight;
		  } else if(HT > 600.0) {
		    //SR2SF
		    htvsmetsr2sf->Fill(HT, MET, weight);
		    mSigPred.at(1)+=weight;
		  }
		} else if(MET > 200.0 && MET < 275.0 && HT > 600.0) {
		  //SR3SF
		  htvsmetsr3sf->Fill(HT, MET, weight);
		  mSigPred.at(2)+=weight;
		}
	      }
	      
	    } else {
	      //OF bin
	      if(MET > 275.0) {
		if(HT > 300.0 && HT < 600.0) {
		  //SR1OF
		  htvsmetsr1of->Fill(HT, MET, weight);
		  mSigPred.at(3)+=weight;
		} else if(HT > 600.0) {
		  //SR2OF
		  htvsmetsr2of->Fill(HT, MET, weight);
		  mSigPred.at(4)+=weight;
		}
	      } else if(MET > 200.0 && MET < 275.0 && HT > 600.0) {
		//SR3OF
		htvsmetsr3of->Fill(HT, MET, weight);
		mSigPred.at(5)+=weight;
	      }	      
	    }
	    
	  //} //for the DR cut?
	  
	}
      }
    }
 
  }

  return;

}
