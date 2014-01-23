#include "analysis_lp5.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
LP::LP(const std::string & name, 
       const std::string & experiment,
       const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
LP::LP(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       //const std::vector<int> & datayields,
       const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

LP::LP(const std::string & name, 
       const std::string & experiment, 
       const unsigned int & numBins,
       const double & intlumi, 
       const std::vector<double> & bgpred,
       const std::vector<double> & bgpreduncert,
       const std::vector<int> & datayields,
       const std::string & fitmode,
       const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


LP::~LP() {}

TSimpleArray<TRootJet> LP::SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta) {

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

TSimpleArray<TRootElectron> LP::SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta, const float & riso) {
  TIter itElec((TCollection*)ELEC);
  TRootElectron *elec;
  itElec.Reset();
  TSimpleArray<TRootElectron> array;
  while( (elec = (TRootElectron*) itElec.Next()) )
    {
      //if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
      //
      //double reliso = (elec->SumEt + elec->SumPt)/elec->PT;
      //std::cout << "reliso: " << reliso << std::endl;
      //reliso > riso
      if(elec->PT < pt || !elec->IsolFlag || fabs(elec->Eta) > eta || (fabs(elec->Eta) > 1.4 && fabs(elec->Eta) < 1.6)) continue;
      array.Add(elec);
    }
  return array;
}

TSimpleArray<TRootMuon> LP::SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta, const float & riso) {
  TIter itMuon((TCollection*)MUON);
  TRootMuon *muon;
  itMuon.Reset();
  TSimpleArray<TRootMuon> array;
  while( (muon = (TRootMuon*) itMuon.Next()) )
    {
      //double reliso = (muon->SumEt + muon->SumPt)/muon->PT;
      //if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
      //reliso > riso
      if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
      array.Add(muon);
    }
  return array;
}

TSimpleArray<TRootETmis> LP::makeETM(const TClonesArray *ETMISS) {
  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();
  TSimpleArray<TRootETmis> array;
  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
    array.Add(etm);
  }
  return array;
}

void LP::initHistos() {
  andir->cd();
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  num_mu_veto = new TH1D("num_mu_veto", ";Muons with P_{T} > 15 GeV;Entries", 10, -0.5, 9.5);
  num_ele_veto = new TH1D("num_ele_veto", ";Electrons with P_{T} > 15 GeV;Entries", 10, -0.5, 9.5);
  hthist = new TH1D("hthist", ";H_{T} [GeV];Entries",250,-5.,2495.);
  lphist_el = new TH1D("lphist_el",";L_{P};Entries",500, -2.505, 2.495);
  sthist_el = new TH1D("sthist_el", ";S_{T} [GeV];Entries",200,-5., 1995.);
  mthist_el = new TH1D("mthist_el",";M_{T} [GeV];Entries",250, -5.,2495.);
  lphist_mu = new TH1D("lphist_mu",";L_{P};Entries",500, -2.505, 2.495);
  sthist_mu = new TH1D("sthist_mu", ";S_{T} [GeV];Entries",200,-5., 1995.);
  mthist_mu = new TH1D("mthist_mu",";M_{T} [GeV];Entries",250, -5.,2495.);
}

void LP::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  TSimpleArray<TRootElectron> ele=SubArrayEl(treereader.Elec(), 20.0, 2.4, 0.065); //the central isolated electrons, pt > PT_ELEC GeV
  TSimpleArray<TRootElectron> eleveto=SubArrayEl(treereader.Elec(), 15.0, 2.4, 0.065); //veto on second lepton with pt > 15 
  TSimpleArray<TRootMuon> mu=SubArrayMu(treereader.Muon(), 20.0, 2.1, 0.10); //the central isolated muons, pt > PT_MUON GeV
  TSimpleArray<TRootMuon> muveto=SubArrayMu(treereader.Muon(), 15.0, 2.1, 0.10); //veto on the existence of a second lepton with pt > 15
  TSimpleArray<TRootJet> goodjets=SubArrayGoodJets(treereader.Jet(), 40.0, 2.4); //check for jets which we should analyse
  TSimpleArray<TRootETmis> etmis=makeETM(treereader.ETMis()); //Missing transverse energy array

  //double calo_met = etmis[0]->ET;
  num_mu_veto->Fill(muveto.GetEntries(), weight);
  num_ele_veto->Fill(eleveto.GetEntries(), weight);


  if(goodjets.GetEntries() >=3 && mu.GetEntries() == 1 && eleveto.GetEntries() == 0 && muveto.GetEntries() == 1) { //should we also check muveto.GetEntries() == 1?
      
    leadingjetpt->Fill(goodjets[0]->PT, weight);
    double total_ht = 0.0;
    unsigned int numbtags = 0;
    for(unsigned int numjets = 0; numjets < goodjets.GetEntries(); numjets++) {
      total_ht += goodjets[numjets]->PT;
      if(goodjets[numjets]->Btag) { numbtags++; }
    }
      
    //btagrate->Fill(numbtags);

    hthist->Fill(total_ht, weight);

    if(total_ht > 500.0) {
	
      double ST = (mu[0]->PT + etmis[0]->ET);
      
      double LP = ((mu[0]->Px * etmis[0]->Px) + (mu[0]->Py * etmis[0]->Py)) / (etmis[0]->ET * etmis[0]->ET);
      
      double MT = TMath::Sqrt(2.0 * etmis[0]->ET * mu[0]->PT * (1.0 -  cos(etmis[0]->Phi -  mu[0]->Phi)));
      
      sthist_mu->Fill(ST, weight);
      lphist_mu->Fill(LP, weight);
      mthist_mu->Fill(MT, weight);
	  
      if(ST > 250.0 && ST < 350.0) {
	if(LP<0.15) {
	  if(total_ht > 500.0 && total_ht <= 750.0) { mSigPred.at(0)+=weight; }
	  if(total_ht > 750.0 && total_ht <= 1000.0) { mSigPred.at(1)+=weight; }
	  if(total_ht > 1000.0) { mSigPred.at(2)+=weight; }
	}
      }
      if(ST > 350.0 && ST < 450.0) {
	if(LP<0.15) {
	  if(total_ht > 500.0 && total_ht <= 750.0) { mSigPred.at(3)+=weight; }
	  if(total_ht > 750.0 && total_ht <= 1000.0) { mSigPred.at(4)+=weight; }
	  if(total_ht > 1000.0) { mSigPred.at(5)+=weight; }
	}
      }
      if(ST > 450.0) {
	if(LP<0.15) {
	  if(total_ht > 500.0 && total_ht <= 750.0) { mSigPred.at(6)+=weight; }
	  if(total_ht > 750.0 && total_ht <= 1000.0) { mSigPred.at(7)+=weight; }
	  if(total_ht > 1000.0) { mSigPred.at(8)+=weight; }
	}
      }
    } 
  }
    
  if(goodjets.GetEntries() >=3 && muveto.GetEntries() == 0 && ele.GetEntries() == 1 && eleveto.GetEntries() == 1) {

    leadingjetpt->Fill(goodjets[0]->PT, weight);
    double total_ht = 0.0;
    unsigned int numbtags=0;
    for(unsigned int numjets = 0; numjets < goodjets.GetEntries(); numjets++) {
      total_ht += goodjets[numjets]->PT;
      if(goodjets[numjets]->Btag) { numbtags++; }
    }

    //btagrate->Fill(numbtags);
      
    hthist->Fill(total_ht, weight);
	
    if(total_ht > 500.0) {
	  
      double ST = (ele[0]->PT + etmis[0]->ET);
      
      double LP = ((ele[0]->Px * etmis[0]->Px) + (ele[0]->Py * etmis[0]->Py)) / (etmis[0]->ET * etmis[0]->ET);
      
      double MT = TMath::Sqrt(2.0 * etmis[0]->ET * ele[0]->PT * (1.0 -  cos(etmis[0]->Phi -  ele[0]->Phi)));
      
      sthist_el->Fill(ST, weight);
      lphist_el->Fill(LP, weight);
      mthist_el->Fill(MT, weight);
      
      if(ST > 250.0 && ST < 350.0) {
	if(LP<0.15) {
	  if(total_ht > 500.0 && total_ht <= 750.0) { mSigPred.at(9)+=weight; }
	  if(total_ht > 750.0 && total_ht <= 1000.0) { mSigPred.at(10)+=weight; }
	  if(total_ht > 1000.0) { mSigPred.at(11)+=weight; }
	}
      }
      if(ST > 350.0 && ST < 450.0) {
	if(LP<0.15) {
	  if(total_ht > 500.0 && total_ht <= 750.0) { mSigPred.at(12)+=weight; }
	  if(total_ht > 750.0 && total_ht <= 1000.0) { mSigPred.at(13)+=weight; }
	  if(total_ht > 1000.0) { mSigPred.at(14)+=weight; }
	}
      }
      if(ST > 450.0) {
	if(LP<0.15) {
	  if(total_ht > 500.0 && total_ht <= 750.0) { mSigPred.at(15)+=weight; }
	  if(total_ht > 750.0 && total_ht <= 1000.0) { mSigPred.at(16)+=weight; }
	  if(total_ht > 1000.0) { mSigPred.at(17)+=weight; }
	}
      }
    }      
  }

  return;
}
