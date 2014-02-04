#include ""

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
CmsSingleLepton20Fb::CmsSingleLepton20Fb(const std::string & name, 
		 const std::string & experiment,
		 const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
CmsSingleLepton20Fb::CmsSingleLepton20Fb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 //const std::vector<int> & datayields,
		 const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

CmsSingleLepton20Fb::CmsSingleLepton20Fb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 const std::vector<double> & bgpred,
		 const std::vector<double> & bgpreduncert,
		 const std::vector<int> & datayields,
		 const std::string & fitmode,
		 const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


CmsSingleLepton20Fb::~CmsSingleLepton20Fb() {}

TSimpleArray<TRootJet> CmsSingleLepton20Fb::SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta) {

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

TSimpleArray<TRootElectron> CmsSingleLepton20Fb::SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta) {
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

TSimpleArray<TRootMuon> CmsSingleLepton20Fb::SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta) {
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

TSimpleArray<TRootETmis> CmsSingleLepton20Fb::makeETM(const TClonesArray *ETMISS) {
  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();
  TSimpleArray<TRootETmis> array;
  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
    array.Add(etm);
  }
  return array;
}

bool CmsSingleLepton20Fb::checkforsecondjetphi(const TSimpleArray<TRootJet> & goodjets, const double & dphisep) {

  if(goodjets.GetEntries() == 2) {
    
    //check if deltaphi j1, j2 < 2.5
    double phijet1 = TMath::ATan2(goodjets[0]->Py, goodjets[0]->Px); //calculate jet phi
    double phijet2 = TMath::ATan2(goodjets[1]->Py, goodjets[1]->Px); //calculate jet phi
    double dphij12 = fabs(phijet1 - phijet2); //calculate the distance in phi between the jet total momentum
    
    //for two vectors in a 2D plane, they cannot be more than pi away from each other
    if(dphij12 > TMath::Pi()) {
      dphij12 = fabs(dphij12 - (2.0 * TMath::Pi()));
    }
    
    if (dphij12 < dphisep) {
      return true;
    } else {
      return false;
    }
  }
  
  return false;
   
}

void CmsSingleLepton20Fb::initHistos() {
  andir->cd();
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  calomet = new TH1D("calomet",";E_{T}^{miss} [GeV];Entries",200,-5.,1995.);
}

void CmsSingleLepton20Fb::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  TSimpleArray<TRootElectron> ele=SubArrayEl(treereader.Elec(), 10.0, 2.4); //the central isolated electrons, pt > PT_ELEC GeV
  TSimpleArray<TRootMuon> mu=SubArrayMu(treereader.Muon(), 10.0, 2.1); //the central isolated muons, pt > PT_MUON GeV
  TSimpleArray<TRootJet> goodjets=SubArrayGoodJets(treereader.Jet(), 30.0, 3.0); //check for jets which we should analyse
  TSimpleArray<TRootETmis> etmis=makeETM(treereader.ETMis()); //Missing transverse energy array

  double calo_met = etmis[0]->ET;

  if(goodjets.GetEntries() <= 2 && goodjets.GetEntries() > 0 && mu.GetEntries() == 0 && ele.GetEntries() == 0) {
      
    if(goodjets[0]->PT > 110.0 && fabs(goodjets[0]->Eta) < 2.4) {

      if(checkforsecondjetphi(goodjets, 2.5) || goodjets.GetEntries()==1) { 
	
	leadingjetpt->Fill(goodjets[0]->PT, weight);
	calomet->Fill(calo_met, weight);
	if(calo_met > 250.0) { mSigPred.at(0)+=weight; }
	if(calo_met > 300.0) { mSigPred.at(1)+=weight; }
	if(calo_met > 350.0) { mSigPred.at(2)+=weight; }
	if(calo_met > 400.0) { mSigPred.at(3)+=weight; }

      }

    }

  }
    

  return;
}
