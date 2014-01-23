#include "analysis_alphatb8.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
AlphaTb8::AlphaTb8(const std::string & name, 
		   const std::string & experiment,
		   const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
AlphaTb8::AlphaTb8(const std::string & name, 
		   const std::string & experiment, 
		   const unsigned int & numBins,
		   const double & intlumi, 
		   //const std::vector<int> & datayields,
		   const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

AlphaTb8::AlphaTb8(const std::string & name, 
		   const std::string & experiment, 
		   const unsigned int & numBins,
		   const double & intlumi, 
		   const std::vector<double> & bgpred,
		   const std::vector<double> & bgpreduncert,
		   const std::vector<int> & datayields,
		   const std::string & fitmode,
		   const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


AlphaTb8::~AlphaTb8() {}

TSimpleArray<TRootJet> AlphaTb8::SubArrayGoodJets(const TClonesArray *JET, const float & pt, const float & eta) {

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

TSimpleArray<TRootElectron> AlphaTb8::SubArrayEl(const TClonesArray *ELEC, const float & pt, const float & eta) {
  TIter itElec((TCollection*)ELEC);
  TRootElectron *elec;
  itElec.Reset();
  TSimpleArray<TRootElectron> array;
  while( (elec = (TRootElectron*) itElec.Next()) )
    {
      //if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
      //double reliso = (elec->SumEt + elec->SumPt)/elec->PT;
      //reliso > 0.065
      if(elec->PT<pt || !elec->IsolFlag || fabs(elec->Eta) > eta) continue;
      array.Add(elec);
    }
  return array;
}

TSimpleArray<TRootMuon> AlphaTb8::SubArrayMu(const TClonesArray *MUON, const float & pt, const float & eta) {
  TIter itMuon((TCollection*)MUON);
  TRootMuon *muon;
  itMuon.Reset();
  TSimpleArray<TRootMuon> array;
  while( (muon = (TRootMuon*) itMuon.Next()) )
    {
      //double reliso = (muon->SumEt + muon->SumPt)/muon->PT;
      //if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
      //reliso > 0.10
      if(muon->PT<pt || !muon->IsolFlag || fabs(muon->Eta) > eta) continue;
      array.Add(muon);
    }
  return array;
}

TSimpleArray<TRootJet> AlphaTb8::SubArrayBadJets(const TClonesArray *JET, const float & pt, const float & eta) {

  TIter itJet((TCollection*)JET);
  TRootJet *jet;
  itJet.Reset();
  TSimpleArray<TRootJet> array;
  while( (jet = (TRootJet*) itJet.Next()) ) {
    //check if any jet has Pt>50 and |eta|>3.0
    if(jet->PT > pt && fabs(jet->Eta) > eta) {
      array.Add(jet);
    }
  }
  return array;
}

TSimpleArray<TRootETmis> AlphaTb8::makeETM(const TClonesArray *ETMISS) {
  TIter itEtMiss((TCollection*)ETMISS);
  TRootETmis *etm;
  itEtMiss.Reset();
  TSimpleArray<TRootETmis> array;
  while( (etm = (TRootETmis*) itEtMiss.Next()) ) {
    array.Add(etm);
  }
  return array;
}


void AlphaTb8::initHistos() {
  andir->cd();
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  hthist = new TH1D("hthist", ";H_{T} [GeV];Entries",250,-5.,2495.);
  mhthist = new TH1D("mhthist", ";Missing H_{T} [GeV];Entries",200,-5.,1995.);
  calomethist = new TH1D("calomethist", ";CALO Missing E_{T} [GeV];Entries",200,-5.,1995.);
  athist = new TH1D("athist", ";#alpha_{T};Normalised",200,-0.005,1.995);
  njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
  mht_over_ht = new TH1D("mht_over_ht",";MH_{T}/H_{T};Entries",200,-0.005, 1.995);
  btagrate = new TH1D("btagrate",";#b-tags;Entries",10,-0.5,9.5);
  calomet_vs_mht = new TH2D("calomet_vs_mht",";caloMET;MHT",200,-5.,1995., 200,-5.,1995.);
  ht_vs_mht_pre_alphaT = new TH2D("ht_vs_mht_pre_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
  ht_vs_mht_post_alphaT = new TH2D("ht_vs_mht_post_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
}

void AlphaTb8::Run(const TreeReader & treereader, const TreeReader & gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  TSimpleArray<TRootElectron> ele=SubArrayEl(treereader.Elec(), 10.0, 2.4);         //the central isolated electrons, pt > PT_ELEC GeV
  TSimpleArray<TRootMuon>     mu=SubArrayMu(treereader.Muon(), 10.0, 2.1);          //the central isolated muons, pt > PT_MUON GeV
  TSimpleArray<TRootJet>      badjets=SubArrayBadJets(treereader.Jet(),50.0,3.0);   //check for jets which we should veto
  TSimpleArray<TRootJet>      goodjets275=SubArrayGoodJets(treereader.Jet(),37.0,3.0); //check for jets which we should analyse
  TSimpleArray<TRootJet>      goodjets325=SubArrayGoodJets(treereader.Jet(),43.0,3.0); //check for jets which we should analyse
  TSimpleArray<TRootJet>      goodjets375=SubArrayGoodJets(treereader.Jet(),50.0,3.0); //check for jets which we should analyse
  TSimpleArray<TRootETmis>    etmis=makeETM(treereader.ETMis()); //Missing transverse energy array

  //TSimpleArray<TRootJet>      goodjets2=SubArrayGoodJets(treereader.Jet(),10.0,3.0); //check for jets which we should analyse
  //TSimpleArray<TRootJet> goodjets2 = SubArrayGoodJets(treereader.Jet(),10.0,3.0); //check for jets which we should analyse

  double calo_met = etmis[0]->ET;
  calomethist->Fill(calo_met, weight);

  energy_sums esums = make_energy_sums(goodjets275, goodjets325, goodjets375);

  if(badjets.GetEntries() == 0 && ele.GetEntries() == 0 && mu.GetEntries() == 0) {

    //some histograms:
    //if(fabs(goodjets275[0]->Eta) < 2.4) { leadingjetpt->Fill(goodjets275[0]->PT); }

    //if(goodjets[0]->PT > 100.0 && fabs(goodjets[0]->Eta) < 2.4) {
    //  if(goodjets[1]->PT> 100.0) {
    if(esums.pass_quality_cuts) { //this cut contains njets>=2, HT check requirements and leading/sub-leading requirements
      njets->Fill(esums.njets, weight);
      btagrate->Fill(esums.nbtags, weight);
      //for alpha_t:
      //std::vector<double> etvec;
      //std::vector<double> pxvec;
      //std::vector<double> pyvec;

      //double total_ht = 0.0;

      //unsigned int numbtags=0;

      //for(unsigned int numjets = 0; numjets < goodjets.GetEntries(); numjets++) {
        //total_ht += goodjets[numjets]->PT;
	//etvec.push_back(goodjets[numjets]->PT);
        //pxvec.push_back(goodjets[numjets]->Px);
	//pyvec.push_back(goodjets[numjets]->Py);
	//if(goodjets[numjets]->Btag) { numbtags++; }
      //}

      //btagrate->Fill(numbtags);

      hthist->Fill(esums.total_ht, weight);
      ht_vs_mht_pre_alphaT->Fill(esums.total_ht, esums.total_mht, weight);

	if(esums.total_ht > 275.0) {
	  //at this point, need to fill mht / met histogram.
	  //double mht_x1 = 0.0;
	  //double mht_y1 = 0.0;

	  //for(unsigned int numjets = 0; numjets < goodjets.GetEntries(); numjets++) {
	  //  mht_x1 += goodjets[numjets]->Px;
	  //  mht_y1 += goodjets[numjets]->Py;
	  //}

	  //double mht_x2 = 0.0;
	  //double mht_y2 = 0.0;

	  //for(unsigned int numjets = 0; numjets < goodjets2.GetEntries(); numjets++) {
	  //  mht_x2 += goodjets2[numjets]->Px;
	  //  mht_y2 += goodjets2[numjets]->Py;
	  //}

	  //double mht = TMath::Sqrt((mht_x1 * mht_x1) +(mht_y1 * mht_y1));
	  //double met = TMath::Sqrt((mht_x2 * mht_x2) +(mht_y2 * mht_y2)); //define met as vector sum of jets with pt>10, |eta|<3

	  mhthist->Fill(esums.total_mht, weight);
	  mht_over_ht->Fill(esums.total_mht/esums.total_ht, weight);
	  //calomethist->Fill(calo_met);
	  calomet_vs_mht->Fill(calo_met, esums.total_mht, weight);

	  if(esums.total_mht/calo_met < 1.25) {
	    std::vector<bool> pseudo;
	    double alpha_t = alphat(esums.etvec, esums.pxvec, esums.pyvec, pseudo, true);
	    if ( pseudo.size() == esums.etvec.size() ) {
	      athist->Fill(alpha_t, weight);
	      if(alpha_t > 0.55) {
		//btagrate->Fill(esums.nbtags, weight);
		ht_vs_mht_post_alphaT->Fill(esums.total_ht, esums.total_mht, weight);
		if(esums.nbtags == 0) {
		  if(esums.njets < 4) { //the >=2 cut is in the make energy sums method above
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(0)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(1)+=weight; }
		    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(2)+=weight; }
		    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(3)+=weight; }
		    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(4)+=weight; }
		    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(5)+=weight; }
		    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(6)+=weight; }
		    if(esums.total_ht > 875.0) { mSigPred.at(7)+=weight; }
		  } else { //must be >=4 jets
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(8)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(9)+=weight; }
		    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(10)+=weight; }
		    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(11)+=weight; }
		    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(12)+=weight; }
		    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(13)+=weight; }
		    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(14)+=weight; }
		    if(esums.total_ht > 875.0) { mSigPred.at(15)+=weight; }
		  }
		} else if(esums.nbtags == 1) {
		  if(esums.njets < 4) {
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(16)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(17)+=weight; }
		    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(18)+=weight; }
		    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(19)+=weight; }
		    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(20)+=weight; }
		    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(21)+=weight; }
		    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(22)+=weight; }
		    if(esums.total_ht > 875.0) { mSigPred.at(23)+=weight; }		   
		  } else {
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(24)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(25)+=weight; }
		    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(26)+=weight; }
		    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(27)+=weight; }
		    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(28)+=weight; }
		    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(29)+=weight; }
		    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(30)+=weight; }
		    if(esums.total_ht > 875.0) { mSigPred.at(31)+=weight; }
		  }
		} else if(esums.nbtags == 2) {
		  if(esums.njets < 4) {
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(32)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(33)+=weight; }
		    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(34)+=weight; }
		    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(35)+=weight; }
		    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(36)+=weight; }
		    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(37)+=weight; }
		    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(38)+=weight; }
		    if(esums.total_ht > 875.0) { mSigPred.at(39)+=weight; }		    
		  } else {
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(40)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(41)+=weight; }
		    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(42)+=weight; }
		    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(43)+=weight; }
		    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(44)+=weight; }
		    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(45)+=weight; }
		    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(46)+=weight; }
		    if(esums.total_ht > 875.0) { mSigPred.at(47)+=weight; }
		  }
		} else if(esums.nbtags == 3) {
		  if(esums.njets >=4) {
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(48)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(49)+=weight; }
		    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(50)+=weight; }
		    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(51)+=weight; }
		    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(52)+=weight; }
		    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(53)+=weight; }
		    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(54)+=weight; }
		    if(esums.total_ht > 875.0) { mSigPred.at(55)+=weight; }
		  }
		} else {
		  if(esums.njets >=4) {
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(56)+=weight; }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(57)+=weight; }
		    if(esums.total_ht > 375.0) { mSigPred.at(58)+=weight; }
		  }
		}
	      }
	    }
	  }
	//}
      }
    }
  }
  return;
}
