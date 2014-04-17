#include "analysis_alphatb.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
AlphaTb::AlphaTb(const std::string & name, 
		 const std::string & experiment,
		 const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
AlphaTb::AlphaTb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 //const std::vector<int> & datayields,
		 const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

AlphaTb::AlphaTb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 const std::vector<double> & bgpred,
		 const std::vector<double> & bgpreduncert,
		 const std::vector<int> & datayields,
		 const std::string & fitmode,
		 const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


AlphaTb::~AlphaTb() {}



void AlphaTb::initHistos() {
  andir->cd();
  cut_sel = new TH1I("cut_selection","cut;Entries",6,-0.5,5.5);
  leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  hthist = new TH1D("hthist", ";H_{T} [GeV];Entries",250,-5.,2495.);
  mhthist = new TH1D("mhthist", ";Missing H_{T} [GeV];Entries",200,-5.,1995.);
  calomethist = new TH1D("calomethist", ";CALO Missing E_{T} [GeV];Entries",200,-5.,1995.);
  athist = new TH1D("athist", ";#alpha_{T};Normalised",200,-0.005,1.995);
  athist2jets = new TH1D("athist2jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
  athist3jets = new TH1D("athist3jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
  athist4jets = new TH1D("athist4jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
  athist5jets = new TH1D("athist5jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
  jet1pt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  jet2pt = new TH1D("secondjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
  deltaphi = new TH1D("deltaphi",";;",72,-0.05,3.55);
  njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
  bjets = new TH1D("bjets", ";N_{jets};Entries",10,-0.5,9.5);
  ejets = new TH1D("ejets", ";N_{jets};Entries",10,-0.5,9.5);
  mjets = new TH1D("mjets", ";N_{jets};Entries",10,-0.5,9.5);
  mht_over_ht = new TH1D("mht_over_ht",";MH_{T}/H_{T};Entries",200,-0.005, 1.995);
  btagrate = new TH1D("btagrate",";#b-tags;Entries",10,-0.5,9.5);
  calomet_vs_mht = new TH2D("calomet_vs_mht",";caloMET;MHT",200,-5.,1995., 200,-5.,1995.);
  ht_vs_mht_pre_alphaT = new TH2D("ht_vs_mht_pre_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
  ht_vs_mht_post_alphaT = new TH2D("ht_vs_mht_post_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
}

void AlphaTb::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> ele=goodleptons(treereader->Elec(), 10.0, 2.4);         //the central isolated electrons, pt > PT_ELEC GeV
  std::vector<jlepton>     mu=goodleptons(treereader->Muon(), 10.0, 2.1);          //the central isolated muons, pt > PT_MUON GeV
  std::vector<jjet>      badjets=badjetsSkim(treereader->Jet(),50.0,3.0);   //check for jets which we should veto
  std::vector<jjet>      goodjets275=goodjetsSkim(treereader->Jet(),37.0,3.0); //check for jets which we should analyse
  std::vector<jjet>      goodjets325=goodjetsSkim(treereader->Jet(),43.0,3.0); //check for jets which we should analyse
  std::vector<jjet>     goodjets375=goodjetsSkim(treereader->Jet(),50.0,3.0); //check for jets which we should analyse
  std::vector<jjet>    etmis=treereader->ETMis(); //Missing transverse energy array
  //TSimpleArray<TRootJet>      goodjets2=SubArrayGoodJets(treereader.Jet(),10.0,3.0); //check for jets which we should analyse
  //TSimpleArray<TRootJet> goodjets2 = SubArrayGoodJets(treereader.Jet(),10.0,3.0); //check for jets which we should analyse
  if (goodjets275.size() != 0)
  {
      jet1pt->Fill(goodjets275[0].Pt());
  } 
  if (goodjets275.size() > 1)
  {
      jet2pt->Fill(goodjets275[1].Pt());
  } 
  double calo_met = etmis[0].Et();
  calomethist->Fill(calo_met, weight);
  bjets->Fill(badjets.size());
  ejets->Fill(ele.size());
  mjets->Fill(mu.size());

  energy_sums esums = make_energy_sums(goodjets275, goodjets325, goodjets375);
  if(badjets.size() == 0 && ele.size() == 0 && mu.size() == 0) {
      cut_sel-> AddBinContent(1);

      //some histograms:
      //if(fabs(goodjets275[0]->Eta) < 2.4) { leadingjetpt->Fill(goodjets275[0]->PT); }

      //if(goodjets[0]->PT > 100.0 && fabs(goodjets[0]->Eta) < 2.4) {
      //  if(goodjets[1]->PT> 100.0) {
      if(esums.pass_quality_cuts) { //this cut contains njets>=2, HT check requirements and leading/sub-leading requirements
	  cut_sel-> AddBinContent(2);
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
	      cut_sel-> AddBinContent(3);
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
		  cut_sel-> AddBinContent(4);
		  std::vector<bool> pseudo;
		  double alpha_t = alphat(esums.etvec, esums.pxvec, esums.pyvec, pseudo, true);
		  if ( pseudo.size() == esums.etvec.size() ) {
		      double phi1=std::atan2(esums.pyvec[0],esums.pxvec[0]);
		      double phi2=std::atan2(esums.pyvec[1],esums.pxvec[1]);
		      switch (esums.njets){
			  case 2:
			      athist2jets->Fill(alpha_t , weight);
			      if (fabs(phi1-phi2) <3.14159265359) 
				  deltaphi->Fill(fabs(phi1-phi2));
			      else 
				  deltaphi->Fill(2.*3.14159265359-fabs(phi2-phi1));
			      break;
			  case 3:
			      athist3jets->Fill(alpha_t , weight);
			      break;
			  case 4:
			      athist4jets->Fill(alpha_t , weight);
			      break;
			  case 5:
			      athist5jets->Fill(alpha_t , weight);
			      break;
		      }
		      cut_sel-> AddBinContent(5);
		      athist->Fill(alpha_t, weight);
		      if(alpha_t > 0.55) {
			  cut_sel-> AddBinContent(6);
			  //btagrate->Fill(esums.nbtags, weight);
			  ht_vs_mht_post_alphaT->Fill(esums.total_ht, esums.total_mht, weight);
			  if(esums.nbtags == 0) {
			      if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(0)+=weight; }
			      if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(1)+=weight; }
			      if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(2)+=weight; }
			      if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(3)+=weight; }
			      if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(4)+=weight; }
			      if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(5)+=weight; }
			      if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(6)+=weight; }
			      if(esums.total_ht > 875.0) { mSigPred.at(7)+=weight; }
			  } else if(esums.nbtags == 1) {
			      if(esums.total_ht > 275.0 && esums.total_ht <= 325.0) { mSigPred.at(8)+=weight; }
			      if(esums.total_ht > 325.0 && esums.total_ht <= 375.0) { mSigPred.at(9)+=weight; }
			      if(esums.total_ht > 375.0 && esums.total_ht <= 475.0) { mSigPred.at(10)+=weight; }
			      if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(11)+=weight; }
			      if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(12)+=weight; }
			      if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(13)+=weight; }
			      if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(14)+=weight; }
			      if(esums.total_ht > 875.0) { mSigPred.at(15)+=weight; }
			  } else if(esums.nbtags == 2) {
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
		      }
		  }
	      }
	      //}
      }
  }
  }
  return;
  }
