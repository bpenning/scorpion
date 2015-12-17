#include "analysis_alphat13T.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
Alphat13T::Alphat13T(const std::string & name, 
                     const std::string & experiment,
                     const unsigned int & numBins) :
    AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
Alphat13T::Alphat13T(const std::string & name, 
                     const std::string & experiment, 
                     const unsigned int & numBins,
                     const double & intlumi, 
                     //const std::vector<int> & datayields,
                     const std::vector<double> & bgpred) :
    AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

Alphat13T::Alphat13T(const std::string & name, 
                     const std::string & experiment, 
                     const unsigned int & numBins,
                     const double & intlumi, 
                     const std::vector<double> & bgpred,
                     const std::vector<double> & bgpreduncert,
                     const std::vector<int> & datayields,
                     const std::string & fitmode,
                     const bool & calculateR) :
    AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


Alphat13T::~Alphat13T() {}



void Alphat13T::initHistos() {
	andir->cd();
	event_weight = new TH1D("event_weight",";weight_num;entries",1000,1940e9,1960e9);
	cut_sel = new TH1D("cut_selection","cut;Entries",7,-0.5,6.5);
	leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
	hthist = new TH1D("hthist", ";H_{T} [GeV];Entries",250,-5.,2495.);
	mhthist = new TH1D("mhthist", ";Missing H_{T} [GeV];Entries",200,-5.,1995.);
	calomethist = new TH1D("calomethist", ";CALO Missing E_{T} [GeV];Entries",200,-5.,1995.);
	athist = new TH1D("athist", ";#alpha_{T};Normalised",200,-0.005,1.995);
	athist2jets = new TH1D("athist2jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
	athist3jets = new TH1D("athist3jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
	athist4jets = new TH1D("athist4jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
	athist5jets = new TH1D("athist5jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
	biasedDPhiHist = new TH1D("bdhpi",";;",36,0.,3.6);
	njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
	bjets = new TH1D("bjets", ";N_{jets};Entries",10,-0.5,9.5);
	ejets = new TH1D("ejets", ";N_{jets};Entries",10,-0.5,9.5);
	mjets = new TH1D("mjets", ";N_{jets};Entries",10,-0.5,9.5);
	mht_over_ht = new TH1D("mht_over_ht",";MH_{T}/H_{T};Entries",200,-0.005, 1.995);
	btagrate = new TH1D("btagrate",";#b-tags;Entries",10,-0.5,9.5);
	calomet_vs_mht = new TH2D("calomet_vs_mht",";caloMET;MHT",200,-5.,1995., 200,-5.,1995.);
	ht_vs_mht_pre_alphaT = new TH2D("ht_vs_mht_pre_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
	ht_vs_mht_post_alphaT = new TH2D("ht_vs_mht_post_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
	TDirectory * alphaStatsInput = andir->mkdir("alphaStatsInput");
    TDirectory * plotSummaryDir = andir->mkdir("plotSummary");
	TString categoryList[32] = {
	    "eq0b_eq1j","eq1b_eq1j",
	"eq0b_eq2j","eq1b_eq2j","eq2b_eq2j",
	"eq0b_eq3j","eq1b_eq3j","eq2b_eq3j","ge3b_eq3j",
	"eq0b_eq4j","eq1b_eq4j","eq2b_eq4j","ge3b_eq4j",
        "eq0b_ge5j","eq1b_ge5j","eq2b_ge5j","ge3b_ge5j",
        "eq0b_eq2a","eq1b_eq2a","eq2b_eq2a",
	"eq0b_eq3a","eq1b_eq3a","eq2b_eq3a","ge3b_eq3a",
	"eq0b_eq4a","eq1b_eq4a","eq2b_eq4a","ge3b_eq4a",
	"eq0b_ge5a","eq1b_ge5a","eq2b_ge5a","ge3b_ge5a"};
    double htBins[9] = {200.,250.,300.,350.,400.,500.,600.,800.,10000.};

    plotSummaryDir->cd();
    plotSummary = new TH2D("h_ht_mht_all","", 8,htBins, 32,0,32);
    plotSummary->Sumw2();
    for (int bin=1; bin<=plotSummary->GetNbinsY(); ++bin)
        plotSummary->GetYaxis()->SetBinLabel(bin,categoryList[bin-1]);

	for (int cat = 0; cat < 32;cat++)
	{
	    TDirectory * catDir = alphaStatsInput->mkdir(categoryList[cat]);
	    catDir->cd();
	    alphaStatsInputHists[cat] = new TH2D("h_mht_"+categoryList[cat],"MHT",120,0.,3000.,60,0.,1500.);
        alphaStatsInputHists[cat]->Sumw2();
	}
	
	event_weight->SetBit(TH1::kCanRebin);
}

void Alphat13T::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {


  andir->cd();

  mCounter+=(weight); //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> ele=goodleptons(treereader->GetElec(), 10.0, 2.5,2.5,2.6);         //the central isolated electrons, pt > PT_ELEC GeV
  std::vector<jlepton> mu=goodleptons(treereader->GetMuon(), 10.0, 2.5,2.5,2.6);          //the central isolated muons, pt > PT_MUON GeV
  std::vector<jphoton> photon=goodphotons(treereader->GetPhoton(), 25.0, 2.5,2.5,2.6);
  std::vector<jtrack> track=trackSkim(treereader->GetIsoChargedTrack(), 10.0);
  std::vector<jjet> badjets=badjetsSkim(treereader->GetJet(),40.0,3.0); //check for jets which we should veto

//std::vector<jjet> goodjets100=goodjetsSkim(treereader->GetJet(),100,3.0); //check for jets which we should analyse
  std::vector<jjet> goodjets40=goodjetsSkim(treereader->GetJet(),40,3.0); //check for jets which we should analyse
  std::vector<jjet> goodjets20=goodjetsSkim(treereader->GetJet(),20,3.0); //check for jets which we should analyse

  std::vector<jjet> etmis=treereader->GetMet(); //Missing transverse energy array
  double calo_met = etmis[0].Et();
  calomethist->Fill(calo_met, weight);
  bjets->Fill(badjets.size());
  ejets->Fill(ele.size());
  mjets->Fill(mu.size());

  cut_sel->Fill(0.,weight);
  event_weight->Fill(weight);
  energy_sums esums = make_energy_sums(goodjets40);
  if(badjets.size() == 0 && ele.size() == 0 && mu.size() == 0) {
      if(track.size() == 0) {   
	  cut_sel->Fill(1.,weight);
	  if(esums.pass_quality_cuts) { //this cut contains njets>=1, HT check requirements and leading/sub-leading requirements
	      cut_sel->Fill(2.,weight);
	      njets->Fill(esums.njets, weight);
	      btagrate->Fill(esums.nbtags, weight);
	      hthist->Fill(esums.total_ht, weight);
              ht_vs_mht_pre_alphaT->Fill(esums.total_ht, esums.total_mht, weight);

	      if(esums.total_ht > 200.0) {
		  cut_sel->Fill(3.,weight);
		  mhthist->Fill(esums.total_mht, weight);
		  mht_over_ht->Fill(esums.total_mht/esums.total_ht, weight);
		  //calomethist->Fill(calo_met);
		  calomet_vs_mht->Fill(calo_met, esums.total_mht, weight);

          if(esums.total_mht > 130.) {

		  if(esums.total_mht/calo_met < 1.25) {
		      double alpha_t;
		      bool pseudoSize;
		      double biasedDPhi;
		      cut_sel->Fill(4.,weight);
		      std::vector<bool> pseudo;

		      if (goodjets40.size() == 1){
			  alpha_t = 10;
			  if (goodjets20.size() == 1){
			      biasedDPhi = 10;
			  }
			  else {
			      biasedDPhi = makeBiasedDPhi(goodjets20);
			  }
			  pseudoSize = pseudo.size() == esums.etvec.size();
		      }
		      else {
			  alpha_t = alphat(esums.etvec, esums.pxvec, esums.pyvec, pseudo, true);
			  biasedDPhi = makeBiasedDPhi(goodjets40);
			  pseudoSize = true;
		      }
		      biasedDPhiHist->Fill(biasedDPhi,weight);
		      ///////////////////////////////
		      //PSEUDO SITV HERE - MUST REMOVE
		      double SITVweight = weight;
		      ///////////////////////////

		      if (pseudoSize) {
			  cut_sel->Fill(5.,weight);

			  if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t < 0.65) return;
			  if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t < 0.60) return;
			  if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t < 0.55) return;
			  if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t < 0.55) return;
			  if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t < 0.53) return;
			  if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t < 0.52) return;
			  if(esums.total_ht > 600.0 && esums.total_ht <= 800.0 && alpha_t < 0.52) return;

			  if (biasedDPhi > 0.5){
			      if (esums.njets == 1){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[0]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,0.5,weight);
				  }
				  else if (esums.nbtags == 1) {
				      alphaStatsInputHists[1]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,1.5,weight);
				  }
			      }
			      else if (esums.njets == 2 && !esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[2]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,2.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[3]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,3.5,weight);
				  }  
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[4]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,4.5,weight);
				  }
			      }
			      else if (esums.njets == 3 && !esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[5]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,5.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[6]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,6.5,weight);
				  }
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[7]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,7.5,weight);
				  }
				  else if (esums.nbtags > 2) {
				      alphaStatsInputHists[8]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,8.5,weight);
				  }
			      }
			      else if (esums.njets == 4 && !esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[9]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,9.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[10]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,10.5,weight);
				  }
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[11]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,11.5,weight);
				  }
				  else if (esums.nbtags > 2) {
				      alphaStatsInputHists[12]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,12.5,weight);
				  }
			      }
			      else if (esums.njets > 4  && !esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[13]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,13.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[14]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,14.5,weight);
				  }
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[15]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,15.5,weight);
				  }
				  else if (esums.nbtags > 2) {
				      alphaStatsInputHists[16]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,16.5,weight);
				  }
			      }
			      else if (esums.njets == 2 && esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[17]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,17.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[18]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,18.5,weight);
				  }
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[19]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,19.5,weight);
				  }
			      }
			      else if (esums.njets == 3 && esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[20]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,20.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[21]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,21.5,weight);
				  }
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[22]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,22.5,weight);
				  }
				  else if (esums.nbtags > 2) {
				      alphaStatsInputHists[23]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,23.5,weight);
				  }
			      }
			      else if (esums.njets == 4 && esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[24]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,24.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[25]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,25.5,weight);
				  }
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[26]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,26.5,weight);
				  }
				  else if (esums.nbtags > 2) {
				      alphaStatsInputHists[27]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,27.5,weight);
				  }
			      }
			      else if (esums.njets > 4  && esums.asym){
				  if (esums.nbtags == 0) {
				      alphaStatsInputHists[28]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,28.5,weight);
				  }
				  if (esums.nbtags == 1) {
				      alphaStatsInputHists[29]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,29.5,weight);
				  }
				  else if (esums.nbtags == 2) {
				      alphaStatsInputHists[30]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,30.5,weight);
				  }
				  else if (esums.nbtags > 2) {
				      alphaStatsInputHists[31]->Fill(esums.total_ht,esums.total_mht,SITVweight);
				      plotSummary->Fill(esums.total_ht,31.5,weight);
				  }
			      }
			      athist->Fill(alpha_t, weight);
			      ht_vs_mht_post_alphaT->Fill(esums.total_ht, esums.total_mht, weight);

			      if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(0)+=(SITVweight); }
			      if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(1)+=(SITVweight); }
			      if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(2)+=(SITVweight); }
			      if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(3)+=(SITVweight); }
			      if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(4)+=(SITVweight); }
			      if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(5)+=(SITVweight); }
			      if(esums.total_ht > 600.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(6)+=(SITVweight); }
			      if(esums.total_ht > 800.0 && alpha_t > 0.0) { mSigPred.at(7)+=(SITVweight); }
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
