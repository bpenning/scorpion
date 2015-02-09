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
	jet1pt = new TH1D("firstjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
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
	event_weight->SetBit(TH1::kCanRebin);
}

void Alphat13T::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader.GetEntries() << std::endl;

  andir->cd();

  mCounter+=(weight); //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> ele=goodleptons(treereader->GetElec(), 10.0, 2.5,2.5,2.6);         //the central isolated electrons, pt > PT_ELEC GeV
  std::vector<jlepton> mu=goodleptons(treereader->GetMuon(), 10.0, 2.5,2.5,2.6);          //the central isolated muons, pt > PT_MUON GeV
  std::vector<jphoton> photon=goodphotons(treereader->GetPhoton(), 25.0, 2.5,2.5,2.6);
  std::vector<jtrack> track=trackSkim(treereader->GetIsoChargedTrack(), 10.0);
  std::vector<jjet> badjets=badjetsSkim(treereader->GetJet(),40.0,3.0); //check for jets which we should veto

  std::vector<jjet> goodjets100=goodjetsSkim(treereader->GetJet(),100,3.0); //check for jets which we should analyse
  std::vector<jjet> goodjets40=goodjetsSkim(treereader->GetJet(),40,3.0); //check for jets which we should analyse

//   std::vector<jjet> goodjets200=goodjetsSkim(treereader->GetJet(),36.7,3.0); //check for jets which we should analyse
//   std::vector<jjet> goodjets275=goodjetsSkim(treereader->GetJet(),36.7,3.0); //check for jets which we should analyse
//   std::vector<jjet> goodjets325=goodjetsSkim(treereader->GetJet(),43.3,3.0); //check for jets which we should analyse
//   std::vector<jjet> goodjets375=goodjetsSkim(treereader->GetJet(),50.0,3.0); //check for jets which we should analyse
  std::vector<jjet> etmis=treereader->GetMet(); //Missing transverse energy array
  if (goodjets100.size() != 0){
    jet1pt->Fill(goodjets100[0].Pt());
  } 
  if (goodjets100.size() > 1)
    {
      jet2pt->Fill(goodjets100[1].Pt());
    } 
  double calo_met = etmis[0].Et();
  calomethist->Fill(calo_met, weight);
  bjets->Fill(badjets.size());
  ejets->Fill(ele.size());
  mjets->Fill(mu.size());

  cut_sel->Fill(0.,weight);
  event_weight->Fill(weight);
  energy_sums esums = make_energy_sums(goodjets40, goodjets100);
  if(badjets.size() == 0 && ele.size() == 0 && mu.size() == 0) {
    if(track.size() == 0) {   
      cut_sel->Fill(1.,weight);

      //some histograms:
      //if(fabs(goodjets275[0]->Eta) < 2.4) { leadingjetpt->Fill(goodjets275[0]->PT); }

      //if(goodjets[0]->PT > 100.0 && fabs(goodjets[0]->Eta) < 2.4) {
      //  if(goodjets[1]->PT> 100.0) {
      if(esums.pass_quality_cuts) { //this cut contains njets>=2, HT check requirements and leading/sub-leading requirements
	      cut_sel->Fill(2.,weight);
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

	      if(esums.total_ht > 200.0) {
          cut_sel->Fill(3.,weight);
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
            cut_sel->Fill(4.,weight);
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
              cut_sel->Fill(5.,weight);
              athist->Fill(alpha_t, weight);
              //cut_sel->Fill(6.,weight);
              //btagrate->Fill(esums.nbtags, weight);
              ht_vs_mht_post_alphaT->Fill(esums.total_ht, esums.total_mht, weight);

                ///////////////////////////////
                //PSEUDO SITV HERE - MUST REMOVE
                double SITVweight = weight * 0.9;
                ///////////////////////////
                if(esums.nbtags == 0) {
                  if(esums.njets == 2 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(1)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(2)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(3)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(4)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(5)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(6)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(7)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(8)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(9)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(10)+=(SITVweight); }
                  }
                  if(esums.njets == 3 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(11)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6)  { mSigPred.at(12)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(13)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(14)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(15)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(16)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(17)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(18)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(19)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(20)+=(SITVweight); }
                  }
                  if(esums.njets == 4 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(21)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(22)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(23)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(24)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(25)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(26)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(27)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(28)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(29)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(30)+=(SITVweight); }
                  }
                  if(esums.njets > 4 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(31)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(32)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(33)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(34)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(35)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(36)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(37)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(38)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(39)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(40)+=(SITVweight); }
                  }
                  
                }
                 else if(esums.nbtags == 1) {
                  if(esums.njets == 2 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(41)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(42)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(43)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(44)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(45)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(46)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(47)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(48)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(49)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(50)+=(SITVweight); }
                  }
                  if(esums.njets == 3 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(51)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6)  { mSigPred.at(52)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(53)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(54)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(55)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(56)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(57)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(58)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(59)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(60)+=(SITVweight); }
                  }
                  if(esums.njets == 4 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(61)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(62)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(63)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(64)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(65)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(66)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(67)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(68)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(69)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(70)+=(SITVweight); }
                  }
                  if(esums.njets > 4 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(71)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(72)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(73)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(74)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(75)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(76)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(77)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(78)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(79)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(80)+=(SITVweight); }
                  }
                 }

                 else if(esums.nbtags == 2) {
                   if(esums.njets == 2 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                     if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(81)+=(SITVweight); }
                     if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(82)+=(SITVweight); }
                     if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(83)+=(SITVweight); }
                     if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(84)+=(SITVweight); }
                     if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(85)+=(SITVweight); }
                     if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(86)+=(SITVweight); }
                     if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(87)+=(SITVweight); }
                     if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(88)+=(SITVweight); }
                     if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(89)+=(SITVweight); }
                     if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(90)+=(SITVweight); }
                  }
                  if(esums.njets == 3 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(91)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6)  { mSigPred.at(92)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(93)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(94)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(95)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(96)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(97)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(98)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(99)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(100)+=(SITVweight); }
                  }
                  if(esums.njets == 4 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(101)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(102)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(103)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(104)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(105)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(106)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(107)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(108)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(109)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(110)+=(SITVweight); }
                  }
                  if(esums.njets > 4 && goodjets100[0].Pt() > 100 && goodjets100[1].Pt() > 100){
                    if(esums.total_ht > 200.0 && esums.total_ht <= 250.0 && alpha_t > 0.65) { mSigPred.at(111)+=(SITVweight); }
                    if(esums.total_ht > 250.0 && esums.total_ht <= 300.0 && alpha_t > 0.6) { mSigPred.at(112)+=(SITVweight); }
                    if(esums.total_ht > 300.0 && esums.total_ht <= 350.0 && alpha_t > 0.55) { mSigPred.at(113)+=(SITVweight); }
                    if(esums.total_ht > 350.0 && esums.total_ht <= 400.0 && alpha_t > 0.55) { mSigPred.at(114)+=(SITVweight); }
                    if(esums.total_ht > 400.0 && esums.total_ht <= 500.0 && alpha_t > 0.53) { mSigPred.at(115)+=(SITVweight); }
                    if(esums.total_ht > 500.0 && esums.total_ht <= 600.0 && alpha_t > 0.52) { mSigPred.at(116)+=(SITVweight); }
                    if(esums.total_ht > 600.0 && esums.total_ht <= 700.0 && alpha_t > 0.52) { mSigPred.at(117)+=(SITVweight); }
                    if(esums.total_ht > 700.0 && esums.total_ht <= 800.0 && alpha_t > 0.52) { mSigPred.at(118)+=(SITVweight); }
                    if(esums.total_ht > 800.0 && esums.total_ht <= 900.0 && alpha_t > 0.505) { mSigPred.at(119)+=(SITVweight); }
                    if(esums.total_ht > 9000.0 && alpha_t > 0.65) { mSigPred.at(120)+=(SITVweight); }
                  }
                 }
                
                
                
                //      else if(esums.nbtags >= 4) {
//                  }
            
            }//if ( pseudo.size() == esums.etvec.size() ) {
          }
        }
      }
    }
  }
  return;
  //}
}
