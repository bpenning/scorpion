#include "analysis_alphat20b.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
AlphaT20b::AlphaT20b(const std::string & name, 
	const std::string & experiment,
	const unsigned int & numBins) :
    AnalysisBase(name, experiment, numBins) {}

    //constructor for object where you also want to run a limit and so other information
    //is necessary
    AlphaT20b::AlphaT20b(const std::string & name, 
	    const std::string & experiment, 
	    const unsigned int & numBins,
	    const double & intlumi, 
	    //const std::vector<int> & datayields,
	    const std::vector<double> & bgpred) :
	AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

	AlphaT20b::AlphaT20b(const std::string & name, 
		const std::string & experiment, 
		const unsigned int & numBins,
		const double & intlumi, 
		const std::vector<double> & bgpred,
		const std::vector<double> & bgpreduncert,
		const std::vector<int> & datayields,
		const std::string & fitmode,
		const bool & calculateR) :
	    AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


	    AlphaT20b::~AlphaT20b() {}



	    void AlphaT20b::initHistos() {
		andir->cd();
		event_weight = new TH1D("event_weight",";weight_num;entries",1000,1940e9,1960e9);
		cut_sel = new TH1D("cut_selection","cut;Entries",9,-0.5,8.5);
		cut_sel23b0 = new TH1D("cut_selection23b0","cut;Entries",9,-0.5,8.5);
		cut_selGe4b0 = new TH1D("cut_selectionGe4b0","cut;Entries",9,-0.5,8.5);
		cut_sel23b1 = new TH1D("cut_selection23b1","cut;Entries",9,-0.5,8.5);
		cut_selGe4b1 = new TH1D("cut_selectionGe4b1","cut;Entries",9,-0.5,8.5);
		leadingjetpt = new TH1D("leadingjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
		hthist = new TH1D("hthist", ";H_{T} [GeV];Entries",250,-5.,2495.);
		mhthist = new TH1D("mhthist", ";Missing H_{T} [GeV];Entries",200,-5.,1995.);
		calomethist = new TH1D("calomethist", ";CALO Missing E_{T} [GeV];Entries",200,-5.,1995.);
		athist = new TH1D("athist", ";#alpha_{T};Normalised",200,-0.005,1.995);
		athist2jets = new TH1D("athist2jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
		athist3jets = new TH1D("athist3jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
		athist4jets = new TH1D("athist4jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
		athist5jets = new TH1D("athist5jets", ";#alpha_{T};Normalised",200,-0.005,1.995);
		boostHisto = new TH1D("boostHisto", ";P_{T} [GeV];Entries",200,-5.,1995.);
		cTrack = new TH1D("cTrack", ";P_{T} [GeV];Entries",2000,-5.,195.);
		jet1pt = new TH1D("firstjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
		jet2pt = new TH1D("secondjetpt", ";P_{T} [GeV];Entries",200,-5.,1995.);
		deltaphi = new TH1D("deltaphi",";;",72,-0.05,3.55);
		njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
		bjets = new TH1D("bjets", ";N_{jets};Entries",10,-0.5,9.5);
		ejets = new TH1D("ejets", ";N_{jets};Entries",10,-0.5,9.5);
		mjets = new TH1D("mjets", ";N_{jets};Entries",10,-0.5,9.5);

		bDPhiHistoGe4At55b0 = new TH1D("bDPhiGe4At55b0", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At55b0 = new TH1D("bDPhi23At55b0", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHistoGe4At53b0 = new TH1D("bDPhiGe4At53b0", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At53b0 = new TH1D("bDPhi23At53b0", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHistoGe4At55b1= new TH1D("bDPhiGe4At55b1", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At55b1 = new TH1D("bDPhi23At55b1", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHistoGe4At53b1 = new TH1D("bDPhiGe4At53b1", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At53b1 = new TH1D("bDPhi23At53b1", ";bDPhi;Entries",11,0.0,3.3);

		bDPhiHistoGe4At55b0EcalCut = new TH1D("bDPhiGe4At55b0EcalCut", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At55b0EcalCut = new TH1D("bDPhi23At55b0EcalCut", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHistoGe4At53b0EcalCut = new TH1D("bDPhiGe4At53b0EcalCut", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At53b0EcalCut = new TH1D("bDPhi23At53b0EcalCut", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHistoGe4At55b1EcalCut= new TH1D("bDPhiGe4At55b1EcalCut", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At55b1EcalCut = new TH1D("bDPhi23At55b1EcalCut", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHistoGe4At53b1EcalCut = new TH1D("bDPhiGe4At53b1EcalCut", ";bDPhi;Entries",11,0.0,3.3);
		bDPhiHisto23At53b1EcalCut = new TH1D("bDPhi23At53b1EcalCut", ";bDPhi;Entries",11,0.0,3.3);

		bDPhiHisto23At53 = new TH1D("bDPhi23At53",";bDPhi;",11,0.0,3.3);
		bDPhiHistoGe4At53 = new TH1D("bDPhiGe4At53",";bDPhi;",11,0.0,3.3);
		bDPhiHisto23At53EcalCut = new TH1D("bDPhi23At53EcalCut",";bDPhi;",11,0.0,3.3);
		bDPhiHistoGe4At53EcalCut = new TH1D("bDPhiGe4At53EcalCut",";bDPhi;",11,0.0,3.3);

		mht_over_ht = new TH1D("mht_over_ht",";MH_{T}/H_{T};Entries",200,-0.005, 1.995);
		btagrate = new TH1D("btagrate",";#b-tags;Entries",10,-0.5,9.5);
		calomet_vs_mht = new TH2D("calomet_vs_mht",";caloMET;MHT",200,-5.,1995., 200,-5.,1995.);
		ht_vs_mht_pre_alphaT = new TH2D("ht_vs_mht_pre_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
		ht_vs_mht_post_alphaT = new TH2D("ht_vs_mht_post_alphaT",";H_{T} [GeV]; Missing H_{T} [GeV]", 250, -5., 2495., 200, -5, 1995.);
		ecal_map = new TH2D("ecalMap",";phi;eta;",100,-3.141,3.141,100,-3.0,3.0);
		event_weight->SetBit(TH1::kCanRebin);
		std::vector<jecal> ecalMap = getMapEcal();
		for (int ecal = 0; ecal < ecalMap.size(); ecal++)
		{
		    ecal_map->Fill(ecalMap[ecal].Phi(),ecalMap[ecal].Eta());
		}
	    }
double AlphaT20b::biasedDPhiMin(std::vector<jjet> inJets)
{
    TLorentzVector mht;
    //jjet biasedMht;
    double biasedDPhiMin= 10;
    TLorentzVector jJet = TLorentzVector();
    for (int jet = 0; jet < inJets.size();jet++)
    {
	jJet.SetPtEtaPhiM(inJets.at(jet).Pt(),0,inJets.at(jet).Phi(),0);
	mht -= jJet;
    }
    for (int jet = 0; jet < inJets.size();jet++)
    {
	jJet.SetPtEtaPhiM(inJets.at(jet).Pt(),0,inJets.at(jet).Phi(),0);
	TLorentzVector biasedMht = mht + jJet;
	double biasedDPhi = fabs(biasedMht.DeltaPhi(jJet));
	if (fabs(biasedDPhi) < biasedDPhiMin) biasedDPhiMin = fabs(biasedDPhi);
    }
    return biasedDPhiMin;

}
void AlphaT20b::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

    //std::cout << "entries: " << treereader.GetEntries() << std::endl;

    andir->cd();

    mCounter+=(weight); //keep a tally of all the files/events we are running over

    //produce subarrays of objects satisfying our criteria
    std::vector<jlepton> ele=goodleptons(treereader->GetElec(), 10.0, 2.5,2.5,2.6);         //the central isolated electrons, pt > PT_ELEC GeV
    std::vector<jlepton> mu=goodleptons(treereader->GetMuon(), 10.0, 2.5,2.5,2.6);          //the central isolated muons, pt > PT_MUON GeV
    std::vector<jphoton> photon=goodphotons(treereader->GetPhoton(), 25.0, 2.5,2.5,2.6);
    std::vector<jtrack> track=trackSkim(treereader->GetIsoChargedTrack(), 10.0);
    std::vector<jtrack> trackTest=trackSkim(treereader->GetIsoChargedTrack(), 0.0);
    std::vector<jjet> badjets=badjetsSkim(treereader->GetJet(),50.0,3.0); //check for jets which we should veto
    std::vector<jjet> goodjets200=goodjetsSkim(treereader->GetJet(),36.7,3.0); //check for jets which we should analyse
    std::vector<jjet> goodjets275=goodjetsSkim(treereader->GetJet(),36.7,3.0); //check for jets which we should analyse
    std::vector<jjet> goodjets325=goodjetsSkim(treereader->GetJet(),43.3,3.0); //check for jets which we should analyse
    std::vector<jjet> goodjets375=goodjetsSkim(treereader->GetJet(),50.0,3.0); //check for jets which we should analyse
    std::vector<jjet> deadEcalJets=goodjetsSkim(treereader->GetJet(),30.0,3.0); //check for jets which we should analyse
    std::vector<jjet> etmis=treereader->GetMet(); //Missing transverse energy array
    std::vector<jparticle> particles = gentreereader->GetGenParticle(); 
    std::vector<jecal> ecalMap = getMapEcal();
    // for (int ecal = 0; ecal < ecalMap.size(); ecal++)
    // {
    // ecal_map->Fill(ecalMap[ecal].Phi(),ecalMap[ecal].Eta());
    // }
    //std::cout << ecalCut << std::endl;
    //bool foundFist = false;
    TLorentzVector stop1;
    TLorentzVector astop1;
    for (int part = 0; part < particles.size(); part++)
    {
	if (particles[part].PID() == 1000006) stop1.SetPtEtaPhiM(particles[part].Pt(),particles[part].Eta(),particles[part].Phi(),0);
	else if (particles[part].PID() == -1000006) astop1.SetPtEtaPhiM(particles[part].Pt(),particles[part].Eta(),particles[part].Phi(),0);
    }
    TLorentzVector boost = stop1+astop1;
    boostHisto->Fill(boost.Pt());

    if (goodjets200.size() != 0)
    {
	jet1pt->Fill(goodjets200[0].Pt());
    } 
    if (goodjets200.size() > 1)
    {
	jet2pt->Fill(goodjets200[1].Pt());
    } 
    double calo_met = etmis[0].Et();
    calomethist->Fill(calo_met, weight);
    bjets->Fill(badjets.size());
    ejets->Fill(ele.size());
    mjets->Fill(mu.size());
    cut_sel->Fill(0.,weight);
    event_weight->Fill(weight);
    energy_sums esums = make_energy_sums_20(goodjets200, goodjets275, goodjets325, goodjets375);
    if (esums.nbtags == 0 && esums.njets  < 4)  cut_sel23b0->Fill(0.,weight);
    if (esums.nbtags == 0 && esums.njets  >= 4)   cut_selGe4b0->Fill(0.,weight);
    if (esums.nbtags == 1 && esums.njets  < 4)   cut_sel23b1->Fill(0.,weight);
    if (esums.nbtags == 1 && esums.njets  >= 4)   cut_selGe4b1->Fill(0.,weight);

    if(badjets.size() == 0 && ele.size() == 0 && mu.size() == 0) {
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
			double bDPhi = biasedDPhiMin(goodjets375); 
			bool ecalCut = deadEcalCut(ecalMap,deadEcalJets);
			if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At53->Fill(bDPhi,weight);
			if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At53->Fill(bDPhi,weight);
			if (esums.nbtags == 0)
			{
			    if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At55b0->Fill(bDPhi,weight);
			    if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At55b0->Fill(bDPhi,weight);
			    if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At53b0->Fill(bDPhi,weight);
			    if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At53b0->Fill(bDPhi,weight);
			}
			if (esums.nbtags == 1)
			{
			    if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At55b1->Fill(bDPhi,weight);
			    if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At55b1->Fill(bDPhi,weight);
			    if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At53b1->Fill(bDPhi,weight);
			    if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At53b1->Fill(bDPhi,weight);
			}
			if(ecalCut)
			{
			    if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At53EcalCut->Fill(bDPhi,weight);
			    if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At53EcalCut->Fill(bDPhi,weight);
			    if (esums.nbtags == 0)
			    {
				if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At55b0EcalCut->Fill(bDPhi,weight);
				if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At55b0EcalCut->Fill(bDPhi,weight);
				if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At53b0EcalCut->Fill(bDPhi,weight);
				if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At53b0EcalCut->Fill(bDPhi,weight);
			    }
			    if (esums.nbtags == 1)
			    {
				if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At55b1EcalCut->Fill(bDPhi,weight);
				if(alpha_t > 0.55 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At55b1EcalCut->Fill(bDPhi,weight);
				if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets >= 4 && goodjets375[1].Pt() > 100) bDPhiHistoGe4At53b1EcalCut->Fill(bDPhi,weight);
				if(alpha_t > 0.53 && esums.total_ht > 375.0 && esums.njets < 4 && goodjets375[1].Pt() > 100) bDPhiHisto23At53b1EcalCut->Fill(bDPhi,weight);
			    }
			}
			else
			{
			    for (int jet = 0; jet < deadEcalJets.size(); jet++)
			    {
				ecal_map->Fill(deadEcalJets[jet].Phi(),deadEcalJets[jet].Eta());
			    }
			}
			if(alpha_t > 0.55) {
			    cut_sel->Fill(6.,weight);
			    if (esums.nbtags == 0 && esums.njets  < 4)  cut_sel23b0->Fill(6.,weight);
			    if (esums.nbtags == 0 && esums.njets  >= 4)   cut_selGe4b0->Fill(6.,weight);
			    if (esums.nbtags == 1 && esums.njets  < 4)   cut_sel23b1->Fill(6.,weight);
			    if (esums.nbtags == 1 && esums.njets  >= 4)   cut_selGe4b1->Fill(6.,weight);
			    if(ecalCut)
			    {
				cut_sel->Fill(7.,weight);

				if (esums.nbtags == 0 && esums.njets  < 4)  cut_sel23b0->Fill(7.,weight);
				if (esums.nbtags == 0 && esums.njets  >= 4)   cut_selGe4b0->Fill(7.,weight);
				if (esums.nbtags == 1 && esums.njets  < 4)   cut_sel23b1->Fill(7.,weight);
				if (esums.nbtags == 1 && esums.njets  >= 4)   cut_selGe4b1->Fill(7.,weight);
				if(trackTest.size() != 0) cTrack->Fill(trackTest.at(0).Pt());
				else {std::cout << "here" << std::endl;}
				if(track.size() == 0) {   
				    cut_sel->Fill(8.,weight);
				    //btagrate->Fill(esums.nbtags, weight);
				    ht_vs_mht_post_alphaT->Fill(esums.total_ht, esums.total_mht, weight);
				    ///////////////////////////////
				    //PSEUDO SITV HERE - MUST REMOVE
				    ///////////////////////////
				    double SITVweight = weight * 0.9;
				    if (esums.nbtags == 0 && esums.njets  < 4)  cut_sel23b0->Fill(8.,weight);
				    if (esums.nbtags == 0 && esums.njets  >= 4) cut_selGe4b0->Fill(8.,weight);
				    if (esums.nbtags == 1 && esums.njets  < 4)  cut_sel23b1->Fill(8.,weight);
				    if (esums.nbtags == 1 && esums.njets  >= 4) cut_selGe4b1->Fill(8.,weight);

				    if(esums.nbtags == 0) {
					if(esums.njets >= 2 && esums.njets <= 3){
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(0)+=(SITVweight); }
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(1)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(2)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(3)+=(SITVweight); }
					    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(4)+=(SITVweight); }
					    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(5)+=(SITVweight); }
					    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(6)+=(SITVweight); }
					    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(7)+=(SITVweight); }
					    if(esums.total_ht > 875.0 && esums.total_ht <= 975.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(8)+=(SITVweight); }
					    if(esums.total_ht > 975.0 && esums.total_ht <= 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(9)+=(SITVweight); }
					    if(esums.total_ht > 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(10)+=(SITVweight); }
					}

					if(esums.njets >= 4){
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(11)+=(SITVweight); }
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(12)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(13)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(14)+=(SITVweight); }
					    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(15)+=(SITVweight); }
					    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(16)+=(SITVweight); }
					    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(17)+=(SITVweight); }
					    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(18)+=(SITVweight); }
					    if(esums.total_ht > 875.0 && esums.total_ht <= 975.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(19)+=(SITVweight); }
					    if(esums.total_ht > 975.0 && esums.total_ht <= 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(20)+=(SITVweight); }
					    if(esums.total_ht > 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(21)+=(SITVweight); }
					}

				    }
				    else if(esums.nbtags == 1) {
					if(esums.njets >= 2 && esums.njets <= 3){
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(22)+=(SITVweight); }
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(23)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(24)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(25)+=(SITVweight); }
					    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(26)+=(SITVweight); }
					    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(27)+=(SITVweight); }
					    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(28)+=(SITVweight); }
					    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(29)+=(SITVweight); }
					    if(esums.total_ht > 875.0 && esums.total_ht <= 975.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(30)+=(SITVweight); }
					    if(esums.total_ht > 975.0 && esums.total_ht <= 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(31)+=(SITVweight); }
					    if(esums.total_ht > 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(32)+=(SITVweight); }
					}
					if(esums.njets >= 4){
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(33)+=(SITVweight); }
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(34)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(35)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(36)+=(SITVweight); }
					    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(37)+=(SITVweight); }
					    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(38)+=(SITVweight); }
					    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(39)+=(SITVweight); }
					    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(40)+=(SITVweight); }
					    if(esums.total_ht > 875.0 && esums.total_ht <= 975.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(41)+=(SITVweight); }
					    if(esums.total_ht > 975.0 && esums.total_ht <= 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(42)+=(SITVweight); }
					    if(esums.total_ht > 1075.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(43)+=(SITVweight); }
					}

				    }

				    else if(esums.nbtags == 2) {
					if(esums.njets >= 2 && esums.njets <= 3){
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(44)+=(SITVweight); }
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(45)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(46)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(47)+=(SITVweight); }
					    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(48)+=(SITVweight); }
					    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(49)+=(SITVweight); }
					    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(50)+=(SITVweight); }
					    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(51)+=(SITVweight); }
					    if(esums.total_ht > 875.0 && esums.total_ht <= 975.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(52)+=(SITVweight); }
					}
					if(esums.njets >= 4){
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(53)+=(SITVweight); }
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(54)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(55)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(56)+=(SITVweight); }
					    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(57)+=(SITVweight); }
					    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(58)+=(SITVweight); }
					    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(59)+=(SITVweight); }
					    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(60)+=(SITVweight); }
					    if(esums.total_ht > 875.0 && esums.total_ht <= 975.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(61)+=(SITVweight); }
					}

				    } else if(esums.nbtags == 3) {
					if(esums.njets >= 4){
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(62)+=(SITVweight); } 
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(63)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(64)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(65)+=(SITVweight); }
					    if(esums.total_ht > 475.0 && esums.total_ht <= 575.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(66)+=(SITVweight); }
					    if(esums.total_ht > 575.0 && esums.total_ht <= 675.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(67)+=(SITVweight); }
					    if(esums.total_ht > 675.0 && esums.total_ht <= 775.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(68)+=(SITVweight); }
					    if(esums.total_ht > 775.0 && esums.total_ht <= 875.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(69)+=(SITVweight); }
					    if(esums.total_ht > 875.0 && esums.total_ht <= 975.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(70)+=(SITVweight); }
					}
				    }
				    else if(esums.nbtags >= 4) {
					if(esums.njets >= 4){ 
					    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(71)+=(SITVweight); }
					    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(72)+=(SITVweight); }
					    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.60 && goodjets325[1].Pt() > 86.7) { mSigPred.at(73)+=(SITVweight); }
					    if(esums.total_ht > 375.0 && esums.total_ht <= 475.0 && goodjets375[1].Pt() > 100.0) { mSigPred.at(74)+=(SITVweight); }
					    //if(esums.total_ht > 475.0 && esums.total_ht <= 575.0) { mSigPred.at(59)+=(SITVweight); }
					    //if(esums.total_ht > 575.0 && esums.total_ht <= 675.0) { mSigPred.at(60)+=(SITVweight); }
					    //if(esums.total_ht > 675.0 && esums.total_ht <= 775.0) { mSigPred.at(61)+=(SITVweight); }
					    //if(esums.total_ht > 775.0 && esums.total_ht <= 875.0) { mSigPred.at(62)+=(SITVweight); }
					    //if(esums.total_ht > 875.0) { mSigPred.at(63)+=(SITVweight); }
					}
				    }
				}

			    }
			}
		    }
		}
	    }
	}
    }
    return;
    //}
    }
