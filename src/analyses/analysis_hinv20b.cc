#include "analysis_hinv20b.hh"


//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
Hinv20b::Hinv20b(const std::string & name, 
		 const std::string & experiment,
		 const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
Hinv20b::Hinv20b(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 //const std::vector<int> & datayields,
		 const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

Hinv20b::Hinv20b(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 const std::vector<double> & bgpred,
		 const std::vector<double> & bgpreduncert,
		 const std::vector<int> & datayields,
		 const std::string & fitmode,
		 const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


Hinv20b::~Hinv20b() {}



void Hinv20b::initHistos() {
  andir->cd();
  //!!MAKE SOME HISTOGRAMS: TAILOR TO HINV ONES
  event_weight = new TH1D("event_weight",";weight_num;entries",1000,1940e9,1960e9);
  jet1pt = new TH1D("firstjetpt", ";P_{T} [GeV];Entries",80,-5.,395.);
  jet2pt = new TH1D("secondjetpt", ";P_{T} [GeV];Entries",60,-5.,295.);
  jet1eta = new TH1D("firstjeteta", ";P_{T} [GeV];Entries",50,-5.,5.);
  jet2eta = new TH1D("secondjeteta", ";P_{T} [GeV];Entries",50,-5.,5.);
  jet1phi = new TH1D("firstjetphi", ";P_{T} [GeV];Entries",16,-3.2,3.2);
  jet2phi = new TH1D("secondjetphi", ";P_{T} [GeV];Entries",16,-3.2,3.2);
  jetmet_mindphi = new TH1D("jetmetmindphi",";;",64,0.,3.2);
  metsignificance = new TH1D("metsignificance",";;",150,0,30);
  met = new TH1D("met",";;",40,0,400.);
  deltaphijj = new TH1D("deltaphijj",";;",72,0.,3.55);
  deltaetajj = new TH1D("deltaetajj",";;",50,0.,10.);
  mjj = new TH1D("mjj",";;",20,0.,2000.);
  njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
  nelectrons = new TH1D("nelectrons", ";N_{jets};Entries",10,-0.5,9.5);
  nmuons = new TH1D("nmuons", ";N_{jets};Entries",10,-0.5,9.5);
  event_weight->SetBit(TH1::kCanRebin);
}

void Hinv20b::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {

  //std::cout << "entries: " << treereader->GetEntries() << std::endl;

  andir->cd();

  mCounter+=(weight); //keep a tally of all the files/events we are running over

  //!!get all the objects, check these are compatible with our definitions
  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> ele=goodleptons(treereader->GetElec(), 10.0, 2.4,1.44,1.57);       //the central isolated electrons, pt > PT_ELEC GeV
  std::vector<jlepton> mu=goodleptons(treereader->GetMuon(), 10.0, 2.1);          //the central isolated muons, pt > PT_MUON GeV
  //!!  std::vector<jlepton> tau=goodleptons(treereader->GetTau(), 20.0, 2.3);          //the central isolated muons, pt > PT_MUON GeV
  std::vector<jjet> vbfjets=goodjetsSkim(treereader->GetJet(),30,4.7);
  std::vector<jjet> triggeremulationjets=goodjetsSkim(treereader->GetJet(),0.,3);
  std::vector<jjet> etmis=treereader->GetMet(); //Missing transverse energy array
  std::vector<double> sumet=treereader->GetSumET(); //scalar sum of transverse energy


  //emulate trigger with mht from jets with eta<3
  double etx=0;
  double ety=0;
  for(unsigned iJet=0;iJet<triggeremulationjets.size();iJet++){
    etx+=triggeremulationjets[iJet].Px();  
    ety+=triggeremulationjets[iJet].Py();  
  }
  double l1met=sqrt(pow(etx,2)+pow(ety,2));

  //Get met significance 
  double met_significance=etmis[0].Et()/sqrt(sumet[0]);
  
  //Loop over jets to find jet 1 and 2 and do jetmetdphi
  double jet1_pt=-1;
  double jet2_pt=-1;
  double jet1_eta=-10000;
  double jet2_eta=-10000;
  double jet1_phi=-1;
  double jet2_phi=-1;
  double jet1_E=-1;
  double jet2_E=-1;
  double jetmetmindphi=10;
  int ijet1=-1,ijet2=-1;
  //std::cout<<"Listing  jet pts:"<<std::endl;                                                                                                        
  for(unsigned iJet=0;iJet<vbfjets.size();iJet++){
    //Find two highest pt jets
    double pt=vbfjets[iJet].Pt();//!!fix for jad framework
    if(pt>jet2_pt){
      jet2_pt=pt;
      ijet2=iJet;
      if(jet2_pt>jet1_pt){
	double tmppt=jet1_pt;
	double tmpijet=ijet1;
	jet1_pt=jet2_pt;
	ijet1=ijet2;
	jet2_pt=tmppt;
	ijet2=tmpijet;
      }
    }

    //Do jetmetdphi calculation
    double thisjetmetdphi=fabs(fabs(fabs(vbfjets[iJet].Phi()-etmis[0].Phi())-3.14159265358979)-3.14159265358979);  
    if(thisjetmetdphi>3.141593)std::cout<<"Warning: Delta phi greater than pi jphi "<<vbfjets[iJet].Phi()<<" met phi "<<etmis[0].Phi()<<std::endl;
    if(thisjetmetdphi<jetmetmindphi){
      jetmetmindphi=thisjetmetdphi;
    }
  }

  //Get first two jets information
  if(vbfjets.size()>=1){
    jet1_eta=vbfjets[ijet1].Eta();
    jet1_phi=vbfjets[ijet1].Phi();
    jet1_E=vbfjets[ijet1].E();
  }
  if(vbfjets.size()>=2){
    jet2_eta=vbfjets[ijet2].Eta();
    jet2_phi=vbfjets[ijet2].Phi();
    jet2_E=vbfjets[ijet2].E();
  }
  int njetscjv=0;
  if(vbfjets.size()>=3){
    for(unsigned iJet=0;iJet<vbfjets.size();iJet++){
      if(jet1_eta>0){
	if((iJet!=ijet1)&&(iJet!=ijet2)&&(vbfjets[iJet].Eta()<jet1_eta)&&(vbfjets[iJet].Eta()>jet2_eta)) njetscjv++;
      }
      else{
	if((iJet!=ijet1)&&(iJet!=ijet2)&&(vbfjets[iJet].Eta()>jet1_eta)&&(vbfjets[iJet].Eta()<jet2_eta)) njetscjv++;
      }
    }
  }

  //Get dijet_deta, dijet_M and dijet_dphi
  double dijet_deta=-1;
  double dijet_dphi=-1;
  double dijet_M=-1;
  if(vbfjets.size()>=2){
    dijet_deta = fabs(vbfjets[ijet1].Eta() - vbfjets[ijet2].Eta());
    dijet_dphi = fmod(fabs(vbfjets[ijet1].Phi()-vbfjets[ijet2].Phi()),3.141593);  
    if(dijet_dphi>3.141593)std::cout<<"Warning: Delta phi greater than pi j1phi "<<vbfjets[ijet1].Phi()<<" j2phi "<<vbfjets[ijet2].Phi()<<std::endl;
    dijet_M = (vbfjets[ijet1]+vbfjets[ijet2]).M();
  }
  if(vbfjets.size()>=2){
    //!!Fill hinv histograms and set msigpred to right value
    //if((jet1_eta*jet2_eta<0)&&jet1_eta<4.7&&jet2_eta<4.74&&dijet_deta>4.2&&jet1_pt>50&&jet2_pt>50&&etmis[0].Et()>130&&dijet_M>1100&&dijet_dphi<1.0&&njetscjv==0){//prompt analysis

    //if((jet1_eta*jet2_eta<0)&&jet1_eta<4.7&&jet2_eta<4.7&&met_significance>3&&dijet_deta>3.6&&jet1_pt>50&&jet2_pt>45&&jetmetmindphi>2.3&&etmis[0].Et()>90){//!!relaxed mjj and metsig cuts
    if((jet1_eta*jet2_eta<0)&&jet1_eta<4.7&&jet2_eta<4.7&&met_significance>4&&dijet_deta>3.6&&jet1_pt>50&&jet2_pt>45&&jetmetmindphi>2.3&&etmis[0].Et()>90){
      //std::cout<<"filling"<<std::endl;
      event_weight->Fill(weight);
      jet1pt->Fill(jet1_pt);
      jet2pt->Fill(jet2_pt);
      jet1eta->Fill(jet1_eta);
      jet2eta->Fill(jet2_eta);
      jet1phi->Fill(jet1_phi);
      jet2phi->Fill(jet2_phi);
      jetmet_mindphi->Fill(jetmetmindphi);
      metsignificance->Fill(met_significance);
      met->Fill(etmis[0].Et());
      deltaphijj->Fill(dijet_dphi);
      deltaetajj->Fill(dijet_deta);
      mjj->Fill(dijet_M);
      njets->Fill(vbfjets.size());
      nelectrons->Fill(ele.size());
      nmuons->Fill(mu.size());
      //}
      //if(jet1_pt>50&&jet2_pt>45&&(jet1_eta*jet2_eta<=0)&&jetmetmindphi>2.3&&met_significance>4&&etmis[0].Et()>90&&dijet_deta>3.6&&dijet_M>1200&&ele.size()==0&&mu.size()==0){
      mSigPred.at(0)+=weight;
    }
  }

  //!!Increment msigpred by eventweight

  //!!residual alphat things below

  /*  if (goodjets200.size() != 0)
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
	      if(alpha_t > 0.55) {
		cut_sel->Fill(6.,weight);
		//btagrate->Fill(esums.nbtags, weight);
		ht_vs_mht_post_alphaT->Fill(esums.total_ht, esums.total_mht, weight);

		///////////////////////////////
		//PSEUDO SITV HERE - MUST REMOVE
		double SITVweight = weight * 0.9;
		///////////////////////////
		if(esums.nbtags == 0) {
		  if(esums.njets >= 2 && esums.njets <= 3){
		    if(esums.total_ht > 200.0 && esums.total_ht <= 275.0 && alpha_t > 0.65 && goodjets200[1].Pt() > 73.3) { mSigPred.at(0)+=(SITVweight); }
		    if(esums.total_ht > 275.0 && esums.total_ht <= 325.0 && alpha_t > 0.60 && goodjets275[1].Pt() > 73.3) { mSigPred.at(1)+=(SITVweight); }
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(2)+=(SITVweight); }
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
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(13)+=(SITVweight); }
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
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(24)+=(SITVweight); }
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
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(35)+=(SITVweight); }
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
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(46)+=(SITVweight); }
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
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(55)+=(SITVweight); }
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
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(64)+=(SITVweight); }
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
		    if(esums.total_ht > 325.0 && esums.total_ht <= 375.0 && alpha_t > 0.55 && goodjets325[1].Pt() > 86.7) { mSigPred.at(73)+=(SITVweight); }
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
  */
  return;
  //}
}
