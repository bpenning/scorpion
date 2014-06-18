#include "analysis_os5.hh"

CmsOs5Fb::CmsOs5Fb(const std::string & name, 
		 const std::string & experiment,
		 const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

CmsOs5Fb::CmsOs5Fb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 //const std::vector<int> & datayields,
		 const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

CmsOs5Fb::CmsOs5Fb(const std::string & name, 
		 const std::string & experiment, 
		 const unsigned int & numBins,
		 const double & intlumi, 
		 const std::vector<double> & bgpred,
		 const std::vector<double> & bgpreduncert,
		 const std::vector<int> & datayields,
		 const std::string & fitmode,
		 const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}


CmsOs5Fb::~CmsOs5Fb() {}

void CmsOs5Fb::initHistos() {
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

void CmsOs5Fb::Run(const Reader * treereader, const Reader * gentreereader, 
        const double & weight) {

  andir->cd();

  mCounter+=weight; //keep a tally of all the files/events we are running over

  //produce subarrays of objects satisfying our criteria
  std::vector<jlepton> poslep = leptonChargeSkim(treereader->GetElec(), 
          treereader->GetMuon(), 10.0, 2.5, 10.0, 2.4, 1.4442,1.566,1);
  std::vector<jlepton> neglep = leptonChargeSkim(treereader->GetElec(), 
          treereader->GetMuon(), 10.0, 2.5, 10.0, 2.4, 1.4442,1.566,-1);
  std::vector<jlepton> leptons = leptonSkim(treereader->GetElec(), 
          treereader->GetMuon(), 10.0, 2.5, 10.0, 2.4, 1.4442,1.566);
  jjet etmis=(treereader->GetMet())[0]; 

  std::vector<jjet> goodjets = goodjetsSkimDRcut(treereader->GetJet(), 30.0, 
          3.0, leptons, 0.4);

  if(goodjets.size() >= 2 ) {
    double HT = 1.05 * getht(goodjets);
    double MET = etmis.Et();
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
