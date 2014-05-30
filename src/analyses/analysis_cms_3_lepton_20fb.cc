#include "analysis_cms_3_lepton_20fb.hh"

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
Cms3Lepton20Fb::Cms3Lepton20Fb(const std::string & name, 
		   const std::string & experiment,
		   const unsigned int & numBins) :
  AnalysisBase(name, experiment, numBins) {}

//constructor for object where you also want to run a limit and so other information
//is necessary
Cms3Lepton20Fb::Cms3Lepton20Fb(const std::string & name, 
		   const std::string & experiment, 
		   const unsigned int & numBins,
		   const double & intlumi, 
		   const std::vector<double> & bgpred) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred) {}

Cms3Lepton20Fb::Cms3Lepton20Fb(const std::string & name, 
		   const std::string & experiment, 
		   const unsigned int & numBins,
		   const double & intlumi, 
		   const std::vector<double> & bgpred,
		   const std::vector<double> & bgpreduncert,
		   const std::vector<int> & datayields,
		   const std::string & fitmode,
		   const bool & calculateR) :
  AnalysisBase(name, experiment, numBins, intlumi, bgpred, bgpreduncert, datayields, fitmode, calculateR) {}

Cms3Lepton20Fb::~Cms3Lepton20Fb() {}

void Cms3Lepton20Fb::initHistos() {
  andir->cd();
  cut_flow_hist = new TH1D("cuts",";cuts;entries",7,0,7);
  nleptons = new TH1D("nleptons", ";N_{jets};Entries",10,-0.5,9.5);
  njets = new TH1D("njets", ";N_{jets};Entries",10,-0.5,9.5);
  nbjets = new TH1D("nbjets", ";N_{bjets};Entries",10,-0.5,9.5);
  nossfpairs = new TH1D("nossfpairs", ";N_{ossfpairs};Entries",10,-0.5,9.5);
  ntaujets = new TH1D("ntaujets", ";N_{taujets};Entries",10,-0.5,9.5);
  met_hist= new TH1D("met",";E_{T}^{miss} [GeV];Entries",10 ,0 ,500);
  ht_hist= new TH1D("ht",";H_{T} [GeV];Entries",10 ,100 ,350);
}

ossf_bools Cms3Lepton20Fb::get_ossf_bools(const std::vector<std::pair<jlepton,jlepton> > & ossf_pairs, double m_ossf_min, double m_ossf_max){
    ossf_bools result;
    result.on_Z=false; 
    double m_ossf;
    std::vector< std::pair<jlepton, jlepton> >::const_iterator ossf_it;
    for (ossf_it=ossf_pairs.begin();ossf_it!=ossf_pairs.end();ossf_it++){
        m_ossf=(ossf_it->first+ossf_it->second).M();
        if (m_ossf_min<m_ossf && m_ossf<m_ossf_max){
            result.on_Z=true;
            break;
        }
    }
    result.below_Z=m_ossf<m_ossf_min;
    result.above_Z=m_ossf>m_ossf_max;
    return result;
}

bool Cms3Lepton20Fb::reject_ossf(const std::vector<jlepton> & leptons, double mlll_min, double mlll_max){
  bool result=false;
  std::vector<jlepton>::const_iterator l1_it;
  std::vector<jlepton>::const_iterator l2_it;
  std::vector<jlepton>::const_iterator l3_it;
  for (l1_it=leptons.begin();l1_it!=leptons.end();l1_it++){
     for (l2_it=l1_it+1;l2_it!=leptons.end();l2_it++){
         //if opposite sign same flavour
         if (l1_it->Charge()==-l2_it->Charge() && l1_it->Flavour()==l2_it->Flavour()){
            for (l3_it=leptons.begin();l3_it!=leptons.end();l3_it++){
                if ((l3_it!=l1_it) && (l3_it!=l2_it)){
                    double mlll=((*l1_it)+(*l2_it)+(*l3_it)).M(); 
                    if (mlll > mlll_min && mlll < mlll_max){
                        //FIXME: kinematic characteristics consistend with background from events with Z boson and jets (Z+jets)
                        if (true){
                            result=true;
                        }
                    }
                }  
            }
         }
     }
  }
  return result;
}

void Cms3Lepton20Fb::Run(const Reader * treereader, const Reader * gentreereader, const double & weight) {
  andir->cd();

  mCounter+=weight; 

  //get electrons, muons, taujets and jets
  //"electron and muon candidates must satisfy pt>10GeV and |eta|<2.4"
  std::vector<jlepton> leptons=leptonSkim(treereader->GetElec(),treereader->GetMuon(), 10.0, 2.4,10.0,2.4);
  //tau candidates must satisfy pt>20 GeV and |eta|<2.3
  std::vector<jjet> taujets=goodjetsSkim(treereader->GetTauJet(),20,2.3); 
  std::vector<jjet> jets=goodjetsSkimDRcut(treereader->GetJet(),30,2.5,leptons,0.3); 
  std::vector<jjet> etmis=treereader->GetMet(); 
  std::vector<std::pair<jlepton,jlepton> > ossf_pairs=get_ossf_pairs(leptons);

  bool selected=false;
  cut_flow_hist->Fill(0.5);
  //at least three reconstructed leptons, where by "lepton" we mean an electron, muon or tau candidate
  if (leptons.size()+taujets.size()>=3){
    cut_flow_hist->Fill(1.5);
    //at most one tau candidate
    if (taujets.size()<=1){
      cut_flow_hist->Fill(2.5);
      //at least one electron or muon candidate must satisfy pt>20 GeV
      if (leptons[0].Pt()>20){
         cut_flow_hist->Fill(3.5);
         //FIXME: can the third lepton be a tau?
         if (!reject_ossf(leptons,75,115)){
            cut_flow_hist->Fill(4.5);
            selected=true;
         }
      }
    }
  }

  if (selected){
    //(kinematic) variables
    double met=etmis[0].E();
    double ht=getht(jets);
    int N_tau=taujets.size();
    int N_b=getnbtags(jets);
    ossf_bools bools=get_ossf_bools(ossf_pairs,75,105);
    bool on_Z=bools.on_Z;
    bool above_Z=bools.above_Z;
    bool below_Z=bools.below_Z;
    // plot some of these variables

    njets->Fill(jets.size());
    nleptons->Fill(leptons.size()+taujets.size());
    nbjets->Fill(N_b);
    nossfpairs->Fill(ossf_pairs.size());
    ntaujets->Fill(N_tau);
    met_hist->Fill(met);
    ht_hist->Fill(ht);
    //definition of signal regions: 
    //Table 2 CMS-SUS-13-002 divisions (HT>200 GeV; HT<200 GeV); rows (OSSF0 ... OSSF2); columns (N_tau=0,N_b=0 ... N_tau=1, N_b>=1)
    if (leptons.size()+taujets.size()>=4){
       if (ht>200){
           if (ossf_pairs.size()==0 ){
               if (met>100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(0)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(1)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(2)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(3)+=weight;
               } else if (50<met && met<100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(4)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(5)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(6)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(7)+=weight;
               } else if (0<met && met<50){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(8)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(9)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(10)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(11)+=weight;
               }
           }else if (ossf_pairs.size()==1){
               if (met>100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(12)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(13)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(14)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(15)+=weight;
               } else if (met>100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(16)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(17)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(18)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(19)+=weight;
               } else if (50<met && met<100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(20)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(21)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(22)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(23)+=weight;
               } else if (50<met && met<100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(24)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(25)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(26)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(27)+=weight;
               } else if (0<met && met<50 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(28)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(29)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(30)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(31)+=weight;
               } else if (0<met && met<50 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(32)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(33)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(34)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(35)+=weight;
               }
           }else if (ossf_pairs.size()==2){
               if (met>100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(36)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(37)+=weight;
               } else if (met>100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(38)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(39)+=weight;
               } else if (50<met && met<100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(40)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(41)+=weight;
               } else if (50<met && met<100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(42)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(43)+=weight;
               } else if (0<met && met<50 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(44)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(45)+=weight;
               } else if (0<met && met<50 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(46)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(47)+=weight;
               }
           }
    }else{
           if (ossf_pairs.size()==0 ){
               if (met>100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(48)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(49)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(50)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(51)+=weight;
               } else if (50<met && met<100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(52)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(53)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(54)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(55)+=weight;
               } else if (0<met && met<50){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(56)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(57)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(58)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(59)+=weight;
               }
           }else if (ossf_pairs.size()==1){
               if (met>100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(60)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(61)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(62)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(63)+=weight;
               } else if (met>100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(64)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(65)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(66)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(67)+=weight;
               } else if (50<met && met<100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(68)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(69)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(70)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(71)+=weight;
               } else if (50<met && met<100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(72)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(73)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(74)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(75)+=weight;
               } else if (0<met && met<50 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(76)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(77)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(78)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(79)+=weight;
               } else if (0<met && met<50 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(80)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(81)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(82)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(83)+=weight;
               }
           
           }else if (ossf_pairs.size()==2){
               if (met>100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(84)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(85)+=weight;
               } else if (met>100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(86)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(87)+=weight;
               } else if (50<met && met<100 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(88)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(89)+=weight;
               } else if (50<met && met<100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(90)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(91)+=weight;
               } else if (0<met && met<50 && on_Z==false){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(92)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(93)+=weight;
               } else if (0<met && met<50 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(94)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(95)+=weight;
               }
           }
        }
     } else {
       if (ht>200){
           if (ossf_pairs.size()==0 ){
               if (met>100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(96)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(97)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(98)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(99)+=weight;
               } else if (50<met && met<100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(100)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(101)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(102)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(103)+=weight;
               } else if (0<met && met<50){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(104)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(105)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(106)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(107)+=weight;
               }
           }else if (ossf_pairs.size()==1){
               if (met>100 && above_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(108)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(109)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(110)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(111)+=weight;
               } else if (met>100 && below_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(112)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(113)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(114)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(115)+=weight;
               } else if (met>100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(116)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(117)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(118)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(119)+=weight;
               } else if (50<met && met<100 && above_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(120)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(121)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(122)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(123)+=weight;
               } else if (50<met && met<100 && below_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(124)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(125)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(126)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(127)+=weight;
               } else if (50<met && met<100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(128)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(129)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(130)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(131)+=weight;
               } else if (0<met && met<50 && above_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(132)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(133)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(134)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(135)+=weight;
               } else if (0<met && met<50 && below_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(136)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(137)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(138)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(139)+=weight;
               } else if (0<met && met<50 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(140)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(141)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(142)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(143)+=weight;
               }
           }
    }else{
           if (ossf_pairs.size()==0 ){
               if (met>100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(144)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(145)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(146)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(147)+=weight;
               } else if (50<met && met<100){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(148)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(149)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(150)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(151)+=weight;
               } else if (0<met && met<50){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(152)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(153)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(154)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(155)+=weight;
               }
           }else if (ossf_pairs.size()==1){
               if (met>100 && above_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(156)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(157)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(158)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(159)+=weight;
               } else if (met>100 && below_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(160)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(161)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(162)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(163)+=weight;
               } else if (met>100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(164)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(165)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(166)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(167)+=weight;
               } else if (50<met && met<100 && above_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(168)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(169)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(170)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(171)+=weight;
               } else if (50<met && met<100 && below_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(172)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(173)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(174)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(175)+=weight;
               } else if (50<met && met<100 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(176)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(177)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(178)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(179)+=weight;
               } else if (0<met && met<50 && above_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(180)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(181)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(182)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(183)+=weight;
               } else if (0<met && met<50 && below_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(184)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(185)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(186)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(187)+=weight;
               } else if (0<met && met<50 && on_Z==true){
                   if ((N_tau==0) && (N_b==0)) mSigPred.at(188)+=weight;
                   else if ((N_tau==1) && (N_b==0)) mSigPred.at(189)+=weight;
                   else if ((N_tau==0) && (N_b>=1)) mSigPred.at(190)+=weight;
                   else if ((N_tau==1) && (N_b>=1)) mSigPred.at(191)+=weight;
               }
           }
        }
     }
  }
  return;
}
