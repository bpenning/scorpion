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
}

bool reject_OSSF(std::vector<jlepton> leptons){
  bool result=false;
  std::vector<jlepton>::iterator l1_it;
  std::vector<jlepton>::iterator l2_it;
  std::vector<jlepton>::iterator l3_it;
  for (l1_it=leptons.begin();l1_it!=leptons.end();l1_it++){
     for (l2_it=l1_it+1;l2_it!=leptons.end();l2_it++){
         //if opposite sign same flavour
         if (l1_it->Charge()==-l1_it->Charge() && l1_it->Flavour()==l2_it->Flavour()){
            for (l3_it=leptons.begin();l3_it!=leptons.end();l3_it++){
                if ((l3_it!=l1_it) && (l3_it!=l2_it)){
                    double mlll=((*l1_it)+(*l2_it)+(*l3_it)).M(); 
                    if (mlll > 75 && mlll < 115){
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
  std::vector<std::pair<jlepton,jlepton> > OSSF_pairs=get_OSSF_pairs(leptons);

  bool selected=false;
  //at least three reconstructed leptons, where by "lepton" we mean an electron, muon or tau candidate
  if (leptons.size()+taujets.size()>=3){
    //at least one electron or muon candidate must satisfy pt>20 GeV
    if (leptons[0].Pt()>20){
      //at most one tau candidate
      if (taujets.size()<=1){
         //"Events with an OSSF pair outside the Z boson mass region (75<m_OSSF<115), but that satisfy 
         //75<m_(OSSF+lepton)<115 are liekly te arise from final-state radiation from Z-boson decay producs,
         //followed by conversion of the photon to a charged lepton pair. Events that meet this condition are
         //rejected if they also exhibit kinematic characteristics consistent with background from events with a Z
         //boson and jets" (p3 CMS-SUS-13-003)
         //FIXME: can the third lepton be a tau?
         if (!reject_OSSF(leptons)){
            selected=true;
         }
      }
    }
  }

  if (selected){
    //(kinematic) variables
    double met=etmis[0].E();
    //FIXME: correct this
    int N_tau=0;
    int N_b=0;
    //definition of signal regions: 
    //Table 2 CMS-SUS-13-002 divisions (HT>200 GeV; HT<200 GeV); rows (OSSF0 ... OSSF2); columns (N_tau=0,N_b=0 ... N_tau=1, N_b>=1)
    if (leptons.size()+taujets.size()>=4){
       if (OSSF_pairs.size()==0){
           if (met>100){
               if ((N_tau==0) && (N_b==0)) mSigPred.at(0)+=weight;
               else if ((N_tau==1) && (N_b==0)) mSigPred.at(1)+=weight;
               else if ((N_tau==0) && (N_b>=1)) mSigPred.at(2)+=weight;
               else if ((N_tau==1) && (N_b>=1)) mSigPred.at(3)+=weight;
           } else if (50<met && met<100){
               if ((N_tau==0) && (N_b==0)) mSigPred.at(5)+=weight;
               else if ((N_tau==1) && (N_b==0)) mSigPred.at(6)+=weight;
               else if ((N_tau==0) && (N_b>=1)) mSigPred.at(7)+=weight;
               else if ((N_tau==1) && (N_b>=1)) mSigPred.at(8)+=weight;
           }
           } else if (0<met && met<50){
               if ((N_tau==0) && (N_b==0)) mSigPred.at(9)+=weight;
               else if ((N_tau==1) && (N_b==0)) mSigPred.at(10)+=weight;
               else if ((N_tau==0) && (N_b>=1)) mSigPred.at(11)+=weight;
               else if ((N_tau==1) && (N_b>=1)) mSigPred.at(12)+=weight;
           }
       
       }else if (OSSF_pairs.size()==1){
       
       }else if (OSSF_pairs.size()==2){
       }
      
    }
  
  return;
}
