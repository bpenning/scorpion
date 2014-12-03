#include "analysis_base.hh"

AnalysisBase::AnalysisBase() {}

//constructor for object where you don't necessarily want to run a limit,
//but rather just run the analysis, look at distributions, see what survives
AnalysisBase::AnalysisBase(const std::string & name, 
			   const std::string & experiment,
			   const unsigned int & numBins) :
  mName(name),
  mExperiment(experiment),
  mNumBins(numBins) {
  mIntLumi = 0.0;
  mDataYields.clear();
  mBGPred.clear();
  //add another constructor with a tfile argument so the manager can call it?
  //mOutputFile = new TFile(sample.Get<std::string>("OutputFile").c_str(), "recreate");
  //aFileOut->cd();
  //TDirectory * mydir = aFileOut->mkdir(folderName, folderName);
  //mydir->cd();
}

//constructor for object where you also want to run a limit and so other information
//is necessary
AnalysisBase::AnalysisBase(const std::string & name, 
			   const std::string & experiment, 
			   const unsigned int & numBins,
			   const double & intlumi, 
			   //const std::vector<int> & datayields,
			   const std::vector<double> & bgpred) :
  mName(name),
  mExperiment(experiment),
  mNumBins(numBins),
  mIntLumi(intlumi),
  //mDataYields(datayields),
  mBGPred(bgpred) {}

//in this constructor we add the bgpreduncert which is necessary to set limits
AnalysisBase::AnalysisBase(const std::string & name, 
			   const std::string & experiment, 
			   const unsigned int & numBins,
			   const double & intlumi, 
			   const std::vector<double> & bgpred,
			   const std::vector<double> & bgpreduncert,
			   const std::vector<int> & datayields,
			   const std::string & fitmode,
			   const bool & calculateR) :
  mName(name),
  mExperiment(experiment),
  mNumBins(numBins),
  mIntLumi(intlumi),
  mBGPred(bgpred),
  mBGPredUncert(bgpreduncert),
  mDataYields(datayields),
  mFitMode(fitmode),
  mCalculateR(calculateR) {}


AnalysisBase::~AnalysisBase() {}

//make pure virtual function for initialising histograms
//MAKE SURE to initialise the TDirectory first
//virtual AnalysisBase::InitHistos()=0;

//copy constructor
AnalysisBase::AnalysisBase(const AnalysisBase & analysisobj) {
  mName = analysisobj.GetName();
  mIntLumi = analysisobj.GetLuminosity();
  mExperiment = analysisobj.GetExperiment();
  mDataYields = analysisobj.GetDataYields();
  mSigPred = analysisobj.GetSignalPrediction();
  mBGPred = analysisobj.GetBackgroundPrediction();
}

void AnalysisBase::initDir(TDirectory * file, const std::string & ifilename ) {
  std::string folderName = ifilename + "_" + this->GetExperiment() + "_" + this->GetName();
  file->cd();
  andir = file->mkdir(folderName.c_str(), folderName.c_str());
  mCurrentDir = folderName;
  return;
}

void AnalysisBase::initTree() {
  andir->cd();
  analysistree = new TTree("AnalysisTree","AnalysisTree");
  analysistree->Branch("sigpred", &mSigPred);
  analysistree->Branch("bkgpred", &mBGPred);
  analysistree->Branch("sigeff", &mSigEff);
  analysistree->Branch("sigefftot", &mSigEffTot);
  analysistree->Branch("sigbr", &mSBR);
  analysistree->Branch("sigbrtot", &mSBRTot);
}

void AnalysisBase::initLimTree() {
  andir->cd();
  limittree = new TTree("LimitTree","LimitTree");
  limittree->Branch("cls", &mCLs);
  limittree->Branch("clserr", &mCLserr);
  limittree->Branch("clb", &mCLb);
  limittree->Branch("clberr", &mCLberr);
  limittree->Branch("clsb", &mCLsb);
  limittree->Branch("clsberr", &mCLsberr);
  limittree->Branch("upperlimitR", &mUpperLim);
  limittree->Branch("upperlimitRerr", &mUpperLimerr);
}

void AnalysisBase::FillTree(const double & xsec) {
  double mtmpxsec = xsec;
  analysistree->Branch("xsec", &mtmpxsec);

  //need to do some stuff here:

  //Turn signpred into efficiency given total events
  mSigEff = mSigPred; //call the copy constructor;
  double totaleff = 0.0;
  for(std::vector<double>::iterator jj=mSigEff.begin(); jj!=mSigEff.end(); jj++) {
    (*jj)/=mCounter;
    totaleff += (*jj);
  }
  mSigEffTot = totaleff;

  //If bgvector exists, make sig/sqrt(bg) vec
  if(mBGPred.size() == mNumBins) {
    mSBR = mSigPred; //call copy constructor on yields
    std::vector<double>::iterator kk=mSBR.begin();
    std::vector<double>::const_iterator mm=mBGPred.begin();
    double tot_sig=0.0;
    double tot_bg=0.0;
    for(; kk!=mSBR.end() && mm!=mBGPred.end(); kk++, mm++) {
      tot_sig+=(*kk);
      tot_bg+=(*mm);
      (*kk)/=TMath::Sqrt((*mm));
    }
    mSBRTot = tot_sig / TMath::Sqrt(tot_bg);
    //    std::cout<<"BPBP "<<tot_sig<<" "<<mSigPred<[mNumBins]<std::endl;
  } else if(mBGPred.size() != 0) {
    std::cerr << "Warning: the background vector was not the same size as the number of bins!" << std::endl;
  } else {
    std::cerr << "Warning: no backgrounds supplied so not filling the s/sqrt(b) branch in the tree" << std::endl;
  }

  analysistree->Fill();
}

void AnalysisBase::FillLimTree(const double & cls, 
			       const double & clserr, 
			       const double & clb, 
			       const double & clberr, 
			       const double & clsb,
			       const double & clsberr,
			       const double & upperlim,
			       const double & upperlimerr) {
  
  mCLs.push_back(cls);
  mCLserr.push_back(clserr);
  mCLb.push_back(clb);
  mCLberr.push_back(clberr);
  mCLsb.push_back(clsb);
  mCLsberr.push_back(clsberr);
  mUpperLim.push_back(upperlim);
  mUpperLimerr.push_back(upperlimerr);

  return;

}

void AnalysisBase::WriteLimTree() {
  limittree->Fill();
  return;
}

void AnalysisBase::reset() {
  mCounter=0.0;
  mSigPred.clear();
  mSigPred.resize(mNumBins,0);
  mSigEff.clear();
  mSigEff.resize(mNumBins,0);
  mSBR.clear();
  mSBR.resize(mNumBins,0);
  return;
}

void AnalysisBase::resetLim() {
  mCLs.clear();
  mCLserr.clear();
  mCLb.clear();
  mCLberr.clear();
  mCLsb.clear();
  mCLsberr.clear();
  mUpperLim.clear();
  mUpperLimerr.clear();
  return;
}

void AnalysisBase::SetName(const std::string & newname) {
  mName = newname;
  return;
}

void AnalysisBase::SetFitMode(const std::string & newFitMode) {
  mFitMode = newFitMode;
  return;
}

void AnalysisBase::SetNumBins(const int & newNumBins) {
  mNumBins = newNumBins;
  return;
}

void AnalysisBase::SetExperiment(const std::string & newexp) {
  mExperiment = newexp;
  return;
}

void AnalysisBase::SetLuminosity(const double & newlumi) {
  mIntLumi = newlumi;
  return;
}

void AnalysisBase::SetDataYields(const std::vector<int> & datayields) {
  mDataYields = datayields;
  return;
}

void AnalysisBase::SetSignalPrediction(const std::vector<double> & signalpred) {
  mSigPred = signalpred;
  return;
}

void AnalysisBase::SetBackgroundPrediction(const std::vector<double> & bgpred) {
  mBGPred = bgpred;
  return;
}

/*
void AnalysisBase::AddExcConfData(const double & data) {
  mExclusionConf.push_back(data);
  return;
}

void AnalysisBase::AddUpperLimData(const double & data) {
  mUpperLim.push_back(data);
  return;
}
*/

//define as const since it doesn't change anything
std::string AnalysisBase::GetName() const {
  return mName;
}

std::string AnalysisBase::GetExperiment() const {
  return mExperiment;
}

unsigned int AnalysisBase::GetNumBins() const {
  return mNumBins;
}

double AnalysisBase::GetLuminosity() const {
  return mIntLumi;
}

std::vector<int> AnalysisBase::GetDataYields() const {
  return mDataYields;
}

std::vector<double> AnalysisBase::GetSignalPrediction() const {
  return mSigPred;
}

std::vector<double> AnalysisBase::GetBackgroundPrediction() const {
  return mBGPred;
}

std::vector<double> AnalysisBase::GetBGUncert() const {
  return mBGPredUncert;
}

std::string AnalysisBase::GetFitMode() const {
  return mFitMode;
}

bool AnalysisBase::CalculateR() const {
  return mCalculateR;
}

std::string AnalysisBase::GetCurrDir() const {
  return mCurrentDir;
}

bool AnalysisBase::operator==(const AnalysisBase & analysis) {
  return(this->GetName() == analysis.GetName());
}

AnalysisBase & AnalysisBase::operator=(const AnalysisBase & analysis) {
  //shallow copy
  if(this != &analysis) { //check for self-assignment
    
    this->SetName(analysis.GetName());
  }

  return *this; //return self-reference
}
