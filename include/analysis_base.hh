#ifndef __ANALYSISBASE__
#define __ANALYSISBASE__

#include <iostream>
#include <string>
#include <vector>

#include "fileobject_class.hh"
#include "Reader.hh"
#include "TDirectory.h"
#include "TTree.h"
#include "TMath.h"
#include "jad_lepton_class.hh"
#include "jad_jet_class.hh"
#include "jad_particle_class.hh"
#include "jad_object_collection.hh"

class AnalysisBase {

public:
  AnalysisBase();

  AnalysisBase(const std::string & name, 
	       const std::string & experiment,
	       const unsigned int & numBins);
  
  AnalysisBase(const std::string & name, 
	       const std::string & experiment, 
	       const unsigned int & numBins,
	       const double & intlumi, 
	       //const std::vector<int> & datayields,
	       const std::vector<double> & bgpred);

  AnalysisBase(const std::string & name, 
	       const std::string & experiment, 
	       const unsigned int & numBins,
	       const double & intlumi, 
	       const std::vector<double> & bgpred,
	       const std::vector<double> & bgpreduncert,
	       const std::vector<int> & datayields,
	       const std::string & fitmode,
	       const bool & calculateR);
  
  ~AnalysisBase();

  AnalysisBase(const AnalysisBase & analysisobj);

  void initDir(TDirectory * file, const std::string & ifilename);

  void SetName(const std::string & newname);
  void SetFitMode(const std::string & newname);
  void SetExperiment(const std::string & newexp);
  void SetLuminosity(const double & newlumi);
  void SetDataYields(const std::vector<int> & datayields);
  void SetSignalPrediction(const std::vector<double> & signalpred);
  void SetBackgroundPrediction(const std::vector<double> & bgpred);
  
  //void AddExcConfData(const double & data);
  //void AddUpperLimData(const double & data);

  std::string GetName() const;
  std::string GetExperiment() const;
  unsigned int GetNumBins() const;
  double GetLuminosity() const;
  std::vector<int> GetDataYields() const;
  std::vector<double> GetSignalPrediction() const;
  std::vector<double> GetBackgroundPrediction() const;
  std::vector<double> GetBGUncert() const;
  std::string GetFitMode() const;
  std::string GetCurrDir() const;
  bool CalculateR() const;

  virtual void initHistos()=0;
  void reset();
  void initTree();
  void initLimTree();
  void resetLim();
  void FillTree(const double & xsec);
  void FillLimTree(const double & cls, 
		   const double & clserr, 
		   const double & clb, 
		   const double & clberr, 
		   const double & clsb,
		   const double & clsberr,
		   const double & upperlim,
		   const double & upperlimerr);

  void WriteLimTree();

  virtual void Run(const Reader *treereader, const Reader *gentreereader, const double & weight)=0;
  AnalysisBase & operator=(const AnalysisBase & analysisobj);
  bool operator==(const AnalysisBase & analysisobj);

protected: //so we can access from derived class...
  TDirectory * andir;  
  std::vector<double> mSigPred; //signal prediction
  double mCounter; //put it here as weight depends on mIntLumi (different for each exp)

private:
  std::string mName;
  std::string mExperiment;
  std::vector<int> mDataYields;
  std::vector<double> mBGPred;
  std::vector<double> mBGPredUncert;
  double mIntLumi;
  unsigned int mNumBins;
  std::string mFitMode;
  //the tree belongs to the base class
  //since unlike the histograms, we always
  //want the same information saved for each analysis
  TTree * analysistree;
  TTree * limittree;
  std::vector<double> mSigEff; //signal efficiency
  double mSigEffTot;
  std::vector<double> mSBR; //s over sqrt(b)
  double mSBRTot;
  std::vector<double> mCLs;
  std::vector<double> mCLserr;
  std::vector<double> mCLb;
  std::vector<double> mCLberr;
  std::vector<double> mCLsb;
  std::vector<double> mCLsberr;
  std::vector<double> mUpperLim; //upper limit on R
  std::vector<double> mUpperLimerr;

  std::string mCurrentDir; //to store the dirpath for the limit code

  bool mCalculateR; //to turn on/off the calculation of the upperlimit on R in each analysis limit
};


#endif
