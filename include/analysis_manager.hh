#ifndef __ANALYSISMANAGER__
#define __ANALYSISMANAGER__

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
//#include "fileobject_class.hh"
#include "filemap_class.hh"
#include "analysis_base.hh"
#include "TreeReader.hh"
#include "TFile.h"
#include "TChain.h"

#include "TH1.h"

#include "CountingModel.h"
#include "CLsLimit.h"
#include "CRandom.h"

struct fitparams {
  double cls;
  double errs;
  double clb;
  double errb;
  double clsb;
  double errsb;
  double rtmp;
  double rtmperr;
};

class AnalysisManager {

public:
  AnalysisManager();
  ~AnalysisManager();
  AnalysisManager(const std::string & outputpath, const bool & loadgeninfo);

  void Add(AnalysisBase & analysis);
  void Run(const FileMap & fileobj);
  void RunMany(const std::vector<FileMap> & fileobjs, const std::string & foldername);
  void Limit(const double & signal_uncertainty, const bool & savestatfile, const bool & combinesearches, const bool & combinecalculater);
  void Write();

private:
  void SetupOutputFile(const std::string & name);
  bool checkDataExists(const FileMap & file);
  bool checkFitMode(const std::string & mode);
  bool inAnalysisArray(const AnalysisBase & analysis);
  void setupAnalyses(const std::string & foldername);
  void SaveStatFile(const std::string & fitmode, 
		    const unsigned int & numbins,
		    const std::vector<double> & signalyields, 
		    const double & sigerror,
		    const std::vector<double> & bgyields,
		    const std::vector<double> & bgerror,
		    const std::vector<int> & datayields,
		    const std::string & currdir);
  std::vector<fitparams> LimitCode(const std::string & fitmode,
				   const std::vector<double> & signalyields,
				   const double & signal_uncertainty,
				   const std::vector<double> & bgyields,
				   const std::vector<double> & bguncert,
				   const std::vector<int> & datayields,
				   const bool & calculateR);
  TFile * mOutputFile;    
  std::vector<AnalysisBase *> mAnalysesToRun; //could replace with vector<struct> if need more info
  std::map<std::string, std::vector<AnalysisBase *> > mAnalysesMap;
  std::string mBasePath;
  std::string mOutputPath;
  std::string mFolderName;
  bool mLoadGenInfo;
};

#endif
