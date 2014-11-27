#include "analysis_manager.hh"

AnalysisManager::AnalysisManager() {
  //resort to a default name - since we
  //recreate the file, nothing is lost
  this->SetupOutputFile("outputfile.root");
}

AnalysisManager::~AnalysisManager() {
  //last thing we do is delete the TFile pointer
  //don't ->Write() here: scope issue...
  delete mOutputFile;
  std::cout << "--ALL DONE--" << std::endl;
}

AnalysisManager::AnalysisManager(const std::string & outputpath, const bool & loadgeninfo, const bool & useEventWeight) {
  this->SetupOutputFile(outputpath+"outputfile.root");
  mOutputPath = outputpath; //set the outputpath as a private member variable
  mLoadGenInfo = loadgeninfo;
  mUseEventWeights = useEventWeight;
  if(mLoadGenInfo) {
    std::cout << "NOTE: Generator information will be read - turn this off if you are not using it (to speed things up!)" << std::endl;
  }
}

void AnalysisManager::SetupOutputFile(const std::string & name) {
  //initialise the TFile pointer
  mOutputFile = new TFile(name.c_str(), "RECREATE");
  return;
}

bool AnalysisManager::inAnalysisArray(const AnalysisBase & analysis) {

  for(std::vector<AnalysisBase *>::iterator an=mAnalysesToRun.begin(); an!=mAnalysesToRun.end(); an++) {
    //first dereference the iterator, then the pointer before comparing with the object
    if(**an == analysis) {
      return true;
    }
  }

  return false;

}

void AnalysisManager::Add(AnalysisBase & analysis) {
  //push_back each analysis into the private data member vector<analysis>
  //but first check to make sure it's not already in!
  if(!inAnalysisArray(analysis)) {
    mAnalysesToRun.push_back(&analysis);
  } else {
    std::cout << "You have already added " << analysis.GetName() << " to the analysis loop. Skipping" << std::endl;
  }

}

void AnalysisManager::Write() {
  //so we have explicit control over when
  //things are written to disk
  //i.e. keep everything in memory until the end!
  //end = after limit code + whatever plugins we dream of
  mOutputFile->Write();
}

bool AnalysisManager::checkDataExists(const FileMap & file) {

  //this method loops over the pushed_back analyses and checks if
  //the experiment required is in the list from the file object

  for(std::vector<AnalysisBase *>::iterator an=mAnalysesToRun.begin(); an!=mAnalysesToRun.end(); an++) {
    if(!file.isExpInMap((*an)->GetExperiment())) {
      return false;
    }
  }
  
  return true;

}

bool AnalysisManager::checkFitMode(const std::string & mode) {

  if(mode == "combined" || mode == "individual" || mode == "strongest") {
    return true;
  }
  return false;

}

void AnalysisManager::setupAnalyses(const std::string & foldername) {

  for(std::vector<AnalysisBase *>::iterator an=mAnalysesToRun.begin(); an!=mAnalysesToRun.end(); an++) {
    //set up analysis by initialising TDirectory and creating histograms
    (*an)->reset(); //make sure existing vectors are reset
    (*an)->initDir(mOutputFile, foldername); //make directory
    (*an)->initHistos(); //make histograms
    (*an)->initTree(); //make tree
    //at the same time, fill a map of experiments <--> analyses to loop over
    mAnalysesMap[(*an)->GetExperiment()].push_back(*an);
  }
  
}

void AnalysisManager::RunMany(const std::vector<FileMap> & files, const std::string & foldername) {

  mAnalysesMap.clear(); //reset to make sure the map is empty

  if(files.size() && mAnalysesToRun.size()) { //check there is an analysis and files

    //determine string to use for folder name
    mFolderName=foldername; //use this so that when setting limits combined over many analyses, can make unique datacard name
    if(files.size() == 1) {
      mFolderName = files[0].GetName(); //i.e. fallback on the filename when only 1 file present
    }

    this->setupAnalyses(mFolderName); //now setup map and initialise the analysis directories etc


    for(std::vector<FileMap>::const_iterator file=files.begin(); file!=files.end(); file++) {
      this->Run((*file));
    }


    //fill the trees (now that analysis is complete)
    for(std::vector<AnalysisBase *>::iterator an=mAnalysesToRun.begin(); an!=mAnalysesToRun.end(); an++) {
      //Loop over all files for each experiment and sum up the xsec. 
      //What isn't valid is the weights IF you were to simply merge the subprocesses together,
      //but the weights are calculated independently for each subprocess
      double xsec=0.0;
      for(std::vector<FileMap>::const_iterator file=files.begin(); file!=files.end(); file++) {
	xsec += (*file).GetCrossSection((*an)->GetExperiment());
      }
      (*an)->FillTree(xsec);
    }


  } else {
    std::cerr << "couldn't find any analyses and/or files to run on in the queue... quitting" << std::endl;
  }

  return;

}

void AnalysisManager::Run(const FileMap & fileobj) {
  //loop over the vector<analysis> and check that all the experiments specified
  //are available in the fileobject, otherwise complain (vector<analysis> must be filled)

  if(checkDataExists(fileobj)) {
    std::cout << "setting up workflow..." << std::endl;
      
    //loop over the map contents (of experiments <--> analyses)
    for( std::map<std::string, std::vector<AnalysisBase *> >::const_iterator ii=mAnalysesMap.begin(); ii!=mAnalysesMap.end(); ii++) {

	std::cout << "running over experiment: " << ii->first << " : " << std::endl;

	//for each experiment, open the files and run over the analyses

	std::vector<std::string> filestoadd = fileobj.GetFileList(ii->first);

	Reader *mytreereader;
	Reader *gentreereader;

	std::cout << fileobj.GetReader() << std::endl;
	TChain chainRec("Analysis");
	TChain chainGen("GEN");
	TChain chainRec1("Delphes");
	TChain chainGen1("Delphes");

	//Use (virtual or abstract?) base pointer class from here 
	if (fileobj.GetReader()==0)
	{

	    for(std::vector<std::string>::const_iterator fil=filestoadd.begin(); fil!=filestoadd.end(); fil++) {
		std::cout << "adding file to TChain: " << *fil << std::endl;
		chainRec.Add((*fil).c_str());
		chainGen.Add((*fil).c_str());
	    }
	    //If statement to choose correct reader here (based on number in filemap):
	    mytreereader= new D2Reader(&chainRec);
	    gentreereader= new D2Reader(&chainGen);

	    //Use (virtu = &d2tempgen;
	}
	else if (fileobj.GetReader()==1)
	{

	    for(std::vector<std::string>::const_iterator fil=filestoadd.begin(); fil!=filestoadd.end(); fil++) {
		std::cout << "adding file to TChain: " << *fil << std::endl;
		chainRec1.Add((*fil).c_str());
		chainGen1.Add((*fil).c_str());
	    }
	    //If statement to choose correct reader here (based on number in filemap):
	    mytreereader= new D3Reader(&chainRec1,false);
	    gentreereader= new D3Reader(&chainGen1,true);

	}

	unsigned int numevents = mytreereader->GetEntries();
	std::cout << "there are " << numevents << " events." << std::endl; 
	double sum = 0.0 ;
    if (mUseEventWeights){
      for(unsigned int counter=0; counter<numevents; counter++) {
           gentreereader->ReadEntry(counter);
           double event_weight = gentreereader->GetWeight();
           sum+=event_weight;
         }
     }
     else{
       sum = numevents;
     }
       std::cout << "The sum over the weights is: " << sum << std::endl;
	//loop over all events:
	for(unsigned int event=0; event<numevents; event++) {
	    mytreereader->ReadEntry(event);	  
	    if(mLoadGenInfo || mUseEventWeights) {
		// JM: added functionality in the manager to turn on/off generator info
		// Added to manager since you will not achieve speedup if adding to individual search (overhead in reading the event, not analysing it)
		// In principle you could NOT require geninfo for e.g. ATLAS analysis but for CMS analysis and here you'd read gen info for both
		// However, best to run twice in this instance if you want speed up
		gentreereader->ReadEntry(event);
	    }

	    double event_weight = 0.0;
        if (mUseEventWeights){
            event_weight = gentreereader->GetWeight();
        }
        else{
            event_weight = 1.0;
        }

	    //loop over all analyses that depend on this experiment
	    for(std::vector<AnalysisBase *>::const_iterator jj=(ii->second).begin(); jj != (ii->second).end(); jj++) {
		//std::cout << "running analysis: " << (*jj)->GetName() << std::endl;
//		double weight = (1.0E+015 * (*jj)->GetLuminosity())  / (numevents / fileobj.GetCrossSection(ii->first));
		double weight = (1.0E+015 * (*jj)->GetLuminosity()*event_weight)  / (sum / fileobj.GetCrossSection(ii->first));
		//Run analysis
		(*jj)->Run(mytreereader, gentreereader, weight);
	    }
	    if (event%1000 == 0) std::cout << std::setprecision(4) << event*100./numevents << "%     " << "\r" <<std::flush;
	}
	delete mytreereader;
	delete gentreereader;

    }

  } else {
      std::cerr << "couldn't find data in file... quitting" << std::endl;
  }
  return;

}

void AnalysisManager::Limit(const double & signal_uncertainty, const bool & savestatfile, const bool & combinesearches, const bool & combinecalculater) {

    std::vector<double> combined_signalyields;
    std::vector<double> combined_bgyields;
    std::vector<double> combined_bguncert;
    std::vector<int> combined_datayields;

    bool combine_searches = combinesearches;
    bool combine_calculateR = combinecalculater;
    std::string combined_fittingmode;

    //check if any analyses in the vector
    if(mAnalysesToRun.size()) {
	//check if those analyses have been run i.e. the size of the mSigPred
	for(std::vector<AnalysisBase *>::iterator an=mAnalysesToRun.begin(); an!=mAnalysesToRun.end(); an++) {
	    if((*an)->GetSignalPrediction().size() == (*an)->GetNumBins() && 
		    //signal_uncertainty.size() == (*an)->GetNumBins() &&
		    (*an)->GetBackgroundPrediction().size() == (*an)->GetNumBins() &&
		    (*an)->GetBGUncert().size() == (*an)->GetNumBins() &&
		    (*an)->GetDataYields().size() == (*an)->GetNumBins() &&
		    checkFitMode((*an)->GetFitMode())) {

		std::cout << "running limit..." << std::endl;
		(*an)->resetLim(); //make sure the vectors are empty
		(*an)->initLimTree(); //set up the limit tree

		//set up vectors
		std::vector<int> datayields = (*an)->GetDataYields();
		std::vector<double> signalyields = (*an)->GetSignalPrediction();
		std::vector<double> bgyields = (*an)->GetBackgroundPrediction();
		std::vector<double> bguncert = (*an)->GetBGUncert();

        std::vector<fitparams> fitresults; 


        if ((*an)->GetFitMode() == "strongest") {
          //Set data yields to background  
          datayields.clear();
          std::vector<double>::iterator bgIt;
          for (bgIt=bgyields.begin(); bgIt!=bgyields.end(); bgIt++){
            datayields.push_back(int(round(*bgIt)));
          }
          //Calculate limits
          std::vector<fitparams> fitresults = LimitCode("individual", 
             signalyields, signal_uncertainty, bgyields, bguncert, 
             datayields, (*an)->CalculateR());
          int maxClsIndex = 0;
          double maxCls = 0.0;
          //find the maximum index
          for (unsigned int i=0; i<fitresults.size(); i++){
            double cls = 1.0 - fitresults[i].cls;
            if (cls > maxCls){
              maxClsIndex = i;
              maxCls = cls;
            }
          }
          std::cout << "Strongest expected CLs " << maxCls << " at index : " << maxClsIndex << std::endl; 
          //set data yields back to real
          datayields.clear();
          datayields = (*an)->GetDataYields();
          //set data to the strongest expected limit
		  datayields[0] = datayields[maxClsIndex];
		  signalyields[0] = signalyields[maxClsIndex];
		  bgyields[0] = bgyields[maxClsIndex];
		  bguncert[0] = bguncert[maxClsIndex];
		  datayields.erase( datayields.begin()+1, datayields.end());
		  signalyields.erase( signalyields.begin()+1, signalyields.end());
		  bgyields.erase( bgyields.begin()+1, bgyields.end());
		  bguncert.erase( bguncert.begin()+1, bguncert.end());
          (*an)->SetFitMode("combined");
          (*an)->SetNumBins(1);
        }
		//check if data file needs to be written out (for debug or use with Higgs Limit code)
		if(savestatfile) {
		    //std::cout << "current dir is " << (*an)->GetCurrDir() << std::endl;
		    SaveStatFile((*an)->GetFitMode(), (*an)->GetNumBins(), 
                signalyields, signal_uncertainty, bgyields, bguncert, datayields, 
                (*an)->GetCurrDir());
		}

		//now run the limit code
		fitresults = LimitCode((*an)->GetFitMode(), 
                signalyields, signal_uncertainty, bgyields, bguncert, 
                datayields, (*an)->CalculateR());

		//now loop over the vector of results
		for(std::vector<fitparams>::const_iterator fp=fitresults.begin(); fp!=fitresults.end(); fp++) {
		    (*an)->FillLimTree((*fp).cls, (*fp).errs, (*fp).clb, (*fp).errb, (*fp).clsb, (*fp).errsb, (*fp).rtmp, (*fp).rtmperr);
		}

		//now write the tree
		(*an)->WriteLimTree();

		//now do some stuff for the combination
		if(combine_searches) {
		    if(combined_fittingmode.empty()) {
			combined_fittingmode = (*an)->GetFitMode();
		    } else if(combined_fittingmode != (*an)->GetFitMode()) {
			std::cout << "There is a mixture of fitting modes defined. Skipping combination..." << std::endl;
			combine_searches = false;
		    }
		    if(mAnalysesToRun.size() > 1) { //no point combining if only one analysis!!
			combined_datayields.insert(combined_datayields.end(), datayields.begin(), datayields.end());
			combined_bguncert.insert(combined_bguncert.end(), bguncert.begin(), bguncert.end());
			combined_bgyields.insert(combined_bgyields.end(), bgyields.begin(), bgyields.end());
			combined_signalyields.insert(combined_signalyields.end(), signalyields.begin(), signalyields.end());
		    } else {
			combine_searches = false;
		    }
		}


	    } else {
		std::cerr << "Either the analysis hasn't been run yet, OR" << std::endl;
		std::cerr << "There are missing/mismatched vectors, OR" << std::endl;
		std::cerr << "The fit mode was not recognised." << std::endl;
		std::cerr << "Please check! Quitting..." << std::endl;
	    }
	}

	if(combine_searches && combined_datayields.size() && combined_bguncert.size() && combined_bgyields.size() && combined_signalyields.size()) {

	    if(savestatfile) {
		SaveStatFile(combined_fittingmode, combined_signalyields.size(), combined_signalyields, signal_uncertainty, combined_bgyields, combined_bguncert, combined_datayields, mFolderName+"_combined_datacard");
	    }


	    std::vector<fitparams> fitresultscombined = LimitCode(combined_fittingmode, combined_signalyields, signal_uncertainty, combined_bgyields, combined_bguncert, combined_datayields, combine_calculateR);

	    //this bit of code is a bit of an afterthought since it was never envisaged that the manager would combine search limits...

	    std::vector<double> mCLs;
	    std::vector<double> mCLserr;
	    std::vector<double> mCLb;
	    std::vector<double> mCLberr;
	    std::vector<double> mCLsb;
	    std::vector<double> mCLsberr;
	    std::vector<double> mUpperLim; //upper limit on R
	    std::vector<double> mUpperLimerr;

	    mOutputFile->cd();
	    TDirectory * combdir = mOutputFile->mkdir(mFolderName.c_str(), mFolderName.c_str());
	    combdir->cd();
	    TTree * limittree = new TTree("LimitTreec","LimitTreec");
	    limittree->Branch("cls", &mCLs);
	    limittree->Branch("clserr", &mCLserr);
	    limittree->Branch("clb", &mCLb);
	    limittree->Branch("clberr", &mCLberr);
	    limittree->Branch("clsb", &mCLsb);
	    limittree->Branch("clsberr", &mCLsberr);
	    limittree->Branch("upperlimitR", &mUpperLim);
	    limittree->Branch("upperlimitRerr", &mUpperLimerr);

	    for(std::vector<fitparams>::const_iterator fp=fitresultscombined.begin(); fp!=fitresultscombined.end(); fp++) {
		mCLs.push_back((*fp).cls);
		mCLserr.push_back((*fp).errs);
		mCLb.push_back((*fp).clb);
		mCLberr.push_back((*fp).errb);
		mCLsb.push_back((*fp).clsb);
		mCLsberr.push_back((*fp).errsb);
		mUpperLim.push_back((*fp).rtmp);
		mUpperLimerr.push_back((*fp).rtmperr);
	    }

	    limittree->Fill();


	} else {
	    std::cout << "Combination not run" << std::endl;
	}

    } else {
	std::cerr << "couldn't find any analyses in the queue... quitting" << std::endl;
    }

    /*
       CRandom *rdm = new CRandom(seed);  //initialise a random generator

    // set up counting model
    lands::CountingModel *cms=new lands::CountingModel();
    cms->SetRdm(rdm);
    cms->SetModelName("alphat");

    std::vector<double> tmpsigbkgs;
    tmpsigbkgs.push_back(1.87372);
    tmpsigbkgs.push_back(787);	

    std::vector<string> tmpprocnn;
    tmpprocnn.push_back("sig");
    tmpprocnn.push_back("bg");
    cms->AddChannel("b1", tmpsigbkgs, 1);
    cms->SetProcessNames(0, tmpprocnn);
    cms->AddObservedData(0, 782);

    tmpsigbkgs.at(0) = 1.73844;
    tmpsigbkgs.at(1) = 310;
    cms->AddChannel("b2", tmpsigbkgs, 1);
    cms->SetProcessNames(1, tmpprocnn);
    cms->AddObservedData(1, 321);

    tmpsigbkgs.at(0) = 2.62457;
    tmpsigbkgs.at(1) = 202;
    cms->AddChannel("b3", tmpsigbkgs, 1);
    cms->SetProcessNames(2, tmpprocnn);
    cms->AddObservedData(2, 196);

    tmpsigbkgs.at(0) = 1.37316;
    tmpsigbkgs.at(1) = 60.4;
    cms->AddChannel("b4", tmpsigbkgs, 1);
    cms->SetProcessNames(3, tmpprocnn);
    cms->AddObservedData(3, 62);

    tmpsigbkgs.at(0) = 0.62232;
    tmpsigbkgs.at(1) = 20.3;
    cms->AddChannel("b5", tmpsigbkgs, 1);
    cms->SetProcessNames(4, tmpprocnn);
    cms->AddObservedData(4, 21);

    tmpsigbkgs.at(0) = 0.304396;
    tmpsigbkgs.at(1) = 7.7;
    cms->AddChannel("b6", tmpsigbkgs, 1);
    cms->SetProcessNames(5, tmpprocnn);
    cms->AddObservedData(5, 6);

    tmpsigbkgs.at(0) = 0.128523;
    tmpsigbkgs.at(1) = 3.2;
    cms->AddChannel("b7", tmpsigbkgs, 1);
    cms->SetProcessNames(6, tmpprocnn);
    cms->AddObservedData(6, 3);

    tmpsigbkgs.at(0) = 0.229988;
    tmpsigbkgs.at(1) = 2.8;
    cms->AddChannel("b8", tmpsigbkgs, 1);
    cms->SetProcessNames(7, tmpprocnn);
    cms->AddObservedData(7, 1);

    //uncertainty on the signal
    cms->AddUncertainty(0, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);
    cms->AddUncertainty(1, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);
    cms->AddUncertainty(2, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);
    cms->AddUncertainty(3, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);
    cms->AddUncertainty(4, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);
    cms->AddUncertainty(5, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);
    cms->AddUncertainty(6, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);
    cms->AddUncertainty(7, 0, 0.15, 0.15, 1, "signal" );
    cms->TagUncertaintyFloatInFit("signal", 1);

    //uncertainty on the backgrounds
    cms->AddUncertainty(0, 1, 0.05, 0.05, 1, "bgb1" );
    cms->TagUncertaintyFloatInFit("bgb1", 1);
    cms->AddUncertainty(1, 1, 0.05, 0.05, 1, "bgb2" );
    cms->TagUncertaintyFloatInFit("bgb2", 1);
    cms->AddUncertainty(2, 1, 0.05, 0.05, 1, "bgb3" );
    cms->TagUncertaintyFloatInFit("bgb3", 1);
    cms->AddUncertainty(3, 1, 0.07, 0.07, 1, "bgb4" );
    cms->TagUncertaintyFloatInFit("bgb4", 1);
    cms->AddUncertainty(4, 1, 0.09, 0.09, 1, "bgb5" );
    cms->TagUncertaintyFloatInFit("bgb5", 1);
    cms->AddUncertainty(5, 1, 0.10, 0.10, 1, "bgb6" );
    cms->TagUncertaintyFloatInFit("bgb6", 1);
    cms->AddUncertainty(6, 1, 0.13, 0.13, 1, "bgb7" );
    cms->TagUncertaintyFloatInFit("bgb7", 1);
    cms->AddUncertainty(7, 1, 0.15, 0.15, 1, "bgb8" );
    cms->TagUncertaintyFloatInFit("bgb8", 1);	



    cms->SetUseSystematicErrors(true);
    cms->ForceSymmetryError(false);
    cms->RemoveChannelsWithExpectedSignal0orBkg0(0);
    cms->MultiSigProcShareSamePDF(false);
    cms->SetMoveUpShapeUncertainties(1);
    cms->SetMass(-1);

    cms->SetTossToyConvention(0);
    cms->SetUseBestEstimateToCalcQ(1);

    // initialize the calculator
    lands::CLsBase frequentist;
    frequentist.SetDebug(debug);
    frequentist.SetRdm(rdm);
    frequentist.SetTestStatistics(testStatistics);
    lands::cms_global= cms;

    double tmp;
    lands::vdata_global=cms->Get_v_data();

    frequentist.SetModel(cms);

    frequentist.prepareLogNoverB();

    frequentist.BuildM2lnQ_data();
    frequentist.BuildM2lnQ_sb(nexps, false, false);
    frequentist.BuildM2lnQ_b(nexps, false, false);
    double errs, errb, errsb;
    double cls = frequentist.CLs(errs);
    double clsb = frequentist.CLsb(errsb);
    double clb = frequentist.CLb(errb);
    std::cout<<"------------------------------------------------------------"<<std::endl;
    std::cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<std::endl;
    std::cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<std::endl;
    std::cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<std::endl;
    std::cout<<"------------------------------------------------------------"<<std::endl;


    lands::CLsLimit clsr95;
    clsr95.SetDebug(debug);
    clsr95.SetRule(rule);
    double rtmp;
    clsr95.SetAlpha(0.05);
    cms->SetSignalScaleFactor(1.);


    rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist,nexps);


    std::cout<<"------------------------------------------------------------"<<std::endl;
    std::cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<" +/- "<<clsr95.LimitErr()<<std::endl;
    std::cout<<"------------------------------------------------------------"<<std::endl;


    delete cms;
    delete rdm;
    */
}

std::vector<fitparams> AnalysisManager::LimitCode(const std::string & fitmode, const std::vector<double> & signalyields, const double & signal_uncertainty, const std::vector<double> & bgyields, const std::vector<double> & bguncert, const std::vector<int> & datayields, const bool & calculateR) {

    //initialise iterators
    std::vector<int>::const_iterator dy = datayields.begin();
    std::vector<double>::const_iterator sy = signalyields.begin();
    std::vector<double>::const_iterator by = bgyields.begin();
    std::vector<double>::const_iterator byunc = bguncert.begin();

    //there are two fitting modes at this point - combined or on their own

    //define some constants for the fit:
    int nexps=10000;
    int seed =1234;
    int testStatistics = 1;
    int rule = 1; // default is CLs

    //setup results
    std::vector<fitparams> results;

    if(fitmode == "combined") {

	CRandom * rdm = new CRandom(seed);  //initialise a random generator

	//set up counting model
	lands::CountingModel * cms = new lands::CountingModel();
	cms->SetRdm(rdm);
	cms->SetModelName("model");

	std::vector<double> tmpsigbkgs(2,0.0);
	std::vector<string> tmpprocnn;
	tmpprocnn.push_back("sig");
	tmpprocnn.push_back("bg");

	unsigned int tmpindex = 0;

	for(; dy != datayields.end() && 
		sy != signalyields.end() &&
		by != bgyields.end() &&
		byunc != bguncert.end()
		; dy++, sy++, by++, byunc++) {

	    stringstream bgident;
	    bgident << "b" << tmpindex+1;
	    stringstream bgbident;
	    bgbident << "bgb" << tmpindex+1;

	    tmpsigbkgs.at(0) = (*sy);
	    tmpsigbkgs.at(1) = (*by);

	    //add channels
	    cms->AddChannel(bgident.str(), tmpsigbkgs, 1);
	    cms->SetProcessNames(tmpindex, tmpprocnn);
	    cms->AddObservedData(tmpindex, (*dy));

	    //add uncertainty on signal - correlated between bins
	    cms->AddUncertainty(tmpindex, 0, signal_uncertainty, signal_uncertainty, 1, "signal");
	    cms->TagUncertaintyFloatInFit("signal", 1);

	    //add uncertainty on background
	    cms->AddUncertainty(tmpindex, 1, (*byunc), (*byunc), 1, bgbident.str() );
	    cms->TagUncertaintyFloatInFit(bgbident.str(), 1);

	    tmpindex++;
	}

	cms->SetUseSystematicErrors(true);
	cms->ForceSymmetryError(false);
	cms->RemoveChannelsWithExpectedSignal0orBkg0(0);
	cms->MultiSigProcShareSamePDF(false);
	cms->SetMoveUpShapeUncertainties(1);
	cms->SetMass(-1);

	cms->SetTossToyConvention(0);
	cms->SetUseBestEstimateToCalcQ(1);

	// initialize the calculator
	lands::CLsBase frequentist;
	frequentist.SetDebug(0);
	frequentist.SetRdm(rdm);
	frequentist.SetTestStatistics(testStatistics);
	lands::cms_global= cms;

	double tmp;
	lands::vdata_global=cms->Get_v_data();

	frequentist.SetModel(cms);

	frequentist.prepareLogNoverB();

	frequentist.BuildM2lnQ_data();
	frequentist.BuildM2lnQ_sb(nexps, false, false);
	frequentist.BuildM2lnQ_b(nexps, false, false);
	double errs, errb, errsb;
	double cls = frequentist.CLs(errs);
	double clsb = frequentist.CLsb(errsb);
	double clb = frequentist.CLb(errb);
	//std::cout<<"------------------------------------------------------------"<<std::endl;
	//std::cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<std::endl;
	//std::cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<std::endl;
	//std::cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<std::endl;
	//std::cout<<"------------------------------------------------------------"<<std::endl;


	lands::CLsLimit clsr95;
	clsr95.SetDebug(0);
	clsr95.SetRule(rule);
	double rtmp;
	clsr95.SetAlpha(0.05);
	cms->SetSignalScaleFactor(1.);

	//JM: 24012013: commented out upperlimit calculation - now it's a function
	//the range for which you try to match CLs with upplimR, and the number of points to interpolate
	if(calculateR) {
	    rtmp = clsr95.LimitOnSignalScaleFactor(cms, 0.1, 3.0, &frequentist, nexps, 50);
	}


	//std::cout<<"------------------------------------------------------------"<<std::endl;
	//std::cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<" +/- "<<clsr95.LimitErr()<<std::endl;
	//std::cout<<"------------------------------------------------------------"<<std::endl;


	delete cms;
	delete rdm;

	fitparams fitresult;
	fitresult.cls = cls;
	fitresult.errs = errs;
	fitresult.clb = clb;
	fitresult.errb = errb;
	fitresult.clsb = clsb;
	fitresult.errsb = errsb;
	//JM: update for calculateR flag
	if(calculateR) {
	    fitresult.rtmp = rtmp;
	    fitresult.rtmperr = clsr95.LimitErr();
	} else {
	    fitresult.rtmp = 0.0;
	    fitresult.rtmperr = 0.0;
	}
	results.push_back(fitresult);

    } else if(fitmode == "individual") {


	for(; dy != datayields.end() && 
		sy != signalyields.end() &&
		by != bgyields.end() &&
		byunc != bguncert.end()
		; dy++, sy++, by++, byunc++) {

	    if((*sy) > 0) { //i.e. only look at bins with signal > 0

		CRandom * rdm = new CRandom(seed);  //initialise a random generator

		//set up counting model
		lands::CountingModel * cms = new lands::CountingModel();
		cms->SetRdm(rdm);
		cms->SetModelName("model");

		std::vector<double> tmpsigbkgs(2,0.0);
		std::vector<string> tmpprocnn;
		tmpprocnn.push_back("sig");
		tmpprocnn.push_back("bg");

		unsigned int tmpindex = 0;

		stringstream bgident;
		bgident << "b" << tmpindex+1;
		stringstream bgbident;
		bgbident << "bgb" << tmpindex+1;

		tmpsigbkgs.at(0) = (*sy);
		tmpsigbkgs.at(1) = (*by);

		//add channels
		cms->AddChannel(bgident.str(), tmpsigbkgs, 1);
		cms->SetProcessNames(tmpindex, tmpprocnn);
		cms->AddObservedData(tmpindex, (*dy));

		//add uncertainty on signal - correlated between bins
		cms->AddUncertainty(tmpindex, 0, signal_uncertainty, signal_uncertainty, 1, "signal");
		cms->TagUncertaintyFloatInFit("signal", 1);

		//add uncertainty on background
		cms->AddUncertainty(tmpindex, 1, (*byunc), (*byunc), 1, bgbident.str() );
		cms->TagUncertaintyFloatInFit(bgbident.str(), 1);


		cms->SetUseSystematicErrors(true);
		cms->ForceSymmetryError(false);
		cms->RemoveChannelsWithExpectedSignal0orBkg0(0);
		cms->MultiSigProcShareSamePDF(false);
		cms->SetMoveUpShapeUncertainties(1);
		cms->SetMass(-1);

		cms->SetTossToyConvention(0);
		cms->SetUseBestEstimateToCalcQ(1);

		// initialize the calculator
		lands::CLsBase frequentist;
		frequentist.SetDebug(0);
		frequentist.SetRdm(rdm);
		frequentist.SetTestStatistics(testStatistics);
		lands::cms_global= cms;

		double tmp;
		lands::vdata_global=cms->Get_v_data();

		frequentist.SetModel(cms);

		frequentist.prepareLogNoverB();

		frequentist.BuildM2lnQ_data();
		frequentist.BuildM2lnQ_sb(nexps, false, false);
		frequentist.BuildM2lnQ_b(nexps, false, false);
		double errs, errb, errsb;
		double cls = frequentist.CLs(errs);
		double clsb = frequentist.CLsb(errsb);
		double clb = frequentist.CLb(errb);
		//std::cout<<"------------------------------------------------------------"<<std::endl;
		//std::cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<std::endl;
		//std::cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<std::endl;
		//std::cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<std::endl;
		//std::cout<<"------------------------------------------------------------"<<std::endl;


		lands::CLsLimit clsr95;
		clsr95.SetDebug(0);
		clsr95.SetRule(rule);
		double rtmp;
		clsr95.SetAlpha(0.05);
		cms->SetSignalScaleFactor(1.);

		//JM: 24012013: commented out upperlimit calculation -- same as above, make configurable
		if(calculateR) {
		    rtmp = clsr95.LimitOnSignalScaleFactor(cms, 0.1, 3.0, &frequentist, nexps, 50);
		}


		//std::cout<<"------------------------------------------------------------"<<std::endl;
		//std::cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<" +/- "<<clsr95.LimitErr()<<std::endl;
		//std::cout<<"------------------------------------------------------------"<<std::endl;


		delete cms;
		delete rdm;

		fitparams fitresult;
		fitresult.cls = cls;
		fitresult.errs = errs;
		fitresult.clb = clb;
		fitresult.errb = errb;
		fitresult.clsb = clsb;
		fitresult.errsb = errsb;
		if(calculateR) {
		    fitresult.rtmp = rtmp;
		    fitresult.rtmperr = clsr95.LimitErr();
		} else {
		    fitresult.rtmp = 0.0;
		    fitresult.rtmperr = 0.0;
		}
		results.push_back(fitresult);

	    } else {
		//the signal expectation was 0 so ignore this bin
		fitparams fitresult;
		fitresult.cls = 0.0;
		fitresult.errs = 0.0;
		fitresult.clb = 0.0;
		fitresult.errb = 0.0;
		fitresult.clsb = 0.0;
		fitresult.errsb = 0.0;
		fitresult.rtmp = 0.0;
		fitresult.rtmperr = 0.0;
		results.push_back(fitresult);
	    }
	}
    }

    return results;

}

void AnalysisManager::SaveStatFile(const std::string & fitmode, 
	const unsigned int & numbins,
	const std::vector<double> & signalyields, 
	const double & sigerror,
	const std::vector<double> & bgyields,
	const std::vector<double> & bguncert,
	const std::vector<int> & datayields,
	const std::string & currdir) {

    std::string filename = mOutputPath + currdir + ".txt"; //mBasePath

    //initialise relevant iterators
    std::vector<int>::const_iterator dy = datayields.begin();
    std::vector<double>::const_iterator sy = signalyields.begin();
    std::vector<double>::const_iterator by = bgyields.begin();
    std::vector<double>::const_iterator byunc = bguncert.begin();


    std::ofstream myfile;
    myfile.open(filename.c_str());

    std::cout << "writing out to file " << filename << std::endl;

    if(fitmode == "combined") {

	myfile << "# Counting experiment with multiple channels" << std::endl;
	myfile << "imax " << numbins << "  number of channels" << std::endl;
	myfile << "jmax 1  number of backgrounds" << std::endl;
	myfile << "kmax " << numbins+1 << "  number of nuisance parameters (sources of systematical uncertainties)" << std::endl;
	myfile << "------------" << std::endl;
	myfile << "# n bins" << std::endl;
//	myfile << "bin            ";
//
//	for(unsigned int i=0; i<numbins; i++) {
//	    myfile << "b" << i+1;
//	    if(i==(numbins - 1)) { myfile << std::endl; }
//	    else { myfile << "  "; }
//	}

	myfile << "Observation    ";
	for(unsigned int i=0; i<numbins; i++) {
	    myfile << datayields.at(i);
	    if(i==(numbins - 1)) { myfile << std::endl; }
	    else { myfile << "  "; }
	}

	myfile << "------------" << std::endl;

	myfile << "# now we list the expected events for signal and all backgrounds in that bin" << std::endl;
	myfile << "# the second 'process' line must have a positive number for backgrounds, and 0 for signal" << std::endl;
	myfile << "# now we list the independent sources of uncertainties, and give their effect (syst. error)" << std::endl;
	myfile << "# on each process and bin" << std::endl;

	myfile << "bin            ";
	for(unsigned int i=0; i<numbins; i++) { 
	    myfile <<  i+1 << "    "  << i+1;
	    if(i==(numbins - 1)) { myfile << std::endl; }
	    else { myfile << "   "; }
	}

//	myfile << "process        ";
//	for(unsigned int i=0; i<numbins; i++) {
//	    myfile << "sig    bg";
//	    if(i==(numbins - 1)) { myfile << std::endl; }
//	    else { myfile << "    "; } 
//	}

	myfile << "process         ";
	for(unsigned int i=0; i<numbins; i++) {
	    myfile << "0    1";
	    if(i==(numbins - 1)) { myfile << std::endl; }
	    else { myfile << "    "; }      
	}

	myfile << "rate           ";
	for(unsigned int i=0; i<numbins; i++) {
	    myfile << signalyields.at(i) <<"   " << bgyields.at(i);
	    if(i==(numbins - 1)) { myfile << std::endl; }
	    else { myfile << "  "; }
	}

	myfile << "------------" << std::endl;
	myfile << "1   lnN   ";
	for(unsigned int i=0; i<numbins; i++) {
	    myfile << sigerror + 1.00 << "   " << " -";
	    if(i==(numbins - 1)) { myfile << "   15% uncertainty on signal" << std::endl; }
	    else { myfile << "   "; }
	}

	unsigned int countme=0; //to place the required bg uncertainty in the right order since number of entries depends on mask
	for(unsigned int i=0; i<numbins; i++) {
	    myfile <<  i+2 << "    lnN     ";
	    for(unsigned int j=0; j<numbins; j++) {
		if(j == countme) { myfile << "-    " << bguncert.at(i) + 1.00 << "    ";}
		else {myfile << "-     -    ";}
	    }
	    myfile << "x% on background in bin" << std::endl;
	    countme++; 
	}



    } else if(fitmode == "individual") {

	for(; dy != datayields.end() && 
		sy != signalyields.end() &&
		by != bgyields.end() &&
		byunc != bguncert.end()
		; dy++, sy++, by++, byunc++) {
	    myfile << (*sy) << "\t" << sigerror << "\t" << (*by) << "\t" << (*byunc) << "\t" << (*dy) <<std::endl;
	}
    }


    myfile.close();

    return;
}
