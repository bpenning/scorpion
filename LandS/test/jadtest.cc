#include <iostream>
#include "CLsLimit.h"
#include "CRandom.h"
//#include "PlotUtilities.h"
//#include "FloridaStyle.C"
//#include <TError.h>
//#include "TSystem.h"
#include "CountingModel.h"
//#include "BayesianBase.h"
//#include "LimitBands.h"
//#include "SignificanceBands.h"
//#include "UtilsROOT.h"
//#include "TGraph.h"
using std::cout;
using std::endl;
using namespace lands;
//void processParameters(int argc, const char* argv[]);
int main(int argc, const char* argv[]){
  
double s=0;  double b=0; double s_err = 0; double b_err = 0; int d=0;
int pdftypeEs = 1; int pdftypeEb = 1; int EsEb_correlated = 0; 
  

// current for FeldmanCousins
bool bAdaptiveSampling = true;
bool bQuickEstimateInitialLimit = true;
double toysFactor_ = 1.;
double rAbsAccuracy_ = 0.01, rRelAccuracy_ = 0.01;
//double rAbsAccuracy_ = 0.1;
//double rRelAccuracy_ = 0.05;
bool lowerLimit_ = false;
double intialRmin = 1., intialRmax = 21;// only when bQuickEstimateInitialLimit==false

const char* fileName;
//int main(int argc, const char* argv[]){
//processParameters(argc, argv);

int calcExpectedMeanLimit=0; 
int nexps=10000;
int seed =1234;
int debug=0;
int testStatistics = 1;
int rule = 1; // default is CLs
int asimov = -1; // -1 for using the input dataset whatever you provide,  0 for using asimov bkg only dataset, 1 for asimov sig+bkg
/*
 if(argc==2) { fileName = argv[1]; }
 else {
   cout<<"please use following format:"<<endl;
   cout<<"./CLs_dataCard.exe inputFileName"<<endl;
   cout<<endl;
   cout<<" or "<<endl;
   cout<<"./CLs_dataCard.exe inputFileName calcExpectedMeanLimit=0 ntoys=100000 seed=1234 debug=0 testStatistics=1 rule=1"<<endl;
   cout<<endl<<" detail of the parameters: "<<endl;
   cout<<" inputFileName:          The input data card designed by Andrey Korytov "<<endl;
   cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
   cout<<" ntoys: 			number of toy experiments to build -2lnQ "<<endl;
   cout<<" seed:			seed for random number generator "<<endl;
   cout<<" debug: 			debug level "<<endl;
   cout<<" testStatistics:		1 for Q_LEP, 2 for Q_TEV, 4 for only profile mu"<<endl;
   cout<<"                         3 for Q_ATLAS, 31 for Q_ATLAS but allowing mu_hat>mu"<<endl;
   cout<<" rule:                   1 for CLs,  2 for CLsb, 3 for FeldmanCousins"<<endl;
   cout<<" asimov:                 0 for bkg only, 1 for sig+bkg, others for the input whatever you provide "<<endl;
   exit(0);
 }
*/
 /*
	// print out  parameters  you just configured:
	cout<<"\n\t Your inputs are the following: "<<endl;	
	cout<<"**********************************************"<<endl;
	cout<<" Input data card is:   "<<fileName<<endl; 
	cout<<(calcExpectedMeanLimit==0?" NOT":"")<<" do calculation of mean limit"<<endl;
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;

	cout<<" testStatistics = "<<testStatistics<<", is ";
	if(testStatistics==2) cout<<" Tevatron type";
	else if(testStatistics==3) cout<<" ATLAS type,  Lamda(mu)";
	else if(testStatistics==31) cout<<" ATLAS type but allowing mu_hat>mu";
	else if(testStatistics==4) cout<<" only profile mu";
	else cout<<" LEP type";
	cout<<endl;

	cout<<" rule = "<< rule<<", is ";
	if(rule==2) cout<<" CLsb";
	else if(rule==1) cout<<" CLs";
	else if(rule==3) cout<<" FeldmanCousins";
	cout<<endl;

	if(asimov==0) cout<<" Using Asimov bkg only dataset"<<endl;
	if(asimov==1) cout<<" Using Asimov signal + bkg  dataset"<<endl;

	cout<<"**********************************************"<<endl;

*/

	//FloridaStyle();
	//TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);

	CRandom *rdm = new CRandom(seed);  //initilize a random generator

	// set up couting model
	CountingModel *cms=new CountingModel();
	cms->SetRdm(rdm);
	//ConfigureModel(cms, -1, fileName, debug); //this function basically parses the datafile
	//look in UtilsRoot.cc:
	//cms->AddChannel
	//cms->SetProcessName
	//cms->AddObservedData
	//cms->AddUncertainty
	//cms->TagUncertaintyFloatinfit

	//replace ConfigureModel with:
	//cms->SetModelName(fileName);
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
	//cms->UseAsimovData(asimov);
	cms->RemoveChannelsWithExpectedSignal0orBkg0(0);
	cms->MultiSigProcShareSamePDF(false);
	cms->SetMoveUpShapeUncertainties(1);
	cms->SetMass(-1);
	//cms->Print();

	cms->SetTossToyConvention(0);
	cms->SetUseBestEstimateToCalcQ(1);

	//cms->SetAllowNegativeSignalStrength(false);

	// initialize the calculator
	CLsBase frequentist;
	frequentist.SetDebug(debug);
	frequentist.SetRdm(rdm);
	frequentist.SetTestStatistics(testStatistics);
	cms_global= cms;
	//vdata_global=cms->Get_v_data();

	double tmp;
	vdata_global=cms->Get_v_data();

	frequentist.SetModel(cms);

	frequentist.prepareLogNoverB();

	//if(rule==1 || rule==2)frequentist.BuildM2lnQ(cms,nexps);
	//else frequentist.BuildM2lnQ(cms, 100);
	frequentist.BuildM2lnQ_data();
	frequentist.BuildM2lnQ_sb(nexps, false, false);
	frequentist.BuildM2lnQ_b(nexps, false, false);
	double errs, errb, errsb;
	double cls = frequentist.CLs(errs);
	double clsb = frequentist.CLsb(errsb);
	double clb = frequentist.CLb(errb);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<endl;
	cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<endl;
	cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<endl;
	cout<<"------------------------------------------------------------"<<endl;


	CLsLimit clsr95;
	clsr95.SetDebug(debug);
	clsr95.SetRule(rule);
	double rtmp;
	clsr95.SetAlpha(0.05);
	cms->SetSignalScaleFactor(1.);


	rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist,nexps);


	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<" +/- "<<clsr95.LimitErr()<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	return 1;
}
/*
void processParameters(int argc, const char* argv[]){
	int npar=1;
	if(argc>=2){
		fileName = argv[1];
		npar++;
		if(argc>=npar+1){
			calcExpectedMeanLimit=atoi(argv[npar]);
			npar++;
			if(argc>=npar+1){
				nexps=atoi( argv[npar] );			
				npar++;
				if(argc>=npar+1){
					seed=atoi( argv[npar] );			
					npar++;
					if(argc>=npar+1){
						debug=atoi( argv[npar] );			
						npar++;
						if(argc>=npar+1){
							testStatistics=atoi( argv[npar] );			
							npar++;
							if(argc>=npar+1){
								rule=atoi( argv[npar] );			
								npar++;
								if(argc>=npar+1){
									asimov=atoi( argv[npar] );			
								}
							}
						}
					}
				}
			}
		}
	}
	if(argc<2) {
		cout<<"please use following format:"<<endl;
		cout<<"./CLs_dataCard.exe inputFileName"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./CLs_dataCard.exe inputFileName calcExpectedMeanLimit=0 ntoys=100000 seed=1234 debug=0 testStatistics=1 rule=1"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" inputFileName:          The input data card designed by Andrey Korytov "<<endl;
		cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
		cout<<" ntoys: 			number of toy experiments to build -2lnQ "<<endl;
		cout<<" seed:			seed for random number generator "<<endl;
		cout<<" debug: 			debug level "<<endl;
		cout<<" testStatistics:		1 for Q_LEP, 2 for Q_TEV, 4 for only profile mu"<<endl;
		cout<<"                         3 for Q_ATLAS, 31 for Q_ATLAS but allowing mu_hat>mu"<<endl;
		cout<<" rule:                   1 for CLs,  2 for CLsb, 3 for FeldmanCousins"<<endl;
		cout<<" asimov:                 0 for bkg only, 1 for sig+bkg, others for the input whatever you provide "<<endl;
		exit(0);
	}
	// print out  parameters  you just configured:
	cout<<"\n\t Your inputs are the following: "<<endl;	
	cout<<"**********************************************"<<endl;
	cout<<" Input data card is:   "<<fileName<<endl; 
	cout<<(calcExpectedMeanLimit==0?" NOT":"")<<" do calculation of mean limit"<<endl;
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;

	cout<<" testStatistics = "<<testStatistics<<", is ";
	if(testStatistics==2) cout<<" Tevatron type";
	else if(testStatistics==3) cout<<" ATLAS type,  Lamda(mu)";
	else if(testStatistics==31) cout<<" ATLAS type but allowing mu_hat>mu";
	else if(testStatistics==4) cout<<" only profile mu";
	else cout<<" LEP type";
	cout<<endl;

	cout<<" rule = "<< rule<<", is ";
	if(rule==2) cout<<" CLsb";
	else if(rule==1) cout<<" CLs";
	else if(rule==3) cout<<" FeldmanCousins";
	cout<<endl;

	if(asimov==0) cout<<" Using Asimov bkg only dataset"<<endl;
	if(asimov==1) cout<<" Using Asimov signal + bkg  dataset"<<endl;

	cout<<"**********************************************"<<endl;
}
*/
