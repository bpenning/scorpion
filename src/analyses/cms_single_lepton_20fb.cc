#include "cms_single_lepton_20fb_functions.hh"


///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// mt2w ////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
double calculateMT2w(std::vector<LorentzVector>& jets, std::vector<bool>& btag, LorentzVector& lep, float met, float metphi){

    // I am asumming that jets is sorted by Pt
    assert ( jets.size() == btag.size() );
    // require at least 2 jets
    if ( jets.size()<2 ) return 99999.; 

    // First we count the number of b-tagged jets, and separate those non b-tagged
    std::vector<int> bjets;
    std::vector<int> non_bjets;
    for( unsigned int i = 0 ; i < jets.size() ; i++ ){
        if( btag.at(i) ) {
            bjets.push_back(i);
        } else {
            non_bjets.push_back(i);
        }
    }	

    int n_btag = (int) bjets.size();
    //	cout << "n_btag = " << n_btag << endl;

    // We do different things depending on the number of b-tagged jets
    // arXiv:1203.4813 recipe

    int nMax=-1;
    if(jets.size()<=3) nMax=non_bjets.size();
    else nMax=3;

    if (n_btag == 0){                  //  0 b-tags
        // If no b-jets select the minimum of the mt2w from all combinations with 
        // the three leading jets
        float min_mt2w = 9999;

        for (int i=0; i<nMax; i++)
            for (int j=0; j<nMax; j++){
                if (i == j) continue;
                float c_mt2w = mt2wWrapper(lep, 
                        jets[non_bjets[i]],
                        jets[non_bjets[j]], met, metphi);
                if (c_mt2w < min_mt2w)
                    min_mt2w = c_mt2w;
            }
        return min_mt2w;

    } else if (n_btag == 1 ){          //  1 b-tags
        // if only one b-jet choose the three non-b leading jets and choose the smaller
        float min_mt2w = 9999;

        for (int i=0; i<nMax; i++){
            float c_mt2w = mt2wWrapper(lep, jets[bjets[0]], jets[non_bjets[i]], met, metphi);
            if (c_mt2w < min_mt2w)
                min_mt2w = c_mt2w;
        }
        for (int i=0; i<nMax; i++){
            float c_mt2w = mt2wWrapper(lep, jets[non_bjets[i]], jets[bjets[0]], met, metphi);
            if (c_mt2w < min_mt2w)
                min_mt2w = c_mt2w;
        }
        return min_mt2w;

    } else if (n_btag >= 2) {          // >=2 b-tags
        // if 3 or more b-jets the paper says ignore b-tag and do like 0-bjets 
        // but we are going to make the combinations with the b-jets
        float min_mt2w = 9999;
        for (int i=0; i<n_btag; i++)
            for (int j=0; j<n_btag; j++){
                if (i == j) continue;
                float c_mt2w = mt2wWrapper(lep, 
                        jets[bjets[i]],
                        jets[bjets[j]], met, metphi);
                if (c_mt2w < min_mt2w)
                    min_mt2w = c_mt2w;
            }
        return min_mt2w;
    }

    return -1.;
}

// This funcion is a wrapper for mt2w_bisect etc that takes LorentzVectors instead of doubles
double mt2wWrapper(LorentzVector& lep, LorentzVector& jet_o, LorentzVector& jet_b, float met, float metphi){

    // same for all MT2x variables
    float metx = met * cos( metphi );
    float mety = met * sin( metphi );

    double pl[4];     // Visible lepton
    double pb1[4];    // bottom on the same side as the visible lepton
    double pb2[4];    // other bottom, paired with the invisible W
    double pmiss[3];  // <unused>, pmx, pmy   missing pT

    pl[0]= lep.E(); pl[1]= lep.Px(); pl[2]= lep.Py(); pl[3]= lep.Pz();
    pb1[1] = jet_o.Px();  pb1[2] = jet_o.Py();   pb1[3] = jet_o.Pz();
    pb2[1] = jet_b.Px();  pb2[2] = jet_b.Py();   pb2[3] = jet_b.Pz();
    pmiss[0] = 0.; pmiss[1] = metx; pmiss[2] = mety;

    pb1[0] = jet_o.E();
    pb2[0] = jet_b.E();

    mt2w_bisect::mt2w mt2w_event;
    mt2w_event.set_momenta(pl, pb1, pb2, pmiss);

    return mt2w_event.get_mt2w();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// chi2 ////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double fc2 (double c1, double m12, double m22, double m02, bool verbose)
{
  if (verbose) {
    printf("c1: %4.2f\n", c1);
    printf("m12: %4.2f\n", m12);
    printf("m22: %4.2f\n", m22);
    printf("m02: %4.2f\n", m02);
  }

  double a = m22;
  double b = (m02 - m12 - m22) * c1;
  double c = m12 * c1 * c1 - PDG_W_MASS * PDG_W_MASS;

  if (verbose) {
    printf("a: %4.2f\n", a);
    printf("b: %4.2f\n", b);
    printf("c: %4.2f\n", c);
  }

  double num = -1. * b + sqrt(b * b - 4 * a * c);
  double den = 2 * a;

  if (verbose) {
    printf("num: %4.2f\n", num);
    printf("den: %4.2f\n", den);
    printf("num/den: %4.2f\n", num/den);
  }

  return (num/den);
}

//--------------------------------------------------------------------
double fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2,
              double m12, double m22, double m02){
  double rat1 = pt1 * (1 - c1) / sigma1;
  double rat2 = pt2 * (1 - fc2(c1, m12, m22, m02)) / sigma2;

  return ( rat1 * rat1 + rat2 * rat2);
}

//--------------------------------------------------------------------
void minuitFunction(int&, double* , double &result, double par[], int){
  result=fchi2(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
}

// This function calculates the hadronic chi2 - SNT version
double calculateChi2(std::vector<LorentzVector>& jets, std::vector<float>& sigma_jets, std::vector<bool>& btag){

	assert(jets.size() == sigma_jets.size());
	assert(jets.size() == btag.size());

	//check at most first 6 jets
	int n_jets = jets.size();
	if (n_jets>6) n_jets = 6;
	//consider at least 3 jets
	if (n_jets<3) return 999999.;

	std::vector<int> v_i, v_j;
	std::vector<double> v_k1, v_k2;
	for ( int i=0; i<n_jets; ++i )
		for ( int j=i+1; j<n_jets; ++j ){

			//
			//  W
			//
			LorentzVector hadW = jets[i] + jets[j];

			//
			//  W Mass Constraint.
			//
			TFitter *minimizer = new TFitter();
			double p1 = -1;

			minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);
			minimizer->SetFCN(minuitFunction);
			minimizer->SetParameter(0 , "c1"     , 1.1             , 1 , 0 , 0);
			minimizer->SetParameter(1 , "pt1"    , 1.0             , 1 , 0 , 0);
			minimizer->SetParameter(2 , "sigma1" , sigma_jets[i]   , 1 , 0 , 0);
			minimizer->SetParameter(3 , "pt2"    , 1.0             , 1 , 0 , 0);
			minimizer->SetParameter(4 , "sigma2" , sigma_jets[j]   , 1 , 0 , 0);
			minimizer->SetParameter(5 , "m12"    , jets[i].mass2() , 1 , 0 , 0);
			minimizer->SetParameter(6 , "m22"    , jets[j].mass2() , 1 , 0 , 0);
			minimizer->SetParameter(7 , "m02"    , hadW.mass2()    , 1 , 0 , 0);

			for (unsigned int k = 1; k < 8; k++)
				minimizer->FixParameter(k);

			minimizer->ExecuteCommand("SIMPLEX", 0, 0);
			minimizer->ExecuteCommand("MIGRAD", 0, 0);

			double c1 = minimizer->GetParameter(0);
			if (c1!=c1) {
				std::cout<<"[PartonCombinatorics::recoHadronicTop] ERROR: c1 parameter is NAN! Skipping this parton combination"
					<<endl;
				continue;
			}
			double c2 = fc2(c1, jets[i].mass2(), jets[j].mass2(), hadW.mass2());

			delete minimizer;


			//     * W Mass check :)
			//     *  Never trust a computer you can't throw out a window.
			//      *  - Steve Wozniak

			// cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = "
			// 	   << ((jets[i] * c1) + (jets[j] * c2)).mass() << endl;

			v_i.push_back(i);
			v_j.push_back(j);
			v_k1.push_back(c1);
			v_k2.push_back(c2);
		}

	//Apply b-consistency requirement
	int n_btag = 0;
	for( int i = 0 ; i < n_jets ; i++ )
		if( btag.at(i) ) n_btag++;

	double chi2min = 99999.;

	//consider b-jet in leading 3 jets
	for ( int b=0; b<n_jets; ++b ) {    

		//if not tagged, consider only 3 leading jets
		if( (!btag.at(b)) && b>2 ) continue;

		//require b-tagging if have more than 1 b-tag
		if( n_btag>1 && (!btag.at(b)) ) continue;
		double pt_b = jets[b].Pt();

		for (unsigned int w = 0; w < v_i.size() ; ++w ) {
			int i = v_i[w];
			int j = v_j[w];
			if ( i==b || j==b ) continue;
			//count number of b-tagged Ws
			int nwb = 0;
			if (btag.at(i)) nwb++;
			if (btag.at(j)) nwb++;
			//no btagged jets in W if have few btags
			if ( n_btag<3  && nwb>0 ) continue;
			//In 3 b-tag case, allow for 1 W jet to be tagged
			// If have more b-tags then btagging information not useful
			if ( n_btag==3 && nwb>1 ) continue;

			double pt_w1 = jets[i].Pt();
			double pt_w2 = jets[j].Pt();

			///
			//  W Mass.
			///
			LorentzVector hadW = jets[i] + jets[j];
			double massW = hadW.mass();

			double c1 = v_k1[w];
			double c2 = v_k2[w];

			///
			// Top Mass.
			///
			LorentzVector hadT = (jets[i] * c1) + (jets[j] * c2) + jets[b];
			double massT = hadT.mass();

			double pt_w = hadW.Pt();
			double sigma_w2 = pow(pt_w1*sigma_jets[i], 2)
				+ pow(pt_w2*sigma_jets[j], 2);
			double smw2 = (1. + 2.*pow(pt_w,2)/pow(massW,2))*sigma_w2;
			double pt_t = hadT.Pt();
			double sigma_t2 = pow(c1*pt_w1*sigma_jets[i],2)
				+ pow(c2*pt_w2*sigma_jets[j],2)
				+ pow(pt_b*sigma_jets[b],2);
			double smtop2 = (1. + 2.*pow(pt_t,2)/pow(massT,2))*sigma_t2;

			double c_chi2 = pow(massT-PDG_TOP_MASS, 2)/smtop2
				+ pow(massW-PDG_W_MASS, 2)/smw2;
			if (c_chi2<chi2min) chi2min = c_chi2;

		}
	}

	return chi2min;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// m2tw_bisect //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
namespace mt2w_bisect
{

mt2w::mt2w(double upper_bound, double error_value, double scan_step)
{
   solved = false;
   momenta_set = false;
   mt2w_b  = 0.;  // The result field.  Start it off at zero.
   this->upper_bound = upper_bound;  // the upper bound of search for MT2W, default value is 500 GeV 
   this->error_value = error_value;  // if we couldn't find any compatible region below the upper_bound, output mt2w = error_value;
   this->scan_step = scan_step;    // if we need to scan to find the compatible region, this is the step of the scan
}

double mt2w::get_mt2w()
{
   if (!momenta_set)
   {
       std::cout <<" Please set momenta first!" << endl;
       return error_value;
   }
        
   if (!solved) mt2w_bisect();
   return mt2w_b;
}


void mt2w::set_momenta(double *pl, double *pb1, double *pb2, double* pmiss)
{
   // Pass in pointers to 4-vectors {E, px, py, px} of doubles.  
   // and pmiss must have [1] and [2] components for x and y.  The [0] component is ignored.
   set_momenta(pl[0],  pl[1],  pl[2],  pl[3],
               pb1[0], pb1[1], pb1[2], pb1[3],
               pb2[0], pb2[1], pb2[2], pb2[3],
               pmiss[1], pmiss[2]);
}



void mt2w::set_momenta(double El,  double plx,  double ply,  double plz,
                       double Eb1, double pb1x, double pb1y, double pb1z,
                       double Eb2, double pb2x, double pb2y, double pb2z,
                       double pmissx, double pmissy)
{
   solved = false;     //reset solved tag when momenta are changed.
   momenta_set = true;
	
	double msqtemp;   //used for saving the mass squared temporarily

//l is the visible lepton
	
	this->El  = El;
	this->plx = plx;
	this->ply = ply;
	this->plz = plz;
	
	Elsq = El*El;
	
	msqtemp = El*El-plx*plx-ply*ply-plz*plz;
	if (msqtemp > 0.0) {mlsq = msqtemp;}
	else {mlsq = 0.0;}                           //mass squared can not be negative
	ml = sqrt(mlsq);                             // all the input masses are calculated from sqrt(p^2)
	
//b1 is the bottom on the same side as the visible lepton
	
	this->Eb1  = Eb1;
	this->pb1x = pb1x;
	this->pb1y = pb1y;
	this->pb1z = pb1z;
	
	Eb1sq = Eb1*Eb1;
	
	msqtemp = Eb1*Eb1-pb1x*pb1x-pb1y*pb1y-pb1z*pb1z;
	if (msqtemp > 0.0) {mb1sq = msqtemp;}
	else {mb1sq = 0.0;}                          //mass squared can not be negative
	mb1 = sqrt(mb1sq);                           // all the input masses are calculated from sqrt(p^2)
	
//b2 is the other bottom (paired with the invisible W)
	
	this->Eb2  = Eb2;
	this->pb2x = pb2x;
	this->pb2y = pb2y;
	this->pb2z = pb2z;
	
	Eb2sq = Eb2*Eb2;
	
	msqtemp = Eb2*Eb2-pb2x*pb2x-pb2y*pb2y-pb2z*pb2z;
	if (msqtemp > 0.0) {mb2sq = msqtemp;}
	else {mb2sq = 0.0;}                          //mass squared can not be negative
	mb2 = sqrt(mb2sq);                           // all the input masses are calculated from sqrt(p^2)

	
//missing pt	
	
	
   this->pmissx = pmissx; 
   this->pmissy = pmissy;
	
//set the values of masses
	
	mv = 0.0;   //mass of neutrino
	mw = 80.4;  //mass of W-boson


//precision?	

   if (ABSOLUTE_PRECISION > 100.*RELATIVE_PRECISION) precision = ABSOLUTE_PRECISION;
   else precision = 100.*RELATIVE_PRECISION;
}


void mt2w::mt2w_bisect()
{
  
   
   solved = true;
   std::cout.precision(11);

	// In normal running, mtop_high WILL be compatible, and mtop_low will NOT.
	double mtop_high = upper_bound; //set the upper bound of the search region
	double mtop_low;                //the lower bound of the search region is best chosen as m_W + m_b

	if (mb1 >= mb2) {mtop_low = mw + mb1;}
	else {mtop_low = mw + mb2;}
	
	// The following if and while deal with the case where there might be a compatable region
	// between mtop_low and 500 GeV, but it doesn't extend all the way up to 500.
	// 
	
	// If our starting high guess is not compatible, start the high guess from the low guess...
    if (teco(mtop_high)==0) {mtop_high = mtop_low;}
	
	// .. and scan up until a compatible high bound is found.
	//We can also raise the lower bound since we scaned over a region that is not compatible
	while (teco(mtop_high)==0 && mtop_high < upper_bound + 2.*scan_step) {

		mtop_low=mtop_high;
		mtop_high = mtop_high + scan_step;
	}
	
	// if we can not find a compatible region under the upper bound, output the error value
	if (mtop_high > upper_bound) {
		mt2w_b = error_value;
		return;
	}
	
    // Once we have an compatible mtop_high, we can find mt2w using bisection method
   while(mtop_high - mtop_low > precision)
   {
      double mtop_mid,teco_mid;
      //bisect
      mtop_mid = (mtop_high+mtop_low)/2.;
      teco_mid = teco(mtop_mid);
      
	   if(teco_mid == 0) {mtop_low  = mtop_mid;}
	   else {mtop_high  = mtop_mid;}
	   
   }
   mt2w_b = mtop_high;   //output the value of mt2w
   return;
}


// for a given event, teco ( mtop ) gives 1 if trial top mass mtop is compatible, 0 if mtop is not.
	
int mt2w::teco(  double mtop)
{
	
//first test if mtop is larger than mb+mw	
	
	if (mtop < mb1+mw || mtop < mb2+mw) {return 0;}

//define delta for convenience, note the definition is different from the one in mathematica code by 2*E^2_{b2}
	  	
	double ETb2sq = Eb2sq - pb2z*pb2z;  //transverse energy of b2
	double delta = (mtop*mtop-mw*mw-mb2sq)/(2.*ETb2sq);
	
	
//del1 and del2 are \Delta'_1 and \Delta'_2 in the notes eq. 10,11
	
	double del1 = mw*mw - mv*mv - mlsq;
	double del2 = mtop*mtop - mw*mw - mb1sq - 2*(El*Eb1-plx*pb1x-ply*pb1y-plz*pb1z);
	
// aa bb cc are A B C in the notes eq.15
	
	double aa = (El*pb1x-Eb1*plx)/(Eb1*plz-El*pb1z);
	double bb = (El*pb1y-Eb1*ply)/(Eb1*plz-El*pb1z);
	double cc = (El*del2-Eb1*del1)/(2.*Eb1*plz-2.*El*pb1z);
	
  
//calculate coefficients for the two quadratic equations (ellipses), which are
//
//  a1 x^2 + 2 b1 x y + c1 y^2 + 2 d1 x + 2 e1 y + f1 = 0 ,  from the 2 steps decay chain (with visible lepton)
//
//  a2 x^2 + 2 b2 x y + c2 y^2 + 2 d2 x + 2 e2 y + f2 <= 0 , from the 1 stop decay chain (with W missing)
//
//  where x and y are px and py of the neutrino on the visible lepton chain

	a1 = Eb1sq*(1.+aa*aa)-(pb1x+pb1z*aa)*(pb1x+pb1z*aa);
	b1 = Eb1sq*aa*bb - (pb1x+pb1z*aa)*(pb1y+pb1z*bb);
	c1 = Eb1sq*(1.+bb*bb)-(pb1y+pb1z*bb)*(pb1y+pb1z*bb);
	d1 = Eb1sq*aa*cc - (pb1x+pb1z*aa)*(pb1z*cc+del2/2.0);
	e1 = Eb1sq*bb*cc - (pb1y+pb1z*bb)*(pb1z*cc+del2/2.0);
	f1 = Eb1sq*(mv*mv+cc*cc) - (pb1z*cc+del2/2.0)*(pb1z*cc+del2/2.0);
	
//  First check if ellipse 1 is real (don't need to do this for ellipse 2, ellipse 2 is always real for mtop > mw+mb)
	
    double det1 = (a1*(c1*f1 - e1*e1) - b1*(b1*f1 - d1*e1) + d1*(b1*e1-c1*d1))/(a1+c1);
	
	if (det1 > 0.0) {return 0;}
	
//coefficients of the ellptical region
	
	a2 = 1-pb2x*pb2x/(ETb2sq);
	b2 = -pb2x*pb2y/(ETb2sq);
	c2 = 1-pb2y*pb2y/(ETb2sq);
	
	// d2o e2o f2o are coefficients in the p2x p2y plane (p2 is the momentum of the missing W-boson)
	// it is convenient to calculate them first and transfer the ellipse to the p1x p1y plane
	d2o = -delta*pb2x;
	e2o = -delta*pb2y;
	f2o = mw*mw - delta*delta*ETb2sq;
	
	d2 = -d2o -a2*pmissx -b2*pmissy;
	e2 = -e2o -c2*pmissy -b2*pmissx;
	f2 = a2*pmissx*pmissx + 2*b2*pmissx*pmissy + c2*pmissy*pmissy + 2*d2o*pmissx + 2*e2o*pmissy + f2o;
	
//find a point in ellipse 1 and see if it's within the ellipse 2, define h0 for convenience
    double x0, h0, y0, r0;
	x0 = (c1*d1-b1*e1)/(b1*b1-a1*c1);
    h0 = (b1*x0 + e1)*(b1*x0 + e1) - c1*(a1*x0*x0 + 2*d1*x0 + f1);
	if (h0 < 0.0) {return 0;}  // if h0 < 0, y0 is not real and ellipse 1 is not real, this is a redundant check.
	y0 = (-b1*x0 -e1 + sqrt(h0))/c1;
	r0 = a2*x0*x0 + 2*b2*x0*y0 + c2*y0*y0 + 2*d2*x0 + 2*e2*y0 + f2;
	if (r0 < 0.0) {return 1;}  // if the point is within the 2nd ellipse, mtop is compatible
	
	
//obtain the coefficients for the 4th order equation 
//devided by Eb1^n to make the variable dimensionless
   long double A4, A3, A2, A1, A0;

	A4 = 
	-4*a2*b1*b2*c1 + 4*a1*b2*b2*c1 +a2*a2*c1*c1 + 
	4*a2*b1*b1*c2 - 4*a1*b1*b2*c2 - 2*a1*a2*c1*c2 + 
	a1*a1*c2*c2;  
	
	A3 =
	(-4*a2*b2*c1*d1 + 8*a2*b1*c2*d1 - 4*a1*b2*c2*d1 - 4*a2*b1*c1*d2 + 
	 8*a1*b2*c1*d2 - 4*a1*b1*c2*d2 - 8*a2*b1*b2*e1 + 8*a1*b2*b2*e1 + 
	 4*a2*a2*c1*e1 - 4*a1*a2*c2*e1 + 8*a2*b1*b1*e2 - 8*a1*b1*b2*e2 - 
     4*a1*a2*c1*e2 + 4*a1*a1*c2*e2)/Eb1;
	
	
	A2 =
	(4*a2*c2*d1*d1 - 4*a2*c1*d1*d2 - 4*a1*c2*d1*d2 + 4*a1*c1*d2*d2 - 
	 8*a2*b2*d1*e1 - 8*a2*b1*d2*e1 + 16*a1*b2*d2*e1 + 
	 4*a2*a2*e1*e1 + 16*a2*b1*d1*e2 - 8*a1*b2*d1*e2 - 
	 8*a1*b1*d2*e2 - 8*a1*a2*e1*e2 + 4*a1*a1*e2*e2 - 4*a2*b1*b2*f1 + 
	 4*a1*b2*b2*f1 + 2*a2*a2*c1*f1 - 2*a1*a2*c2*f1 + 
     4*a2*b1*b1*f2 - 4*a1*b1*b2*f2 - 2*a1*a2*c1*f2 + 2*a1*a1*c2*f2)/Eb1sq;
	
	A1 =
	(-8*a2*d1*d2*e1 + 8*a1*d2*d2*e1 + 8*a2*d1*d1*e2 - 8*a1*d1*d2*e2 - 
	 4*a2*b2*d1*f1 - 4*a2*b1*d2*f1 + 8*a1*b2*d2*f1 + 4*a2*a2*e1*f1 - 
	 4*a1*a2*e2*f1 + 8*a2*b1*d1*f2 - 4*a1*b2*d1*f2 - 4*a1*b1*d2*f2 - 
     4*a1*a2*e1*f2 + 4*a1*a1*e2*f2)/(Eb1sq*Eb1);
	
	A0 =
	(-4*a2*d1*d2*f1 + 4*a1*d2*d2*f1 + a2*a2*f1*f1 + 
	 4*a2*d1*d1*f2 - 4*a1*d1*d2*f2 - 2*a1*a2*f1*f2 + 
     a1*a1*f2*f2)/(Eb1sq*Eb1sq);
	
   long  double A0sq, A1sq, A2sq, A3sq, A4sq;
   A0sq = A0*A0;
   A1sq = A1*A1;
   A2sq = A2*A2;
   A3sq = A3*A3;
   A4sq = A4*A4;
   
   long double B3, B2, B1, B0;
   B3 = 4*A4;
   B2 = 3*A3;
   B1 = 2*A2;
   B0 = A1;
   
   long double C2, C1, C0;
   C2 = -(A2/2 - 3*A3sq/(16*A4));
   C1 = -(3*A1/4. -A2*A3/(8*A4));
   C0 = -A0 + A1*A3/(16*A4);
   
   long double D1, D0;
   D1 = -B1 - (B3*C1*C1/C2 - B3*C0 -B2*C1)/C2;
   D0 = -B0 - B3 *C0 *C1/(C2*C2)+ B2*C0/C2;
   
   long double E0;
   E0 = -C0 - C2*D0*D0/(D1*D1) + C1*D0/D1;
   
   long  double t1,t2,t3,t4,t5;
//find the coefficients for the leading term in the Sturm sequence  
   t1 = A4;
   t2 = A4;
   t3 = C2;
   t4 = D1;
   t5 = E0;
 

//The number of solutions depends on diffence of number of sign changes for x->Inf and x->-Inf
   int nsol;
   nsol = signchange_n(t1,t2,t3,t4,t5) - signchange_p(t1,t2,t3,t4,t5);

//Cannot have negative number of solutions, must be roundoff effect
   if (nsol < 0) nsol = 0;
	
    int out;
	if (nsol == 0) {out = 0;}  //output 0 if there is no solution, 1 if there is solution
	else {out = 1;}

   return out;
  
}  

inline int mt2w::signchange_n( long double t1, long double t2, long double t3, long double t4, long double t5)
{
   int nsc;
   nsc=0;
   if(t1*t2>0) nsc++;
   if(t2*t3>0) nsc++;
   if(t3*t4>0) nsc++;
   if(t4*t5>0) nsc++;
   return nsc;
}
inline int mt2w::signchange_p( long double t1, long double t2, long double t3, long double t4, long double t5)
{
   int nsc;
   nsc=0;
   if(t1*t2<0) nsc++;
   if(t2*t3<0) nsc++;
   if(t3*t4<0) nsc++;
   if(t4*t5<0) nsc++;
   return nsc;
}

}//end namespace mt2w_bisect
