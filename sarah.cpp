#include <iostream>
#include <fstream>
#include "analysis_manager.hh"

int main(){

    const std::string & fitmode("combined");
    const std::vector<double> signalyields(in_signalyields,in_signalyields+1);
    const double  signal_uncertainty(0.1);
    const std::vector<double>  bgyields(in_signalyields,in_signalyields+1);
    const std::vector<double>  bguncert(in_signalyields,in_signalyields+1);
    const std::vector<int>  datayields(in_signalyields,in_signalyields+1);
    const bool  calculateR(true);
    AnalysisManager jad;
    std::vector<fitparams> my_params=jad.LimitCode( fitmode,
				    signalyields,
				    signal_uncertainty,
				    bgyields,
				    bguncert,
				    datayields,
				    calculateR);
    std::cout << my_params[0].cls << std::endl;
    return 0;
}
