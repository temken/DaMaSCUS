#ifndef __ER_LUX_hpp_
#define __ER_LUX_hpp_

#include <vector>

//Analytic
	//Simplified  calculation of event rates at a LUX-type detector
		extern double R_LUX_A(double rhoDM,double mChi,double sigma,double vE);
//Monte Carlo
		extern std::vector<double> R_LUX_MC(std::vector<std::vector<double>> &drdehisto);
#endif