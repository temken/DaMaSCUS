#ifndef __ER_CRESSTII_hpp_
#define __ER_CRESSTII_hpp_

#include <vector>

// Detection efficiency and response function for a CRESST-II type detector
double Efficiency(double E);
double DetectorResponse(double ER, double E);

// Analytic
// Energy spectrum
extern double dRdE_CRESSTII_A(double E, double rhoDM, double mChi, double sigma, double vE);
// Event rates at CRESST-II type detector
extern double R_CRESSTII_A(double rhoDM, double mChi, double sigma, double vE);

// Monte Carlo
extern std::vector<double> dRdE_CRESSTII_MC(double E, std::vector<std::vector<double>>& drdehistoO, std::vector<std::vector<double>>& drdehistoCa, std::vector<std::vector<double>>& drdehistoW, double vE);
extern std::vector<double> R_CRESSTII_MC(std::vector<std::vector<double>>& drdehistoO, std::vector<std::vector<double>>& drdehistoCa, std::vector<std::vector<double>>& drdehistoW, double mChi, double vE);

#endif