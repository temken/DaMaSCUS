// Contains non-detector-specific functions related with direct detection rates.

#ifndef __ER_General_hpp_
#define __ER_General_hpp_

#include <vector>

// Analytic

// Velocity necessary for a DM particle of mass mChi to recoil a nucleus of mass number A with a given recoil energy.
extern double vMinimal(double RecoilEnergy, double mChi, double A);
// Maximum recoil energy a DM particle of mass mChi and speed v can cause on a nucleus of mass number A
extern double ERMax(double v, double mChi, double A);
// Eta-Function: Integral of f(v)/v with lower integration limit v_Min
extern double EtaFunctionA(double vMin, double vE);
// Direct detection recoil spectrum
extern double dRdErA(double ER, double X, double A, double rhoDM, double mChi, double sigma, double vE);

// Monte Carlo

// Eta function
extern std::vector<std::vector<double>> EtaHistogram(std::vector<std::vector<double>>& velocityhistogram);
extern std::vector<double> EtaFunctionMC(double vMin, std::vector<std::vector<double>>& etaarray);
// Recoil spectrum
extern std::vector<std::vector<double>> dRdEHistogram(double X, double A, std::vector<double> rhoDM, double mChi, double sigma, std::vector<std::vector<double>>& etahistogram);
extern std::vector<double> dRdEMC(double ER, std::vector<std::vector<double>>& drdehisto);
#endif