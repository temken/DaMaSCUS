#include "ER_LUX.hpp"

#include <iostream>

#include "Physical_Parameters.hpp"
#include "ER_General.hpp"

//Analytic
double R_LUX_A(double rhoDM,double mChi,double sigma,double vE)
		{
			double minDev=0.1;
			//Mean value from Detector Efficiency taken from fig. 2 from https://arxiv.org/pdf/1608.07648.pdf)
			double efficiency = 0.580638; 
			double A=131.0;
			double Ethreshold = 1.1*keV;
			double Emax = 60*keV;
			//Numerical integration with trapezoidal rule
				double integral_old=1.0;
				double integral_new=0.0;
				double interval = Emax-Ethreshold;
				double dE = interval;
				while(100.0*abs(integral_old-integral_new)/abs(integral_old)>minDev)
				{
					dE/=10;
					integral_old=integral_new;
					int n = interval/dE;
					double sum=0.0;
					for(int k=1;k<n;k++)
					{
						double E = Ethreshold+k*interval/n;
						sum+=dRdErA(E,1,A,rhoDM,mChi,sigma,vE);
					}
					integral_new=interval/n*efficiency*(dRdErA(Ethreshold,1,A,rhoDM,mChi,sigma,vE)/2.0+dRdErA(Emax,1,A,rhoDM,mChi,sigma,vE)/2.0+sum);
				}
				
			return integral_new;
		}

//Monte Carlo
	std::vector<double> R_LUX_MC(std::vector<std::vector<double>> &drdehisto)
	{
		std::vector<double> output;
		double minDev=0.1;
		double efficiency = 0.580638; //Mean value from Detector Efficiency taken from fig. 2 from https://arxiv.org/pdf/1608.07648.pdf)
		double Ethreshold = 1.1*keV;
		double Emax = 60*keV;
		//Numerical integration with trapezoidal rule
			double integral_old=1.0;
			double integral_new=0.0;
			double error;
			double interval = Emax-Ethreshold;
			double dE = interval;
			while(100.0*abs(integral_old-integral_new)/abs(integral_old)>minDev)
			{
				dE/=2;
				integral_old=integral_new;
				int n = interval/dE;
				double sum=0.0;
				double UpperBoundSum=0.0;
				for(int k=1;k<n;k++)
				{
					double E = Ethreshold+k*interval/n;
					std::vector<double> dR = dRdEMC(E,drdehisto);
					sum 			+=efficiency*dR[0];
					UpperBoundSum 	+=efficiency*dR[1];
				}
				integral_new=	interval/n*(efficiency*dRdEMC(Ethreshold,drdehisto)[0]/2.0+efficiency*dRdEMC(Emax,drdehisto)[0]/2.0+sum);
				error = 		interval/n*(efficiency*dRdEMC(Ethreshold,drdehisto)[1]/2.0+efficiency*dRdEMC(Emax,drdehisto)[1]/2.0+UpperBoundSum);
			}
			output.push_back(integral_new);
			output.push_back(error);
		return output;
	}