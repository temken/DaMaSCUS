#include "ER_CRESSTII.hpp"

#include <iostream>
#include <cmath>

#include "General_Utilities.hpp"
#include "Physical_Parameters.hpp"
#include "ER_General.hpp"

//Detector response functions
	double Efficiency(double E)
	{
		double output=0.0;
		double epsilon[28][2]={
			{0.3*keV	,0.122},
			{0.351*keV	,0.171},
			{0.400*keV	,0.223},
			{0.448*keV	,0.274},
			{0.498*keV	,0.322},
			{0.549*keV	,0.369},
			{0.597*keV	,0.409},
			{0.647*keV	,0.444},
			{0.700*keV	,0.474},
			{0.748*keV	,0.498},
			{0.796*keV	,0.523},
			{0.847*keV	,0.543},
			{0.896*keV	,0.561},
			{0.946*keV	,0.578},
			{0.995*keV	,0.591},
			{1.097*keV	,0.613},
			{1.196*keV	,0.633},
			{1.394*keV	,0.656},
			{1.598*keV	,0.672},
			{1.799*keV	,0.681},
			{1.998*keV	,0.687},
			{2.196*keV	,0.692},
			{2.396*keV	,0.694},
			{2.595*keV	,0.695},
			{2.800*keV	,0.698},
			{3.000*keV	,0.698},
			{3.497*keV	,0.697},
			{3.999*keV	,0.697}
		};
		if(E<epsilon[0][0] || E>epsilon[27][0]) output=0.0;
		else
		{
			for(int i=0;i<28;i++)
			{
				if(epsilon[i][0]>=E)
				{
					double slope =(epsilon[i][1]-epsilon[i-1][1])/(epsilon[i][0]-epsilon[i-1][0]); 
					output = epsilon[i-1][1] + slope * (E-epsilon[i-1][0]);
					break;
				}
			}
		}
		return output;
	}

	double DetectorResponse(double ER,double E)
	{
		double sigmaE2=pow(62*eV,2);
		return 0.5*Efficiency(E)/sqrt(2*M_PI*sigmaE2)*exp(-(ER-E)*(ER-E)/2.0/sigmaE2);
	}

//Analytic
	double dRdE_CRESSTII_A(double E,double rhoDM,double mChi,double sigma,double vE)
	{
		double minDev=0.1;
		double ERmin = 0.0;
		double sigmaE=62*eV;
		int A[]={16,40,184};
		double X[]={0.22,0.14,0.64};
		//Find the max. recoil energy
			std::vector<double> ERmaxV;
			for(int i=0;i<3;i++)	ERmaxV.push_back(2*Mu(mChi,A[i]*mNucleon)*Mu(mChi,A[i]*mNucleon)*pow(vE+vesc,2.0)/A[i]/mNucleon+3.0*sigmaE);
			// for(int i=0;i<3;i++)	cout <<ERmaxV[i] <<endl;
			double ERmax = *std::max_element( ERmaxV.begin(), ERmaxV.end() );
			// cout <<ERmax <<endl;
		//Numerical integration with trapezoidal rule
			double integral_old=1.0;
			double integral_new=0.0;
			double interval = ERmax-ERmin;
			double dER = interval;
			while(100.0*abs(integral_old-integral_new)/abs(integral_old)>minDev)
			{
				dER/=2;
				integral_old=integral_new;
				int n = interval/dER;
				double sum=0.0;
				for(int k=1;k<n;k++)
				{
					double ER = ERmin+k*interval/n;
					for(int i=0;i<3;i++)
					{
						sum+=DetectorResponse(ER,E)*dRdErA(ER,X[i],A[i],rhoDM,mChi,sigma,vE);
					}
					
				}
				integral_new=interval/n*sum;
				for(int i=0;i<3;i++)
				{
					integral_new+=interval/n*(DetectorResponse(ERmin,E)*dRdErA(ERmin,X[i],A[i],rhoDM,mChi,sigma,vE)/2.0+DetectorResponse(ERmax,E)*dRdErA(ERmax,X[i],A[i],rhoDM,mChi,sigma,vE)/2.0);
				}
			}
			
		return integral_new;
	}

	double R_CRESSTII_A(double rhoDM,double mChi,double sigma,double vE)
	{
		double minDev=0.1;
		double integral_new=0.0;
		//Integration Limits
			double Emin = 0.3*keV;
			double Emax = 4*keV;
			//Integration
			double integral_old=1.0;
			double interval = Emax-Emin;
			double dE = interval;
			while(100.0*abs(integral_old-integral_new)/abs(integral_old)>minDev)
			{
				dE/=2;
				integral_old=integral_new;
				int n = interval/dE;
				double sum=0.0;
				for(int k=1;k<n;k++)
				{
					double E = Emin+k*dE;
					sum+=dRdE_CRESSTII_A(E,rhoDM,mChi,sigma,vE);
				}
				integral_new=interval/n*(dRdE_CRESSTII_A(Emin,rhoDM,mChi,sigma,vE)/2.0+dRdE_CRESSTII_A(Emax,rhoDM,mChi,sigma,vE)/2.0+sum);
			}
		return integral_new;
	}

//Monte Carlo
	std::vector<double> dRdE_CRESSTII_MC(double E,std::vector<std::vector<double>> &drdehistoO,std::vector<std::vector<double>> &drdehistoCa,std::vector<std::vector<double>> &drdehistoW,double vE)
	{
		std::vector<double> output;
		double minDev=0.1;
		double ERmin = 0.0;
		int A[]={16,40,184};
		double sigmaE=62*eV;
		//Find the max. recoil energy
			std::vector<double> ERmaxV;
			for(int i=0;i<3;i++)	ERmaxV.push_back(2*Mu(mChi,A[i]*mNucleon)*Mu(mChi,A[i]*mNucleon)*pow(vE+vesc,2.0)/A[i]/mNucleon+3.0*sigmaE);
			double ERmax = *std::max_element( ERmaxV.begin(), ERmaxV.end() );
		//Numerical integration with trapezoidal rule
			double integral_old=1.0;
			double integral_new=0.0;
			double error;
			double interval = ERmax-ERmin;
			double dER = interval;
			double Dev=100.0*abs(integral_old-integral_new)/abs(integral_old);
			while(Dev>minDev)
			{
				dER/=2;
				integral_old=integral_new;
				int n = ceil(interval/dER);
				double sum=0.0;
				double UpperBoundSum=0.0;
				for(int k=1;k<n;k++)
				{
					double ER = ERmin+k*interval/n;
					double response=DetectorResponse(ER,E);
					std::vector<double> dRo = dRdEMC(ER,drdehistoO);
					std::vector<double> dRca = dRdEMC(ER,drdehistoCa);					
					std::vector<double> dRw = dRdEMC(ER,drdehistoW);					
					sum 			+=	response*(dRo[0]+dRca[0]+dRw[0]);
					UpperBoundSum 	+=	response*(dRo[1]+dRca[1]+dRw[1]);
				}
				std::vector<double> dRminO = dRdEMC(ERmin,drdehistoO);
				std::vector<double> dRmaxO = dRdEMC(ERmax,drdehistoO);
				std::vector<double> dRminCa = dRdEMC(ERmin,drdehistoCa);
				std::vector<double> dRmaxCa = dRdEMC(ERmax,drdehistoCa);
				std::vector<double> dRminW = dRdEMC(ERmin,drdehistoW);
				std::vector<double> dRmaxW = dRdEMC(ERmax,drdehistoW);
				integral_new=	dER*( DetectorResponse(ERmin,E)*(dRminO[0]/2.0+dRminCa[0]/2.0+dRminW[0]/2.0)+DetectorResponse(ERmax,E)*(dRmaxO[0]/2.0+dRmaxCa[0]/2.0+dRmaxW[0]/2.0) +sum);
				error=			dER*( DetectorResponse(ERmin,E)*(dRminO[1]/2.0+dRminCa[1]/2.0+dRminW[1]/2.0)+DetectorResponse(ERmax,E)*(dRmaxO[1]/2.0+dRmaxCa[1]/2.0+dRmaxW[1]/2.0) +UpperBoundSum);
	
				Dev=100.0*abs(integral_old-integral_new)/abs(integral_old);
			}
		output.push_back(integral_new);
		output.push_back(error);
			return output;
	}

	std::vector<double> R_CRESSTII_MC(std::vector<std::vector<double>> &drdehistoO,std::vector<std::vector<double>> &drdehistoCa,std::vector<std::vector<double>> &drdehistoW,double mChi,double vE)
	{
		std::vector<double> output;
		double minDev=0.1;
		double integral_new=0.0;
		//Integration Limits
			double Emin = 0.3*keV;
			double Emax = 4*keV;
		//Integration
			double error;
			double integral_old=1.0;
			double interval = Emax-Emin;
			double dE = interval;
			while(100.0*abs(integral_old-integral_new)/abs(integral_old)>minDev)
			{
				dE/=2;
				integral_old=integral_new;
				int n = ceil(interval/dE);
				double sum=0.0;
				double UpperBoundSum=0.0;
				for(int k=1;k<n;k++)
				{
					double E = Emin+k*dE;
					std::vector<double> dR=dRdE_CRESSTII_MC(E,drdehistoO,drdehistoCa,drdehistoW,vE);
					sum 			+=	dR[0];
					UpperBoundSum 	+=	dR[1];
				}
				std::vector<double> dRmin=dRdE_CRESSTII_MC(Emin,drdehistoO,drdehistoCa,drdehistoW,vE);
				std::vector<double> dRmax=dRdE_CRESSTII_MC(Emax,drdehistoO,drdehistoCa,drdehistoW,vE);
				integral_new 	= 	dE*(dRmin[0]/2.0+dRmax[0]/2.0+sum);
				error 			=	dE*(dRmin[1]/2.0+dRmax[1]/2.0+UpperBoundSum);
			}
			output.push_back(integral_new);
			output.push_back(error);
		return output;
	}