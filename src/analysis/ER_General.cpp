#include "ER_General.hpp"

#include <iostream>

#include "Physical_Parameters.hpp"

// Analytic
double vMinimal(double RecoilEnergy, double mChi, double A)
{
	return sqrt(mNucleon * A * RecoilEnergy / 2.0 / pow(Mu(mChi, A * mNucleon), 2.0));
}

double ERMax(double v, double mChi, double A)
{
	return 2.0 * v * v * pow(Mu(mChi, A * mNucleon), 2.0) / A / mNucleon;
}

double EtaFunctionA(double vMin, double vE)
{
	double xMin = vMin / v0;
	double xEsc = vesc / v0;
	double xE	= vE / v0;
	if(xMin > xE + xEsc)
		return 0.0;
	else if(xMin > abs(xE - xEsc))
		return pow(M_PI, 1.5) * v0 * v0 / 2.0 / Nesc / xE * (erf(xEsc) - erf(xMin - xE) - 2.0 / sqrt(M_PI) * (xE + xEsc - xMin) * exp(-xEsc * xEsc));
	else if(xEsc > xE)
		return pow(M_PI, 1.5) * v0 * v0 / 2.0 / Nesc / xE * (erf(xMin + xE) - erf(xMin - xE) - 4.0 / sqrt(M_PI) * xE * exp(-xEsc * xEsc));
	else
		return 1.0 / v0 / xE;
}

double dRdErA(double ER, double X, double A, double rhoDM, double mChi, double sigma, double vE)
{
	return X / 2.0 * rhoDM / mChi * sigma * pow(A / Mu(mChi, mNucleon), 2.0) * EtaFunctionA(vMinimal(ER, mChi, A), vE);
}

// Monte Carlo
// MC Eta Function
std::vector<std::vector<double>> EtaHistogram(std::vector<std::vector<double>>& velocityhistogram)
{
	// Initialize Histogram vector
	double h = velocityhistogram[1][0] - velocityhistogram[0][0];
	int bins = velocityhistogram.size();
	std::vector<std::vector<double>> histogram(bins, std::vector<double>(4));
	for(int i = 0; i < bins; i++)
	{
		double vMin = i * h;
		// Errors
		histogram[i][0] = vMin;		 // x coordinate / bin position
		histogram[i][1] = 0.0;		 // y coordinate / bin height
		histogram[i][2] = h / 2.0;	 // error on x
		histogram[i][3] = 0.0;		 // error on y
		// Calculate the eta function by adding up the PDF histogram bins
		for(int j = i; j < bins; j++)
		{
			histogram[i][1] += velocityhistogram[j][1] / (j + 0.5);	  // the h cancels
			histogram[i][3] += pow(velocityhistogram[j][3] / (j + 0.5), 2.0) + pow(velocityhistogram[j][1] / (j + 0.5) / (j + 0.5), 2.0);
		}
		histogram[i][3] = sqrt(histogram[i][3]);
	}
	return histogram;
}
std::vector<double> EtaFunctionMC(double vMin, std::vector<std::vector<double>>& etahistogram)
{
	std::vector<double> output;
	double vMax = etahistogram.back()[0];
	if(vMin > vMax)
	{
		output.push_back(0.0);
		output.push_back(0.0);
	}
	else
	{
		double h		 = etahistogram[1][0] - etahistogram[0][0];
		unsigned int bin = vMin / h;
		output.push_back(etahistogram[bin][1]);
		output.push_back(etahistogram[bin][3]);
	}
	return output;
}
// MC Recoil Spectrum
std::vector<std::vector<double>> dRdEHistogram(double X, double A, std::vector<double> rhoDM, double mChi, double sigma, std::vector<std::vector<double>>& etahistogram)
{
	double prefactor = X / 2.0 / mChi * sigma * pow(A / Mu(mChi, mNucleon), 2.0);
	// Histogram parameter
	int bins = etahistogram.size();
	std::vector<std::vector<double>> histo(bins, std::vector<double>(4));
	// Compute the histogram
	for(int i = 0; i < bins; i++)
	{
		double h;
		if(i == 0)
			h = ERMax(etahistogram[i][0], mChi, A) / 2.0;
		else
			h = ERMax(etahistogram[i][0], mChi, A) - ERMax(etahistogram[i - 1][0], mChi, A);
		histo[i][0] = ERMax(etahistogram[i][0], mChi, A);
		histo[i][1] = prefactor * rhoDM[0] * etahistogram[i][1];
		histo[i][2] = h / 2.0;
		histo[i][3] = prefactor * sqrt(pow(rhoDM[1] * etahistogram[i][1], 2.0) + pow(rhoDM[0] * etahistogram[i][3], 2.0));
	}
	return histo;
}

std::vector<double> dRdEMC(double ER, std::vector<std::vector<double>>& drdehisto)
{
	std::vector<double> output;
	double ERmax = drdehisto.back()[0];
	if(ER > ERmax)
	{
		output.push_back(0.0);
		output.push_back(0.0);
	}
	else if(ER < 0)
	{
		cout << "Critical Error in dRdEMC: ER is negative." << endl;
		return output;
	}
	else
	{
		// Find the bin of ER
		unsigned int bin, leftbin, rightbin;
		for(bin = 0; bin < drdehisto.size(); bin++)
		{
			double e = drdehisto[bin][0] + drdehisto[bin][2];
			if(e > ER)
				break;
		}
		// Find the two bins to use for the interpolation
		if(ER < drdehisto.front()[0] || (ER > drdehisto[bin][0] && ER < drdehisto.back()[0]))
		{
			leftbin	 = bin;
			rightbin = bin + 1;
		}
		else if(ER > drdehisto.back()[0] || (ER < drdehisto[bin][0]))
		{
			leftbin	 = bin - 1;
			rightbin = bin;
		}
		else
		{
			cout << "Critical error in dRdEMC." << endl;
			return output;
		}
		// Linear Interpolation
		double slope	= (drdehisto[rightbin][1] - drdehisto[leftbin][1]) / (drdehisto[rightbin][0] - drdehisto[leftbin][0]);
		double varslope = (pow(drdehisto[rightbin][3], 2.0) + pow(drdehisto[leftbin][3], 2.0)) / pow(drdehisto[rightbin][0] - drdehisto[leftbin][0], 2.0);
		double drde		= drdehisto[bin][1] + (ER - drdehisto[bin][0]) * slope;
		double variance = pow(drdehisto[bin][3], 2.0) + pow(ER - drdehisto[bin][0], 2.0) * varslope;
		output.push_back(drde);
		output.push_back(sqrt(variance));
		// output.push_back(drdehisto[bin][3]);
	}
	return output;
}