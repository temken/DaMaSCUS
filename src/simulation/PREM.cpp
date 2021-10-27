#include "PREM.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

#include "General_Utilities.hpp"
#include "RN_Generators.hpp"

// PREM
// Radii
//							0			1			2			3			4			5			6			7			8			9
const double PREM_Radii[10] = {1221.5 * km, 3480 * km, 5701 * km, 5771 * km, 5971 * km, 6151 * km, 6346.6 * km, 6356 * km, 6368 * km, 6371 * km};
// Density Coefficients, we omit the ocean layer with density 1.02, since experiments are located under rocks, not water.
//					0			1			2			3			4			5			6			7		8		9		10
//					I.C.		O.C.		L.Mantle 	Trans. 1	Trans. 2	Trans. 3	LVZ&LID		Crust1	Crust2	Ocean	Space
const double a[11] = {13.0885, 12.5815, 7.9565, 5.3197, 11.2494, 7.1089, 2.6910, 2.9, 2.6, 2.6, 0};
const double b[11] = {0, -1.2638, -6.4761, -1.4836, -8.0298, -3.8045, 0.6924, 0, 0, 0, 0};
const double c[11] = {-8.8381, -3.6426, 5.5283, 0, 0, 0, 0, 0, 0, 0, 0};
const double d[11] = {0, -5.5281, -3.0807, 0, 0, 0, 0, 0, 0, 0, 0};
// PREM Functions
int PREM_Layer(double r)
{
	// The 10 mechanical Layers
	for(int i = 9; i >= 0; i--)
	{
		if(r > PREM_Radii[i])
		{
			return i + 1;
		}
	}
	// Inner Core:
	return 0;
}
double PREM_Mass_Density(double r)
{
	int i	 = PREM_Layer(r);
	double x = r / rEarth;
	return (a[i] + b[i] * x + c[i] * pow(x, 2) + d[i] * pow(x, 3)) * gram * pow(cm, -3);
}

// Earth Composition
const double Elements_Core[9][3] =
	{
		{26, 56, 0.855},   // Iron			Fe
		{14, 28, 0.06},	   // Silicon		Si
		{28, 58, 0.052},   // Nickel		Ni
		{16, 32, 0.019},   // Sulfur		S
		{24, 52, 0.009},   // Chromium		Cr
		{25, 55, 0.003},   // Manganese		Mn
		{15, 31, 0.002},   // Phosphorus	P
		{6, 12, 0.002},	   // Carbon		C
		{1, 1, 0.0006}	   // Hydrogen		H
};
const double Elements_Mantle[14][3] =
	{
		{8, 16, 0.440},		// Oxygen		O
		{12, 24, 0.228},	// Magnesium		Mg
		{14, 28, 0.21},		// Silicon		Si
		{26, 56, 0.0626},	// Iron			Fe
		{20, 40, 0.0253},	// Calcium		Ca
		{13, 27, 0.0235},	// Aluminium		Al
		{11, 23, 0.0027},	// Natrium		Na
		{24, 52, 0.0026},	// Chromium		Cr
		{28, 58, 0.002},	// Nickel		Ni
		{25, 55, 0.001},	// Manganese		Mn
		{16, 32, 0.0003},	// Sulfur		S
		{6, 12, 0.0001},	// Carbon		C
		{1, 1, 0.0001},		// Hydrogen		H
		{15, 31, 0.00009}	// Phosphorus	P
};

// Probability to scatter on a certain element. These are calculated in Initialize_PREM().
double Scatter_Probability_Core[9][3];
double Scatter_Probability_Mantle[14][3];

// calculates the space-independent prefactor and the scatternucleus probabilities for a given parameter point
// Mean Free Path PreFactor for F(q^2)=1. Including a form factor requires the calculation of these factors at every step of a simulation, because everything becomes velocity depending.
double g_PREM[2] = {0.0, 0.0};	 // Core Mantle
void Initialize_PREM(double mX, double sigma0)
{
	g_PREM[0] = 0.0;
	g_PREM[1] = 0.0;
	// 1. MFP Prefactor
	// core
	for(int i = 0; i < 9; i++)
	{
		g_PREM[0] += Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * sigmaSI(mX, sigma0, Elements_Core[i][1]);
	}
	// mantle
	for(int i = 0; i < 14; i++)
	{
		g_PREM[1] += Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * sigmaSI(mX, sigma0, Elements_Mantle[i][1]);
	}
	// 2. Probabilities. The probabilities are velocity depending if form factors are included.
	double total = 0.0;
	// core
	total = 0.0;
	for(int i = 0; i < 9; i++)
	{
		Scatter_Probability_Core[i][0] = Elements_Core[i][0];
		Scatter_Probability_Core[i][1] = Elements_Core[i][1];
		Scatter_Probability_Core[i][2] = Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * sigmaSI(mX, sigma0, Elements_Core[i][1]);
		total += Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * sigmaSI(mX, sigma0, Elements_Core[i][1]);
	}
	for(int i = 0; i < 9; i++)
	{
		Scatter_Probability_Core[i][2] /= total;
		// cout <<"Core \t" <<Scatter_Probability_Core[i][0] <<"\t"<<Scatter_Probability_Core[i][1]<<"\t"<<Scatter_Probability_Core[i][2] <<endl;
	}

	// mantle
	total = 0.0;
	for(int i = 0; i < 14; i++)
	{
		Scatter_Probability_Mantle[i][0] = Elements_Mantle[i][0];
		Scatter_Probability_Mantle[i][1] = Elements_Mantle[i][1];
		Scatter_Probability_Mantle[i][2] = Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * sigmaSI(mX, sigma0, Elements_Mantle[i][1]);
		total += Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * sigmaSI(mX, sigma0, Elements_Mantle[i][1]);
	}
	for(int i = 0; i < 14; i++)
	{
		Scatter_Probability_Mantle[i][2] /= total;
		// cout <<"Mantle \t" <<Scatter_Probability_Mantle[i][0] <<"\t"<<Scatter_Probability_Mantle[i][1]<<"\t"<<Scatter_Probability_Mantle[i][2] <<endl;
	}
}

void Update_PREM(double mX, double sigma0, double velocity)
{
	if(FormFactor == "HelmApproximation")
	{
		g_PREM[0] = 0.0;
		g_PREM[1] = 0.0;
		// 1. MFP Prefactor
		// core
		for(int i = 0; i < 9; i++)
		{
			g_PREM[0] += Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * TotalsigmaSI(mX, sigma0, Elements_Core[i][1], velocity);
		}
		// mantle
		for(int i = 0; i < 14; i++)
		{
			g_PREM[1] += Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * TotalsigmaSI(mX, sigma0, Elements_Mantle[i][1], velocity);
		}
		// 2. Probabilities. The probabilities are velocity depending if form factors are included.
		double total = 0.0;
		// core
		total = 0.0;
		for(int i = 0; i < 9; i++)
		{
			Scatter_Probability_Core[i][2] = Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * TotalsigmaSI(mX, sigma0, Elements_Core[i][1], velocity);
			total += Elements_Core[i][2] / NucleusMass(Elements_Core[i][1]) * TotalsigmaSI(mX, sigma0, Elements_Core[i][1], velocity);
		}
		for(int i = 0; i < 9; i++)
		{
			Scatter_Probability_Core[i][2] /= total;
			// cout <<"Core \t" <<Scatter_Probability_Core[i][0] <<"\t"<<Scatter_Probability_Core[i][1]<<"\t"<<Scatter_Probability_Core[i][2] <<endl;
		}
		// mantle
		total = 0.0;
		for(int i = 0; i < 14; i++)
		{
			Scatter_Probability_Mantle[i][2] = Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * TotalsigmaSI(mX, sigma0, Elements_Mantle[i][1], velocity);
			total += Elements_Mantle[i][2] / NucleusMass(Elements_Mantle[i][1]) * TotalsigmaSI(mX, sigma0, Elements_Mantle[i][1], velocity);
		}
		for(int i = 0; i < 14; i++)
		{
			Scatter_Probability_Mantle[i][2] /= total;
		}
	}
	else
	{
		cout << "Error in Update_PREM." << endl;
	}
}

// Return Prefactor:
double g_Factor(int layer)
{
	if(layer <= 1)
		return g_PREM[0];
	else if(layer <= 9)
		return g_PREM[1];
	else
		return 0.0;
}

// Free Path Vector
// When does a particle leave its current layer?
double tExit(Eigen::Vector3d& position, Eigen::Vector3d& vel)
{
	double r  = position.norm();
	double v  = vel.norm();
	double rv = position.dot(vel);
	int layer = PREM_Layer(r);
	// Check if the particle is well inside the earth.
	if(layer == 10)
	{
		cout << "Error in tExit(): Position argument lies outside the earth." << endl;
		return 0;
	}
	// If not we have to distinguish two cases.
	double R[2];
	// Case A:Most often position lies directly on the border between two layers.
	if(r == PREM_Radii[layer])
	{
		// Position directly at the Earth surface:
		if(layer == 9)
		{
			// Leaving the Earth
			if(rv / r / v >= 0)
				return 1 * km / v;
			// Entering the Earth.
			R[0] = PREM_Radii[9];
			R[1] = PREM_Radii[8];
		}
		// Position directly at core border:
		else if(layer == 0)
		{
			R[0] = PREM_Radii[1];
			R[1] = PREM_Radii[0];
		}
		else
		{
			R[0] = PREM_Radii[layer + 1];
			R[1] = PREM_Radii[layer - 1];
		}
	}
	// Case B: Position is not at the layer border:
	else
	{
		// Core
		if(layer == 0)
			return (-rv + sqrt(pow(rv, 2) - pow(v * r, 2) + pow(v * PREM_Radii[layer], 2))) / pow(v, 2) + 10 * cm / v;
		else
		{
			R[0] = PREM_Radii[layer];
			R[1] = PREM_Radii[layer - 1];
		}
	}
	// With the two radii we can find tExit.
	double texit = 0.0;
	// Outer radius:
	texit = (-rv + sqrt(pow(rv, 2) - pow(v * r, 2) + pow(v * R[0], 2))) / pow(v, 2);
	// Inner Radius
	double radic = pow(rv, 2) - pow(v * r, 2) + pow(v * R[1], 2);
	if(radic >= 0)
	{
		double t1 = (-rv - sqrt(radic)) / pow(v, 2);
		double t2 = (-rv + sqrt(radic)) / pow(v, 2);
		if(t1 > 0 && t1 < texit)
			texit = t1;
		if(t2 > 0 && t2 < texit)
			texit = t2;
	}
	return texit + 10 * cm / v;
}

// Exponent of the scattering probability
double ProbabilityExponent(Eigen::Vector3d& position, Eigen::Vector3d& vel, double L)
{
	double r0		= position.norm();
	double v		= vel.norm();
	double cosalpha = position.dot(vel) / r0 / v;
	int i			= PREM_Layer(r0);
	double Ltilde	= sqrt(L * L + 2 * L * r0 * cosalpha + r0 * r0);
	// Coefficients
	double C1, C2, C3, C4;
	C1			   = L * rEarth * rEarth;
	C3			   = L * (r0 * r0 + r0 * L * cosalpha + L * L / 3);
	double epsilon = 1e-14;
	if(abs(cosalpha + 1) > epsilon)
	{
		C2 = rEarth / 2 * (Ltilde * (L + r0 * cosalpha) - r0 * r0 * cosalpha + (1 - cosalpha * cosalpha) * r0 * r0 * log((L + Ltilde + r0 * cosalpha) / ((1 + cosalpha) * r0)));
		C4 = 1.0 / 8 / rEarth * ((5 - 3 * cosalpha * cosalpha) * (cosalpha * pow(r0, 3) * Ltilde - pow(r0, 4) * cosalpha) + 2 * L * L * Ltilde * (L + 3 * r0 * cosalpha) + L * Ltilde * r0 * r0 * (5 + cosalpha * cosalpha) + 3 * pow(r0, 4) * pow(1 - cosalpha * cosalpha, 2) * log((L + Ltilde + r0 * cosalpha) / ((1 + cosalpha) * r0)));
	}
	else
	{
		C2 = (L * (L - 2 * r0) * sqrt(pow(L - r0, 2)) * rEarth) / (2. * (L - r0));
		C4 = (L * (L - 2 * r0) * sqrt(pow(L - r0, 2)) * (pow(L, 2) - 2 * L * r0 + 2 * pow(r0, 2))) / (4. * (L - r0) * rEarth);
	}
	return g_Factor(i) / rEarth / rEarth * (a[i] * C1 + b[i] * C2 + c[i] * C3 + d[i] * C4) * gram * pow(cm, -3);
}
double ScatterProbability(Eigen::Vector3d& position, Eigen::Vector3d& vel, double L)
{
	return 1 - exp(-ProbabilityExponent(position, vel, L));
}
double NewtonMethod(Eigen::Vector3d& position, Eigen::Vector3d& vel, double xi, double Lambda_Total, double NewtonStart = 500 * km)
{
	double r0		= position.norm();
	double v		= vel.norm();
	double cosalpha = position.dot(vel) / r0 / v;
	int i			= PREM_Layer(r0);
	double delta	= 1e-10;
	double L		= NewtonStart;
	double f, dfdL, Ltilde;
	while(abs(Lambda_Total + ProbabilityExponent(position, vel, L) + log(1 - xi)) > delta)
	{
		Ltilde = sqrt(L * L + 2 * L * r0 * cosalpha + r0 * r0);
		f	   = Lambda_Total + ProbabilityExponent(position, vel, L) + log(1 - xi);
		dfdL   = g_Factor(i) * (a[i] + Ltilde / rEarth * b[i] + pow(Ltilde, 2) / (rEarth * rEarth) * c[i] + pow(Ltilde, 3) / pow(rEarth, 3) * d[i]) * gram / pow(cm, 3);
		L -= f / dfdL;
	}
	return L;
}
// Find the free path vector
Eigen::Vector3d FreePathVector(double mX, double sigma, Eigen::Vector3d& position, Eigen::Vector3d& vel, std::mt19937& PRNG)
{
	if(FormFactor != "None")
		Update_PREM(mX, sigma, vel.norm());
	Eigen::Vector3d r = position;
	double v		  = vel.norm();
	double xi		  = ProbabilitySample(PRNG);
	double logxi	  = -log(1 - xi);
	// Add up terms in the exponent of P(L)
	double Lambda = 0;
	int layer	  = PREM_Layer(r.norm());
	double t, L, Lambda_l, LambdaNew;
	while(layer < 10)
	{
		t		  = tExit(r, vel);
		L		  = t * v;
		Lambda_l  = ProbabilityExponent(r, vel, L);
		LambdaNew = Lambda + Lambda_l;
		// No scattering in this layer
		if(LambdaNew < logxi)
		{
			Lambda = LambdaNew;
			r += t * vel;
			layer = PREM_Layer(r.norm());
		}
		// Scattering in this layer
		else
		{
			// With Newton Method:
			double l = NewtonMethod(r, vel, xi, Lambda, (logxi - Lambda) / Lambda_l * L);
			r += l / v * vel;
			layer = 10;	  // exit while loop
		}
	}
	return (r - position);
}

// Which nucleus does the DM scatter on?
int ScatterNucleus(Eigen::Vector3d& position, std::mt19937& PRNG)
{
	int layer  = PREM_Layer(position.norm());
	double xi  = ProbabilitySample(PRNG);
	double sum = 0.0;
	// The two compositional layers
	// Mantle
	if(layer > 1 && layer < 10)
	{
		for(int i = 0; i < 14; i++)
		{
			sum += Scatter_Probability_Mantle[i][2];
			if(sum > xi)
			{
				return Scatter_Probability_Mantle[i][1];
			}
		}
	}
	// Core
	else if(layer < 2)
	{
		for(int i = 0; i < 9; i++)
		{
			sum += Scatter_Probability_Core[i][2];
			if(sum > xi)
			{
				return Scatter_Probability_Core[i][1];
			}
		}
	}
	// Space
	return 0;
}
