#include "Physical_Parameters.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

#include "General_Utilities.hpp"
// Units:
// Energy
const double GeV = 1.0;
const double eV	 = 1.0E-9 * GeV;
const double keV = 1.0E-6 * GeV;
const double MeV = 1.0E-3 * GeV;
const double TeV = 1.0E3 * GeV;
// Mass
const double gram = 5.617977528089887E23 * GeV;
const double kg	  = 1E3 * gram;
// Length
const double cm			 = 5.067E13 / GeV;
const double mm			 = 0.1 * cm;
const double meter		 = 100 * cm;
const double km			 = 1000 * meter;
const double fm			 = 1E-15 * meter;
const double pb			 = 1E-36 * pow(cm, 2);
const double parsec		 = 3.0857E16 * meter;
const double kpc		 = 1E3 * parsec;
const double Mpc		 = 1E6 * parsec;
const double Bohr_Radius = 5.291772083e-11 * meter;
// Time
const double sec	= 299792458 * meter;
const double minute = 60 * sec;
const double hour	= 60 * minute;
const double day	= 24 * hour;
const double year	= 365 * day;
// Temperature
const double Kelvin = 8.62E-14 * GeV;
// Others
const double erg = gram * pow(cm / sec, 2);

// Specific Parameters:
// Masses
const double mPlanck   = 1.2209E19 * GeV;
const double GNewton   = pow(mPlanck, -2);
const double mProton   = 0.938 * GeV;
const double mElectron = 0.511 * MeV;
const double mNucleon  = 0.932 * GeV;
// Geographic Parameters
const double mEarth	  = 5.972E24 * kg;
const double rEarth	  = 6371 * km;
const double rhoEarth = 5.51 * gram * pow(cm, -3);
// Solar Parameters
const double mSun	= 1.989E30 * kg;
const double rSun	= 6.957E8 * meter;
const double rhoSun = 151 * gram * pow(cm, -3);
// Dark Matter Halo Parameters
const double v0	  = 220 * km / sec;
const double vesc = 544 * km / sec;

const double Nesc = M_PI * v0 * v0 * (sqrt(M_PI) * v0 * erf(vesc / v0) - 2 * vesc * exp(-vesc * vesc / v0 / v0));
const double ymax = 4 * M_PI / Nesc * v0 * v0 * exp(-1);
// degree to rad
const double deg = M_PI / 180.0;
// Unit Conversion
double InUnits(double quantity, double dimension)
{
	return quantity / dimension;
}

// Average inverse velocity
double AverageInvVelocity(double vE)
{
	if(vE == 0.0)
	{
		return 2 * M_PI * v0 * v0 / Nesc * (1 - exp(-vesc * vesc / v0 / v0));
	}
	return pow(M_PI, 3.0 / 2.0) * pow(v0, 3) / Nesc / vE * (erf(vE / v0) - 2 / sqrt(M_PI) * vE / v0 * exp(-vesc * vesc / v0 / v0));
}

// Physical Functions

// Nucleus Mass
double NucleusMass(double A)
{
	return A * mNucleon;
}

// Reduced mass
double Mu(double m1, double m2)
{
	return m1 * m2 / (m1 + m2);
}

// cross section for wimp nucleus scattering with zero momentum transfer
double sigmaSI(double mX, double sigman0, double A)
{
	double Z = (A == 1) ? 1 : 0.5 * A;
	return sigman0 * pow(Mu(mX, NucleusMass(A)), 2) / pow(Mu(mX, mNucleon), 2) * pow(Z, 2);
}

double TotalsigmaSI(double mX, double sigman0, double A, double vX)
{
	if(FormFactor != "None" && vX == 0)
		cout << "Error in TotalsigmaSI." << endl;
	// form factor
	long double ff = 1.0;
	if(FormFactor == "HelmApproximation")
	{
		double mA		  = NucleusMass(A);
		double RA2		  = pow((0.3 + 0.9 * pow(mA / GeV, 1.0 / 3.0)) * fm, 2.0);
		long double alpha = 3.0 / 4.0 / RA2 / pow(Mu(mX, mA), 2.0) / vX / vX;
		// Taylor expansion of alpha*(1.0-exp(-1.0/alpha)) around alpha=inf
		if(alpha < 1.0e6)
			ff = alpha * (1.0 - exp(-1.0 / alpha));
		else
			ff = 1.0 - 2.0 / alpha;
		if(ff == 0.0)
			cout << "Error in TotalsigmaSI. ff=0" << endl;
	}
	else if(FormFactor == "ChargeScreening" || FormFactor == "LightMediator")
	{
		int Z			= (A == 1) ? 1 : A / 2;
		double a		= Thomas_Fermi_Radius(Z);
		double mNucleus = NucleusMass(A);
		double q2max	= 4.0 * std::pow(Mu(mX, mNucleus), 2.0) * vX * vX;
		double x		= a * a * q2max;
		ff				= (FormFactor == "ChargeScreening") ? 1.0 + 1.0 / (1.0 + x) - 2.0 / x * log(1.0 + x) : x * x / (1.0 + x);
	}
	return sigmaSI(mX, sigman0, A) * ff;
}

// Form Factor
double FF_HelmApproximation(double qSquared, double A)
{
	double RA = (0.3 + 0.9 * pow(NucleusMass(A) / GeV, 1.0 / 3.0)) * fm;
	return exp(-qSquared * pow(RA, 2.0) / 3.0);
}
// Integrated from 0 to qmax
double FF_HelmApproximation_Integrated(double mX, double vDM, double A)
{
	double mA	 = NucleusMass(A);
	double qMax2 = 4.0 * pow(Mu(mX, mA) * vDM, 2.0);
	double RA2	 = pow((0.3 + 0.9 * pow(mA / GeV, 1.0 / 3.0)) * fm, 2.0);
	return 3.0 / RA2 * (1 - exp(-qMax2 * RA2 / 3.0));
}

double Thomas_Fermi_Radius(int Z)
{
	return 0.25 * std::pow(9.0 * M_PI * M_PI / 2.0 / Z, 1.0 / 3.0) * Bohr_Radius;
}

// Coordinate Systems
Eigen::Vector3d SphericalCoordinates(double r, double theta, double phi)
{
	Eigen::Vector3d v(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
	return v;
}

Eigen::Vector3d Equat2Gal(Eigen::Vector3d& v, double T)
{
	double arcsecond = deg / 3600;
	double zetaA	 = 2306.083227 * arcsecond * T + 0.298850 * arcsecond * T * T;
	double zA		 = 2306.077181 * arcsecond * T + 1.092735 * arcsecond * T * T;
	double thetaA	 = 2004.191903 * arcsecond * T - 0.429493 * arcsecond * T * T;
	double lCP		 = 122.932 * deg;
	double alphaGP	 = 192.85948 * deg;
	double deltaGP	 = 27.12825 * deg;
	Eigen::Matrix3d P;
	Eigen::Matrix3d M;
	P << cos(zetaA) * cos(thetaA) * cos(zA) - sin(zetaA) * sin(zA), -sin(zetaA) * cos(thetaA) * cos(zA) - cos(zetaA) * sin(zA), -sin(thetaA) * cos(zA),
		cos(zetaA) * cos(thetaA) * sin(zA) + sin(zetaA) * cos(zA), -sin(zetaA) * cos(thetaA) * sin(zA) + cos(zetaA) * cos(zA), -sin(thetaA) * sin(zA),
		cos(zetaA) * sin(thetaA), -sin(zetaA) * sin(thetaA), cos(thetaA);
	M << -sin(lCP) * sin(alphaGP) - cos(lCP) * cos(alphaGP) * sin(deltaGP), sin(lCP) * cos(alphaGP) - cos(lCP) * sin(alphaGP) * sin(deltaGP), cos(lCP) * cos(deltaGP),
		cos(lCP) * sin(alphaGP) - sin(lCP) * cos(alphaGP) * sin(deltaGP), -cos(lCP) * cos(alphaGP) - sin(lCP) * sin(alphaGP) * sin(deltaGP), sin(lCP) * cos(deltaGP),
		cos(alphaGP) * cos(deltaGP), sin(alphaGP) * cos(deltaGP), sin(deltaGP);
	return M * P.inverse() * v;
}
Eigen::Vector3d GeoEcl2Gal(Eigen::Vector3d& v, double T)
{
	double arcsecond = deg / 3600;
	double zetaA	 = 2306.083227 * arcsecond * T + 0.298850 * arcsecond * T * T;
	double zA		 = 2306.077181 * arcsecond * T + 1.092735 * arcsecond * T * T;
	double thetaA	 = 2004.191903 * arcsecond * T - 0.429493 * arcsecond * T * T;
	double lCP		 = 122.932 * deg;
	double alphaGP	 = 192.85948 * deg;
	double deltaGP	 = 27.12825 * deg;
	Eigen::Matrix3d P;
	Eigen::Matrix3d M;
	Eigen::Matrix3d R;
	P << cos(zetaA) * cos(thetaA) * cos(zA) - sin(zetaA) * sin(zA), -sin(zetaA) * cos(thetaA) * cos(zA) - cos(zetaA) * sin(zA), -sin(thetaA) * cos(zA),
		cos(zetaA) * cos(thetaA) * sin(zA) + sin(zetaA) * cos(zA), -sin(zetaA) * cos(thetaA) * sin(zA) + cos(zetaA) * cos(zA), -sin(thetaA) * sin(zA),
		cos(zetaA) * sin(thetaA), -sin(zetaA) * sin(thetaA), cos(thetaA);
	M << -sin(lCP) * sin(alphaGP) - cos(lCP) * cos(alphaGP) * sin(deltaGP), sin(lCP) * cos(alphaGP) - cos(lCP) * sin(alphaGP) * sin(deltaGP), cos(lCP) * cos(deltaGP),
		cos(lCP) * sin(alphaGP) - sin(lCP) * cos(alphaGP) * sin(deltaGP), -cos(lCP) * cos(alphaGP) - sin(lCP) * sin(alphaGP) * sin(deltaGP), sin(lCP) * cos(deltaGP),
		cos(alphaGP) * cos(deltaGP), sin(alphaGP) * cos(deltaGP), sin(deltaGP);
	double eps = 23.4393 * deg - 0.0130 * deg * T;
	R << 1, 0, 0,
		0, cos(eps), -sin(eps),
		0, sin(eps), cos(eps);

	return M * P.inverse() * R * v;
}
Eigen::Vector3d HelEcl2Gal(Eigen::Vector3d& v, double T)
{
	double arcsecond = deg / 3600;
	double zetaA	 = 2306.083227 * arcsecond * T + 0.298850 * arcsecond * T * T;
	double zA		 = 2306.077181 * arcsecond * T + 1.092735 * arcsecond * T * T;
	double thetaA	 = 2004.191903 * arcsecond * T - 0.429493 * arcsecond * T * T;
	double lCP		 = 122.932 * deg;
	double alphaGP	 = 192.85948 * deg;
	double deltaGP	 = 27.12825 * deg;
	Eigen::Matrix3d P;
	Eigen::Matrix3d M;
	Eigen::Matrix3d R;
	P << cos(zetaA) * cos(thetaA) * cos(zA) - sin(zetaA) * sin(zA), -sin(zetaA) * cos(thetaA) * cos(zA) - cos(zetaA) * sin(zA), -sin(thetaA) * cos(zA),
		cos(zetaA) * cos(thetaA) * sin(zA) + sin(zetaA) * cos(zA), -sin(zetaA) * cos(thetaA) * sin(zA) + cos(zetaA) * cos(zA), -sin(thetaA) * sin(zA),
		cos(zetaA) * sin(thetaA), -sin(zetaA) * sin(thetaA), cos(thetaA);
	M << -sin(lCP) * sin(alphaGP) - cos(lCP) * cos(alphaGP) * sin(deltaGP), sin(lCP) * cos(alphaGP) - cos(lCP) * sin(alphaGP) * sin(deltaGP), cos(lCP) * cos(deltaGP),
		cos(lCP) * sin(alphaGP) - sin(lCP) * cos(alphaGP) * sin(deltaGP), -cos(lCP) * cos(alphaGP) - sin(lCP) * sin(alphaGP) * sin(deltaGP), sin(lCP) * cos(deltaGP),
		cos(alphaGP) * cos(deltaGP), sin(alphaGP) * cos(deltaGP), sin(deltaGP);
	double eps = 23.4393 * deg - 0.0130 * deg * T;
	R << 1, 0, 0,
		0, cos(eps), -sin(eps),
		0, sin(eps), cos(eps);

	return -M * P.inverse() * R * v;
}

Eigen::Vector3d Gal2Equat(Eigen::Vector3d& v, double T)
{
	double arcsecond = deg / 3600;
	double zetaA	 = 2306.083227 * arcsecond * T + 0.298850 * arcsecond * T * T;
	double zA		 = 2306.077181 * arcsecond * T + 1.092735 * arcsecond * T * T;
	double thetaA	 = 2004.191903 * arcsecond * T - 0.429493 * arcsecond * T * T;
	double lCP		 = 122.932 * deg;
	double alphaGP	 = 192.85948 * deg;
	double deltaGP	 = 27.12825 * deg;
	Eigen::Matrix3d P;
	Eigen::Matrix3d M;
	P << cos(zetaA) * cos(thetaA) * cos(zA) - sin(zetaA) * sin(zA), -sin(zetaA) * cos(thetaA) * cos(zA) - cos(zetaA) * sin(zA), -sin(thetaA) * cos(zA),
		cos(zetaA) * cos(thetaA) * sin(zA) + sin(zetaA) * cos(zA), -sin(zetaA) * cos(thetaA) * sin(zA) + cos(zetaA) * cos(zA), -sin(thetaA) * sin(zA),
		cos(zetaA) * sin(thetaA), -sin(zetaA) * sin(thetaA), cos(thetaA);
	M << -sin(lCP) * sin(alphaGP) - cos(lCP) * cos(alphaGP) * sin(deltaGP), sin(lCP) * cos(alphaGP) - cos(lCP) * sin(alphaGP) * sin(deltaGP), cos(lCP) * cos(deltaGP),
		cos(lCP) * sin(alphaGP) - sin(lCP) * cos(alphaGP) * sin(deltaGP), -cos(lCP) * cos(alphaGP) - sin(lCP) * sin(alphaGP) * sin(deltaGP), sin(lCP) * cos(deltaGP),
		cos(alphaGP) * cos(deltaGP), sin(alphaGP) * cos(deltaGP), sin(deltaGP);
	return P * M.inverse() * v;
}

// Time-related Functions
// fractional days since J2000.0
double FractionalDays(int date[3], int time[3])
{
	double n = 0.0;
	if(date[1] == 1 || date[1] == 2)
	{
		n += floor(365.25 * (date[2] - 1)) + floor(30.61 * (date[1] + 13));
	}
	else
	{
		n += floor(365.25 * date[2]) + floor(30.61 * (date[1] + 1));
	}
	n += date[0] - 730563.5 + time[0] / 24. + time[1] / 1440. + time[2] / 86400.;
	return n;
}
// Local Apparent Sidereal Time for a given Longitude and date. Returns GAST without longitude input
double LASTinSeconds(double nJ2000, double longitude)
{
	double T = nJ2000 / 36525.;
	// GMST:
	double t		= 86400. * (0.779057273264 + 0.00273781191135448 * nJ2000 + fmod(nJ2000, 1)) + 0.00096707 + 307.47710227 * T + 0.092772113 * pow(T, 2);
	double epsilonA = 23.439279444 * deg - 0.01301021361 * deg * T;
	double Omega	= 125.04455501 * deg - 0.05295376 * deg * nJ2000;
	double L		= 280.47 * deg - 0.98565 * deg * nJ2000;
	double DeltaPsi = -1.1484 * sin(Omega) - 0.0864 * cos(2 * L);
	// //GAST:
	t += DeltaPsi * cos(epsilonA) + 0.000176 * sin(Omega) + 0.000004 * sin(2 * Omega);
	// LAST:
	t += longitude / 2 / M_PI * 86400.;
	return fmod(t, 86400);
}

// Earth's Velocity in the galactic frame
Eigen::Vector3d EarthVelocity(double n)
{
	// 1. Sun's rotation around galactic center:
	Eigen::Vector3d vGal(0, 220 * km / sec, 0);
	// 2. Sun's peculiar motion:
	Eigen::Vector3d vPec(11.1 * km / sec, 12.2 * km / sec, 7.3 * km / sec);
	// 3. Earth's rotation around galactic center:
	double e	 = 0.01671;	  // Ellipticity of earth's orbit
	double L	 = fmod(280.46 * deg + n * 0.9856474 * deg, 2 * M_PI);
	double omega = fmod(282.932 * deg + n * 0.0000471 * deg, 2 * M_PI);
	double T	 = n / 36525.0;
	Eigen::Vector3d ex(1, 0, 0);
	Eigen::Vector3d ey(0, 1, 0);
	// Basis vectors in heliocentric ecliptic coordinate system
	Eigen::Vector3d exEcl = HelEcl2Gal(ex, T);	 //(0.054876-0.024232*T,-0.494109-0.002689*T,0.867666+1.546e-6*T);
	Eigen::Vector3d eyEcl = HelEcl2Gal(ey, T);	 //(0.993824+0.001316*T,0.110992-0.011851*T,0.000352+0.021267*T);
	double ve			  = 29.79 * km / sec;
	Eigen::Vector3d uE	  = -ve * (sin(L) + e * sin(2 * L - omega)) * exEcl + ve * (cos(L) + e * cos(2 * L - omega)) * eyEcl;
	// Return the sum of all components
	return vGal + vPec + uE;
}
