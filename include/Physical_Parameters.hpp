#ifndef __Physical_Parameters_hpp_
#define __Physical_Parameters_hpp_

#include <Eigen/Geometry>
#include <vector>

using namespace std;

// Units:
// Energy
extern const double GeV, eV, keV, MeV, TeV;
// Mass
extern const double gram, kg;
// Length
extern const double cm, mm, meter, km, fm, pb, parsec, kpc, Mpc, Bohr_Radius;
// Time
extern const double sec, minute, hour, day, year;
// Temperature
extern const double Kelvin;
// Others
extern const double erg;

// Specific Parameters:
// Masses
extern const double mPlanck, GNewton, mProton, mElectron, mNucleon;
// Geographic Parameters
extern const double mEarth, rEarth, rhoEarth;
// Solar Parameters
extern const double mSun, rSun, rhoSun;
// Dark Matter Halo Parameters
extern const double v0, vesc;
extern const double Nesc;
extern const double ymax;
// degree to rad
extern const double deg;

// Unit Conversion
extern double InUnits(double quantity, double dimension);
// Reduced Mass
extern double Mu(double m1, double m2);

// Average inverse velocity
extern double AverageInvVelocity(double vE);

// Nucleus Mass
extern double NucleusMass(double A);

// Wimp Nucleus Cross-section with zero momentum transfer:
extern double sigmaSI(double mX, double sigman0, double A);
// Total Wimp Nucleus Cross-section:
extern double TotalsigmaSI(double mX, double sigman0, double A, double vX = 0.0);
// Helm Form Factor (approximation)
extern double FF_HelmApproximation(double qSquared, double A);
extern double FF_HelmApproximation_Integrated(double mX, double vDM, double A);

// Charge screening
extern double Thomas_Fermi_Radius(int Z);

// Coordinate System Change
extern Eigen::Vector3d SphericalCoordinates(double r, double theta, double phi);
extern Eigen::Vector3d Equat2Gal(Eigen::Vector3d& v, double T = 0.0);
extern Eigen::Vector3d GeoEcl2Gal(Eigen::Vector3d& v, double T = 0.0);
extern Eigen::Vector3d HelEcl2Gal(Eigen::Vector3d& v, double T = 0.0);
extern Eigen::Vector3d Gal2Equat(Eigen::Vector3d& v, double T = 0.0);

// fractional days
extern double FractionalDays(int date[3], int time[3]);
extern double LASTinSeconds(double nJ2000, double longitude = 0);
// Earth's Velocity in the galactic frame
extern Eigen::Vector3d EarthVelocity(double n = 0.0);

#endif
