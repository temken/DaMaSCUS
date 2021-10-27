#include "IC_Generator.hpp"

#include <fstream>
#include <iostream>

#include "General_Utilities.hpp"
#include "RN_Generators.hpp"

// Two definitions for the default values
Event InitialCondition(std::mt19937& PRNG)
{
	return InitialCondition(0, vEarth, PRNG);
}
Event InitialCondition(double t, Eigen::Vector3d& vearth, std::mt19937& PRNG, double R)
{
	Eigen::Vector3d IniPosi, IniVeli;
	// Initial Velocity
	double vnorm = MaxwellSample(PRNG);
	IniVeli		 = SphericalCoordinates(vnorm, ThetaSample(PRNG), PhiSample(PRNG));
	// Boost into the geocentric frame:
	IniVeli -= vearth;
	// Random Point in a circle at distance R
	Eigen::Vector3d ez = -IniVeli.normalized();
	Eigen::Vector3d ex(0, ez[2] / sqrt(ez[1] * ez[1] + ez[2] * ez[2]), -ez[1] / sqrt(ez[1] * ez[1] + ez[2] * ez[2]));
	Eigen::Vector3d ey = ez.cross(ex);
	double phi		   = PhiSample(PRNG);
	double xi		   = ProbabilitySample(PRNG);
	IniPosi			   = R * ez + pow(xi, 1.0 / 2.0) * rEarth * (cos(phi) * ex + sin(phi) * ey);
	// Projection onto a sphere around earth (optional, but I have to do the projection also in time!)
	//  double tProject = 1/pow(IniVeli.norm(),2)*(-IniPosi.dot(IniVeli)-sqrt(pow(IniPosi.dot(IniVeli),2)+pow(IniVeli.norm(),2)*(R*R-pow(IniPosi.norm(),2))));
	//  IniPosi+=tProject*IniVeli;
	//  t+=tProject;

	if(acos(-IniPosi.normalized().dot(IniVeli.normalized())) > asin(rEarth / IniPosi.norm()))
	{
		cout << "Error in InitialCondition(): Send in particle misses the Earth." << endl;
	}
	// Return the result

	return Event(t, IniPosi, IniVeli);
}
