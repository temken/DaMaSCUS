#include "Trajectory_Simulation.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "IC_Generator.hpp"
#include "PREM.hpp"
#include "RN_Generators.hpp"

// This function calculates the final velocity of a DM particle of mass mX after it scattered on a nucleus A with initial velocity vini.

double Sample_Scattering_Angle(double mDM, double mNucleus, Eigen::Vector3d& vini, std::string form_factor, std::mt19937& PRNG)
{
	double cos_chi = 0.0;
	if(form_factor == "None")
		cos_chi = 2.0 * ProbabilitySample(PRNG) - 1.0;
	else if(form_factor == "HelmApproximation")
	{
		double xi	 = ProbabilitySample(PRNG);
		double rA2	 = pow((0.3 + 0.9 * pow(mNucleus / GeV, 1.0 / 3.0)) * fm, 2.0);
		double qmax2 = 4.0 * Mu(mDM, mNucleus) * Mu(mDM, mNucleus) * vini.dot(vini);
		if(qmax2 == 0.0)
			cout << "Error in Scatter with HelmApproximation. (qmax2==0)" << endl;
		double q2;
		// If the argument of the exp becomes too small exp() returns just 1. and no random angle is coming out. In this case we use a taylor expanded result
		if(rA2 * qmax2 > 1.0e-10)
			q2 = -3.0 / rA2 * log(1.0 - xi * (1.0 - exp(-1.0 / 3.0 * rA2 * qmax2)));
		else
			q2 = xi * qmax2;
		if(q2 == 0.0)
			cout << "Error in Scatter with HelmApproximation. (q2==0)" << endl;
		cos_chi = 1.0 - 2.0 * q2 / qmax2;
	}
	else if(form_factor == "ChargeScreening")
	{
	}
	else if(form_factor == "LightMediator")
	{
	}
	else
	{
		cout << "Error in Sample_Scattering_Angle(): Form factor " << form_factor << " not recognized." << endl;
		std::exit(EXIT_FAILURE);
	}
	std::cout << "Scattering angle: " << cos_chi << std::endl;
	return cos_chi;
}

Eigen::Vector3d Scatter(Eigen::Vector3d& vini, double mX, double A, std::mt19937& PRNG)
{
	double mNucleus = NucleusMass(A);

	double cos_scattering_angle = Sample_Scattering_Angle(mX, mNucleus, vini, FormFactor, PRNG);

	// Construction of n, the unit vector pointing into the direction of vfinal.
	double sin_scattering_angle = sqrt(1.0 - cos_scattering_angle * cos_scattering_angle);
	double phi					= PhiSample(PRNG);
	double cos_phi				= cos(phi);
	double sin_phi				= sin(phi);

	Eigen::Vector3d ev = vini.normalized();
	double aux		   = sqrt(1.0 - pow(ev[2], 2.0));
	Eigen::Vector3d n({cos_scattering_angle * ev[0] + sin_scattering_angle / aux * (ev[0] * ev[2] * cos_phi - ev[1] * sin_phi),
					   cos_scattering_angle * ev[1] + sin_scattering_angle / aux * (ev[1] * ev[2] * cos_phi + ev[0] * sin_phi),
					   cos_scattering_angle * ev[2] - aux * cos_phi * sin_scattering_angle});

	Eigen::Vector3d vfinal = mNucleus / (mX + mNucleus) * vini.norm() * n + mX / (mX + mNucleus) * vini;
	return vfinal;
}

double VelocityCut(double mX, double Ethreshold, double A)
{
	return (mX + NucleusMass(A)) / mX * sqrt(Ethreshold / 2 / NucleusMass(A));
}

// This function performs the trajectory simulation for a single DM particle and saves the events of scattering in one list.
Trajectory ParticleTrack(double mX, double sigman0, Event IniCondi, double vcut, std::mt19937& PRNG, double rFinal)
{
	// We start outside the earth
	bool InsideEarth = false;
	// Vector for the output:
	vector<Event> EventList;
	// Initial Conditions
	EventList.push_back(IniCondi);
	double t		  = IniCondi.Time();
	Eigen::Vector3d x = IniCondi.Position();
	Eigen::Vector3d v = IniCondi.Velocity();
	// First we check, if the particle even passes the surface of the earth. This is usually never called upon, since the IC generator doesn't generate particles pointing away from the earth.
	double alphaCone;
	if(abs(x.norm() - rEarth) < 1.0 * meter)
		alphaCone = M_PI / 2.0;	  // this is necessary in the case that the particles start on the earth surface. then rearth/x.norm can be >1 for numerical reasons.-> sin(1.00000000000001)=nan.
	else
		alphaCone = asin(rEarth / x.norm());
	if(acos(-x.normalized().dot(v.normalized())) < alphaCone)
	{
		double TimeOfEntry = (-x.dot(v) - sqrt(pow(x.dot(v), 2) - v.dot(v) * (x.dot(x) - pow(rEarth, 2)))) / (v.dot(v));
		t += TimeOfEntry;
		x += v * TimeOfEntry;
		// FreePathVector() might interpret sitting on top the surface as being outside the Earth due to numerical imprecision.
		// This would result in a particle track without a single scattering, which is why we move the particle 1mm underground at the start.
		if(x.norm() > rEarth)
			x += mm * v.normalized();
		// Save point of entry and continue
		EventList.push_back(Event(t, x, v));
		InsideEarth = true;
	}
	// If the particle doesnt pass the earth surface for whatever reason, we just append an additional point and do not enter the loop.
	else
	{
		t += 4 * rEarth / v.norm();
		x += 4 * rEarth * v.normalized();
		EventList.push_back(Event(t, x, v));
	}
	// Simulation Underground:
	while(InsideEarth)
	{
		// New Position
		Eigen::Vector3d xnew = x + FreePathVector(mX, sigman0, x, v, PRNG);
		// Still inside the planet?
		if(xnew.norm() > rEarth)
		{
			InsideEarth		  = false;
			double TimeOfExit = (-x.dot(v) + sqrt(pow(x.dot(v), 2) - v.dot(v) * (x.dot(x) - pow(rEarth, 2)))) / (v.dot(v));
			// Save Point of exit
			t += TimeOfExit;
			x += TimeOfExit * v;
			EventList.push_back(Event(t, x, v));
			// Save final point
			t += rFinal / v.norm();
			x += rFinal * v.normalized();
			EventList.push_back(Event(t, x, v));
		}
		// Still in Earth! ->Scattering occurs.
		else
		{
			t += (xnew - x).norm() / v.norm();
			x = xnew;
			v = Scatter(v, mX, ScatterNucleus(x, PRNG), PRNG);
			EventList.push_back(Event(t, x, v));
			// are we still above our velocity cut? if not, we abort the loop.
			if(v.norm() < vcut)
			{
				InsideEarth = false;
			}
		}
	}
	return Trajectory(EventList);
}
