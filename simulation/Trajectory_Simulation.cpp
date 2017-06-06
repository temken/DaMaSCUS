#include "Trajectory_Simulation.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "PREM.hpp"
#include "RN_Generators.hpp"
#include "IC_Generator.hpp"


//This function calculates the final velocity of a DM particle of mass mX after it scattered on a nucleus A with initial velocity vini.
Eigen::Vector3d Scatter(Eigen::Vector3d& vini,double mX,double A,std::mt19937& PRNG)
{
	Eigen::Vector3d vfinal;
	double mNucleus=NucleusMass(A);
	if(FormFactor=="HelmApproximation")
	{
		//Unit Vector in direction of velocity
			Eigen::Vector3d ev=vini.normalized();
		//Determination of n, the unit vector pointing into the direction of vfinal, see Landau Lifshitz.
			Eigen::Vector3d ep(ev(2),ev(2),-(ev(0)+ev(1)));
			ep.normalize();
		//scattering angle is not uniform anymore and depends on the velocity via the form factor
			//double chi = ThetaSample(PRNG);
			long double qmax2 = 4.0*Mu(mX,mNucleus)*Mu(mX,mNucleus)*vini.dot(vini);
			if(qmax2==0.0)cout<<"Error in Scatter with HelmApproximation. (qmax2==0)" <<endl;
			double xi=ProbabilitySample(PRNG);
			long double rA2 =pow( (0.3+0.9*pow(mNucleus/GeV,1.0/3.0))*fm ,2.0);
			long double q2;
			//If the argument of the exp becomes too small exp() returns just 1. and no random angle is coming out. In this case we use a taylor expanded result
			if(rA2*qmax2>1.0e-10)q2 = -3.0/rA2*log(1.0-xi*(1.0-exp(-1.0/3.0*rA2*qmax2)));
			else 	q2 = xi*qmax2;
			if(q2==0.0)cout<<"Error in Scatter with HelmApproximation. (q2==0)" <<endl;
			
			double chi = acos(1.0-2.0*q2/qmax2);

			
		//construct vector mit scattering angle chi
			Eigen::Vector3d n;
			if(chi<M_PI/2)
			{
				n=ev+tan(chi)*ep;
			}	
			else
			{
				n=-ev-tan(chi)*ep;
			}
			n.normalize();
			//We chose a particular ep, therefore we finally rotate n around ev by a random angle.
				double alpha=PhiSample(PRNG);
				n=(n.dot(ev))*ev+cos(alpha)*(ev.cross(n)).cross(ev)+sin(alpha)*ev.cross(n);
			//Output
				vfinal=mNucleus/(mX+mNucleus)*vini.norm()*n+mX/(mX+mNucleus)*vini;	
	}
	else
	{
		//Unit vector with uniformly distributed random direction:
			Eigen::Vector3d n=SphericalCoordinates(1,ThetaSample(PRNG),PhiSample(PRNG));
		//Output
			vfinal=mNucleus/(mX+mNucleus)*vini.norm()*n+mX/(mX+mNucleus)*vini;
	}

	return vfinal;
}

double VelocityCut(double mX,double Ethreshold,double A)
{
	return (mX+NucleusMass(A))/mX*sqrt(Ethreshold/2/NucleusMass(A));
}

//This function performs the trajectory simulation for a single DM particle and saves the events of scattering in one list.
Trajectory ParticleTrack(double mX,double sigman0,Event IniCondi, double vcut, std::mt19937& PRNG,double rFinal)
{
	//We start outside the earth
		bool InsideEarth=false;
	//Vector for the output:
		vector<Event> EventList;
	//Initial Conditions
		EventList.push_back(IniCondi);
		double t=IniCondi.Time();
		Eigen::Vector3d x=IniCondi.Position();
		Eigen::Vector3d v=IniCondi.Velocity();
	//First we check, if the particle even passes the surface of the earth. This is usually never called upon, since the IC generator doesn't generate particles pointing away from the earth.
		double alphaCone;
		if(abs(x.norm()-rEarth)<1.0*meter) alphaCone=M_PI/2.0;//this is necessary in the case that the particles start on the earth surface. then rearth/x.norm can be >1 for numerical reasons. And this fucks shit up. sin(1.00000000000001)=nan.
		else alphaCone=asin(rEarth/x.norm());
		if (acos(-x.normalized().dot(v.normalized()))<alphaCone)
		{
			double TimeOfEntry =(-x.dot(v)-sqrt(pow(x.dot(v),2)-v.dot(v)*(x.dot(x)-pow(rEarth,2))))/(v.dot(v));
			t+=TimeOfEntry;
			x+=v*TimeOfEntry;
			EventList.push_back(Event(t,x,v));
			InsideEarth=true;
		}
		else //if (acos(-x.normalized().dot(v.normalized()))>asin(rEarth/x.norm())) //If not we just append an additional point and do not enter the loop.
		{
			t+=4*rEarth/v.norm();
			x+=4*rEarth*v.normalized();
			EventList.push_back(Event(t,x,v));
		}
	//Simulation Underground:
		while(InsideEarth)
		{
			//New Position
				Eigen::Vector3d xnew=x+FreePathVector(mX,sigman0,x,v,PRNG);
			//Still inside the planet?
				if (xnew.norm()>rEarth)
				{
					InsideEarth=false;
					double TimeOfExit =(-x.dot(v)+sqrt(pow(x.dot(v),2)-v.dot(v)*(x.dot(x)-pow(rEarth,2))))/(v.dot(v));
					//Save Point of exit
					t+=TimeOfExit;
					x+=TimeOfExit*v;
					EventList.push_back(Event(t,x,v));
					//Save final point
					t+=rFinal/v.norm();
					x+=rFinal*v.normalized();
					EventList.push_back(Event(t,x,v));
				}
			//Still in Earth! ->Scattering occurs.
				else
				{	
					t+=(xnew-x).norm()/v.norm();
					x=xnew;
					v=Scatter(v,mX,ScatterNucleus(x,PRNG),PRNG);
					EventList.push_back(Event(t,x,v));
					//are we still above our velocity cut? if not, we abort the loop.
					if (v.norm()<vcut)
					{
						InsideEarth=false;
					}
				}
		}
	return Trajectory(EventList);
}






