#ifndef __PREM_hpp_
#define __PREM_hpp_

#include <Eigen/Geometry>
#include <random>

#include "Physical_Parameters.hpp"

//PREM Functions 
	extern int PREM_Layer(double r);


//Mean Free Path PreFactor
	extern double g_PREM[2];
	extern double g_Factor(int layer);
	extern void Initialize_PREM(double mChi,double sigma0);
	extern void Update_PREM(double mChi,double sigma0,double velocity);

//tExit Function: Calculates when the particle exits the current layer
	extern double tExit(Eigen::Vector3d& position, Eigen::Vector3d& vel);


//Free Path Vector
	extern Eigen::Vector3d FreePathVector(double mX,double sigma,Eigen::Vector3d& position,Eigen::Vector3d& vel,std::mt19937& PRNG);
//Which nucleus does the DM scatter on?
	extern int ScatterNucleus(Eigen::Vector3d& position,std::mt19937& PRNG);

#endif