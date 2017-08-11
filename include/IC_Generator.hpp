#ifndef __IC_Generator_hpp_
#define __IC_Generator_hpp_

#include <Eigen/Geometry>
#include <cmath>
#include <random>

#include "Trajectory_Class.hpp"
#include "Physical_Parameters.hpp"


//Initial Condition Generator
extern Event InitialCondition(std::mt19937& PRNG);
extern Event InitialCondition(double t,Eigen::Vector3d &vearth,std::mt19937& PRNG,double R = 1.1*rEarth);


#endif