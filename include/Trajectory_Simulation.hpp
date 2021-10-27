#ifndef __Trajectory_Simulation_hpp_
#define __Trajectory_Simulation_hpp_

#include <Eigen/Geometry>

#include <random>

#include "General_Utilities.hpp"
#include "Physical_Parameters.hpp"
#include "Trajectory_Class.hpp"

extern Eigen::Vector3d Scatter(Eigen::Vector3d& vini, double mX, double A, std::mt19937& PRNG);

extern double VelocityCut(double mX, double Ethreshold = 1 * eV, double A = 131);

extern Trajectory ParticleTrack(double mX, double sigman0, Event IniCondi, double vcut, std::mt19937& PRNG, double rFinal = 1 * meter);

#endif