#ifndef __Simulation_Utilities_hpp_
#define __Simulation_Utilities_hpp_

// Headers from libphysica
#include "Linear_Algebra.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"

namespace DaMaSCUS
{
// 1. Event class
struct Event
{
	double time;
	libphysica::Vector position, velocity;

	// Constructors
	Event();
	Event(double t, const libphysica::Vector& pos, const libphysica::Vector& vel);

	double Radius() const;
	double Speed() const;
	double Angular_Momentum() const;

	double Isodetection_Angle(const libphysica::Vector& vel_earth) const;
	int Isodetection_Ring(const libphysica::Vector& vel_earth, unsigned int number_of_rings) const;

	Event In_Units(double unit_distance, double unit_time) const;

	//Overloading the output operator <<
	friend std::ostream& operator<<(std::ostream& output, const Event& event);
};

// 2. Generator of initial conditions
extern Event Initial_Conditions(obscura::DM_Distribution& halo_model);

}	// namespace DaMaSCUS

#endif