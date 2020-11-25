#ifndef __Simulation_Trajectory_hpp_
#define __Simulation_Trajectory_hpp_

#include <vector>

// Headers from libphysica
#include "Natural_Units.hpp"

// Headers from obscura
#include "DM_Particle.hpp"

#include "Earth_Model.hpp"

namespace DaMaSCUS
{

// 1. Result of one trajectory

class Trajectory
{
  private:
	std::vector<Event> events;

  public:
	Trajectory(std::vector<Event>& event_list);

	unsigned int Number_of_Scatterings() const;
	double Deepest_Depth() const;

	Event Interpolate_Trajectory(double t) const;

	std::vector<Event> Points_of_Depth_Crossing(double underground_depth) const;

	void Save_to_File(std::string& file_path) const;
	void Print_Summary(int mpi_rank = 0) const;
};

// 2. Simulator
Trajectory Simulate_Trajectory(Event initial_conditions, Earth_Model& earth_model, obscura::DM_Particle& DM, std::mt19937& PRNG, double minimal_speed = libphysica::natural_units::cm / libphysica::natural_units::sec);

}	// namespace DaMaSCUS

#endif