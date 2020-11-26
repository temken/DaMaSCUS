#ifndef __Data_Generation_hpp_
#define __Data_Generation_hpp_

#include <random>
#include <string>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"
#include "DM_Particle.hpp"

#include "Earth_Model.hpp"

namespace DaMaSCUS
{
class Simulation_Data
{
  private:
	// Configuration
	unsigned int min_sample_size_above_threshold;
	double minimum_speed_threshold;
	unsigned int isodetection_rings;
	double initial_and_final_radius = 1.1 * libphysica::natural_units::rEarth;
	double underground_depth;

	std::mt19937 PRNG;

	// Results
	unsigned long int number_of_trajectories;
	unsigned long int number_of_free_particles;
	unsigned long int number_of_stuck_particles;
	double average_number_of_scatterings;
	double computing_time;

	std::vector<unsigned long int> number_of_data_points;

	// MPI
	int mpi_rank, mpi_processes;
	void Perform_MPI_Reductions();

  public:
	std::vector<std::vector<libphysica::DataPoint>> data;

	Simulation_Data(unsigned int sample_size, double depth, double v_min = 0.0, unsigned int iso_rings = 1);

	void Generate_Data(obscura::DM_Particle& DM, Earth_Model& earth_model, obscura::DM_Distribution& halo_model);

	double Average_Scatterings() const;
	double Average_Weight(unsigned int iso_ring) const;

	double Free_Ratio() const;
	double Stuck_Ratio() const;

	double Lowest_Speed(unsigned int iso_ring = 0) const;
	double Highest_Speed(unsigned int iso_ring = 0) const;

	std::string Summary();
};
}	// namespace DaMaSCUS

#endif