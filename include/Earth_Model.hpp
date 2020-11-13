#ifndef __Earth_Model_hpp_
#define __Earth_Model_hpp_

#include <random>

// Headers from obscura
#include "DM_Particle.hpp"
#include "Target_Nucleus.hpp"

#include "Simulation_Utilities.hpp"

namespace DaMaSCUS
{

class Earth_Model
{
  protected:
	std::string name;

	// Density layers
	std::vector<double> density_layer_radii;
	std::vector<std::vector<double>> density_coefficients;

	// Nuclear composition
	std::vector<double> composition_layer_radii;
	std::vector<std::vector<obscura::Element>> elements;
	std::vector<std::vector<double>> nuclear_abundances;

	unsigned int Current_Density_Layer(double r) const;
	unsigned int Current_Composition_Layer(double r) const;

	double Time_of_Layer_Exit(const Event& current_event) const;

  public:
	Earth_Model();

	double Mass_Density(double r) const;

	double Mean_Free_Path(obscura::DM_Particle& DM, double r, double vDM) const;

	obscura::Isotope Sample_Target_Isotope(obscura::DM_Particle& DM, double r, double vDM, std::mt19937& PRNG) const;
	Event Sample_Next_Event(Event& current_event, obscura::DM_Particle& DM, std::mt19937& PRNG) const;

	void Print_Summary(int mpi_rank = 0) const;
};

}	// namespace DaMaSCUS

#endif