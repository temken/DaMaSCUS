#ifndef __Earth_Model_hpp_
#define __Earth_Model_hpp_

#include <random>

// Headers from libphysica
#include "Numerics.hpp"

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
	std::vector<double> layer_transitions;
	std::vector<std::vector<double>> density_coefficients;

	// Nuclear composition
	std::vector<std::vector<obscura::Element>> elements;
	std::vector<std::vector<double>> nuclear_abundances;

	// Interpolated speed-dependent prefactor of the mean free path
	std::vector<libphysica::Interpolation> mfp_prefactors;

  public:
	Earth_Model(obscura::DM_Particle& DM, double vMax);
	unsigned int Current_Density_Layer(double r) const;
	unsigned int Current_Composition_Layer(double r) const;

	double Time_of_Layer_Exit(const Event& current_event) const;

	double Mass_Density(double r) const;

	double Mean_Free_Path(obscura::DM_Particle& DM, double r, double vDM);
	double Mean_Free_Path_Interpolated(obscura::DM_Particle& DM, double r, double vDM);

	obscura::Isotope Sample_Target_Isotope(obscura::DM_Particle& DM, double r, double vDM, std::mt19937& PRNG) const;
	double Lambda(Event& event, double distance);
	Event Sample_Next_Event(Event& current_event, obscura::DM_Particle& DM, std::mt19937& PRNG);

	void Print_Summary(int mpi_rank = 0) const;
};

}	// namespace DaMaSCUS

#endif