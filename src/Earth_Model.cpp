#include "Earth_Model.hpp"

#include <cmath>

// Headers from libphysica
#include "Natural_Units.hpp"

namespace DaMaSCUS
{

using namespace libphysica::natural_units;

Earth_Model::Earth_Model()
: name("Preliminary Reference Earth Model")
{
	density_layer_radii	 = {1221.5 * km, 3480 * km, 5701 * km, 5771 * km, 5971 * km, 6151 * km, 6346.6 * km, 6356 * km, 6368 * km, 6371 * km};
	density_coefficients = {{13.0885, 0.0, -8.8381, 0.0},
							{12.5815, -1.2638, -3.6426, -5.5281},
							{7.9565, -6.4761, 5.5283, -3.0807},
							{5.3197, -1.4836, 0.0, 0.0},
							{11.2494, -8.0298, 0.0, 0.0},
							{11.2494, -3.8045, 0.0, 0.0},
							{2.6910, 0.6924, 0.0, 0.0},
							{2.9, 0.0, 0.0, 0.0},
							{2.6, 0.0, 0.0, 0.0},
							{2.6, 0.0, 0.0, 0.0},
							{0.0, 0.0, 0.0, 0.0}};

	composition_layer_radii = {1221.5 * km, 6368 * km, 6371 * km};
	elements				= {{},
				   {},
				   {}};
	nuclear_abundances		= {{},
						   {},
						   {}};
}

unsigned int Earth_Model::Current_Density_Layer(double r) const
{
	for(unsigned int i = 0; i < density_layer_radii.size(); i++)
		if(r < density_layer_radii[i])
			return i;
	return 10;
}

unsigned int Earth_Model::Current_Composition_Layer(double r) const
{
	for(unsigned int i = 0; i < composition_layer_radii.size(); i++)
		if(r < composition_layer_radii[i])
			return i;
	return 3;
}

double Earth_Model::Time_of_Layer_Exit(const Event& current_event) const
{
	return 0.0;
}

double Earth_Model::Mass_Density(double r) const
{
	int layer  = Current_Density_Layer(r);
	double x   = r / rEarth;
	double rho = 0.0;
	for(unsigned int i = 0; i < 4; i++)
		rho += density_coefficients[layer][i] * pow(x, i);
	return rho * gram / cm / cm / cm;
}

double Earth_Model::Mean_Free_Path(obscura::DM_Particle& DM, double r, double vDM) const
{
	unsigned int composition_layer = Current_Composition_Layer(r);
	double rho					   = Mass_Density(r);
	double lambda_inverse		   = 0.0;
	for(unsigned int i = 0; i < elements[composition_layer].size(); i++)
	{
		double abundance = nuclear_abundances[composition_layer][i];
		for(auto& isotope : elements[composition_layer][i].isotopes)
			lambda_inverse += abundance * isotope.abundance * rho / isotope.mass * DM.Sigma_Nucleus(isotope, vDM);
	}
	return 1.0 / lambda_inverse;
}

obscura::Isotope Earth_Model::Sample_Target_Isotope(obscura::DM_Particle& DM, double r, double vDM, std::mt19937& PRNG) const
{
}

Event Earth_Model::Sample_Next_Event(Event& current_event, obscura::DM_Particle& DM, std::mt19937& PRNG) const
{
}

void Earth_Model::Print_Summary(int mpi_rank) const
{
}

}	// namespace DaMaSCUS