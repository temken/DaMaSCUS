#include "Earth_Model.hpp"

#include <cmath>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"

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
	elements				= {
		   {obscura::Get_Element(26), obscura::Get_Element(14), obscura::Get_Element(28), obscura::Get_Element(16), obscura::Get_Element(24), obscura::Get_Element(25), obscura::Get_Element(15), obscura::Get_Element(6), obscura::Get_Element(1)},
		   {obscura::Get_Element(8), obscura::Get_Element(12), obscura::Get_Element(14), obscura::Get_Element(26), obscura::Get_Element(20), obscura::Get_Element(13), obscura::Get_Element(11), obscura::Get_Element(24), obscura::Get_Element(28), obscura::Get_Element(25), obscura::Get_Element(16), obscura::Get_Element(6), obscura::Get_Element(1), obscura::Get_Element(15)},
		   {obscura::Get_Element(8), obscura::Get_Element(14), obscura::Get_Element(13), obscura::Get_Element(26), obscura::Get_Element(20), obscura::Get_Element(19), obscura::Get_Element(11), obscura::Get_Element(12)}};
	nuclear_abundances = {
		{0.855, 0.06, 0.052, 0.019, 0.009, 0.003, 0.002, 0.002, 0.0006},
		{0.44, 0.228, 0.21, 0.0626, 0.0253, 0.0235, 0.0027, 0.0026, 0.002, 0.001, 0.0003, 0.0001, 0.0001, 0.00009},
		{0.466, 0.277, 0.081, 0.05, 0.036, 0.028, 0.026, 0.021}};
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
	double x  = current_event.Radius();
	double v  = current_event.Speed();
	double xv = current_event.position.Dot(current_event.velocity);

	unsigned int current_layer = Current_Density_Layer(x);
	if(current_layer > 9)
	{
		std::cerr << "Error in Earth_Model::Time_of_Layer_Exit(): Event lies outside the Earth with r = " << In_Units(x, km) << " km." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	double R_next = density_layer_radii[current_layer];
	if(current_layer > 0)
	{
		double time_of_min_radial_distance = -1.0 * xv / v / v;
		double min_radial_distance		   = x * x - xv * xv / v / v;
		if(time_of_min_radial_distance > 0 && min_radial_distance < density_layer_radii[current_layer - 1])
			R_next = density_layer_radii[current_layer - 1];
	}

	double t_layer_exit = (-xv + sqrt(R_next * R_next * v * v - v * v * x * x + xv * xv)) / v / v;
	return t_layer_exit + cm / v;
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
	// double xi = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	return obscura::Isotope();
}

Event Earth_Model::Sample_Next_Event(Event& current_event, obscura::DM_Particle& DM, std::mt19937& PRNG) const
{
	return Event();
}

void Earth_Model::Print_Summary(int mpi_rank) const
{
}

}	// namespace DaMaSCUS