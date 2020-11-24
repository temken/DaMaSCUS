#include "Earth_Model.hpp"

#include <cmath>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"
#include "Utilities.hpp"

namespace DaMaSCUS
{

using namespace libphysica::natural_units;

Earth_Model::Earth_Model(obscura::DM_Particle& DM, double vMax)
: name("Preliminary Reference Earth Model")
{
	layer_transitions	 = {1221.5 * km, 3480 * km, 5701 * km, 5771 * km, 5971 * km, 6151 * km, 6346.6 * km, 6356 * km, 6368 * km, 6371 * km};
	density_coefficients = {{13.0885, 0.0, -8.8381, 0.0},
							{12.5815, -1.2638, -3.6426, -5.5281},
							{7.9565, -6.4761, 5.5283, -3.0807},
							{5.3197, -1.4836, 0.0, 0.0},
							{11.2494, -8.0298, 0.0, 0.0},
							{7.1089, -3.8045, 0.0, 0.0},
							{2.6910, 0.6924, 0.0, 0.0},
							{2.9, 0.0, 0.0, 0.0},
							{2.6, 0.0, 0.0, 0.0},
							{2.6, 0.0, 0.0, 0.0},
							{0.0, 0.0, 0.0, 0.0}};

	elements = {
		{obscura::Get_Element(26), obscura::Get_Element(14), obscura::Get_Element(28), obscura::Get_Element(16), obscura::Get_Element(24), obscura::Get_Element(25), obscura::Get_Element(15), obscura::Get_Element(6), obscura::Get_Element(1)},
		{obscura::Get_Element(8), obscura::Get_Element(12), obscura::Get_Element(14), obscura::Get_Element(26), obscura::Get_Element(20), obscura::Get_Element(13), obscura::Get_Element(11), obscura::Get_Element(24), obscura::Get_Element(28), obscura::Get_Element(25), obscura::Get_Element(16), obscura::Get_Element(6), obscura::Get_Element(1), obscura::Get_Element(15)},
		{obscura::Get_Element(8), obscura::Get_Element(14), obscura::Get_Element(13), obscura::Get_Element(26), obscura::Get_Element(20), obscura::Get_Element(19), obscura::Get_Element(11), obscura::Get_Element(12)}};
	nuclear_abundances = {
		{0.855, 0.06, 0.052, 0.019, 0.009, 0.003, 0.002, 0.002, 0.0006},
		{0.44, 0.228, 0.21, 0.0626, 0.0253, 0.0235, 0.0027, 0.0026, 0.002, 0.001, 0.0003, 0.0001, 0.0001, 0.00009},
		{0.466, 0.277, 0.081, 0.05, 0.036, 0.028, 0.026, 0.021}};

	// Interpolate lambda^-1 / rho as a function of the DM speed.
	std::vector<double> speeds = libphysica::Linear_Space(0.0, vMax, 200);
	for(unsigned int i = 0; i < elements.size(); i++)
	{
		std::vector<double> prefactors;
		double r = 0.5 * rEarth;   //cancels
		for(auto& v : speeds)
			prefactors.push_back(1.0 / Mean_Free_Path(DM, r, v) / Mass_Density(r));
		mfp_prefactors.push_back(libphysica::Interpolation(speeds, prefactors));
	}
}

unsigned int Earth_Model::Current_Density_Layer(double r) const
{
	for(unsigned int i = 0; i < layer_transitions.size(); i++)
		if(r < layer_transitions[i])
			return i;
	return 10;
}

unsigned int Earth_Model::Current_Composition_Layer(double r) const
{
	unsigned int density_layer = Current_Density_Layer(r);
	if(density_layer == 0)
		return 0;
	else if(density_layer < 9)
		return 1;
	else if(density_layer == 9)
		return 2;
	else
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

	double R_next		   = layer_transitions[current_layer];
	double time_layer_exit = (-xv + sqrt(R_next * R_next * v * v - v * v * x * x + xv * xv)) / v / v;
	if(current_layer > 0)
	{
		double time_of_min_radial_distance = -1.0 * xv / v / v;
		double min_radial_distance		   = x * x - xv * xv / v / v;
		if(time_of_min_radial_distance > 0 && min_radial_distance < layer_transitions[current_layer - 1])
		{
			R_next			= layer_transitions[current_layer - 1];
			time_layer_exit = (-xv - sqrt(R_next * R_next * v * v - v * v * x * x + xv * xv)) / v / v;
		}
	}

	return time_layer_exit + cm / v;
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

double Earth_Model::Mean_Free_Path(obscura::DM_Particle& DM, double r, double vDM)
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

double Earth_Model::Mean_Free_Path_Interpolated(obscura::DM_Particle& DM, double r, double vDM)
{
	unsigned int composition_layer = Current_Composition_Layer(r);
	double rho					   = Mass_Density(r);
	return 1.0 / (mfp_prefactors[composition_layer](vDM) * rho);
}

obscura::Isotope Earth_Model::Sample_Target_Isotope(obscura::DM_Particle& DM, double r, double vDM, std::mt19937& PRNG) const
{
	// double xi = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	return obscura::Isotope();
}

double Earth_Model::Lambda(Event& event, double distance)
{
	double r		   = event.Radius();
	double v		   = event.Speed();
	double cos_alpha   = event.position.Dot(event.velocity) / r / v;
	unsigned int layer = Current_Density_Layer(r);

	// Parameters and coefficients
	double L_tilde = sqrt(distance * distance + 2.0 * distance * r * cos_alpha + r * r);
	double C_1, C_2, C_3, C_4;
	C_1 = distance;
	C_3 = distance * (r * r + r * distance * cos_alpha + distance * distance / 3.0) / rEarth / rEarth;
	if(std::fabs(cos_alpha + 1) > 1e-14)
	{
		C_2 = 1.0 / 2.0 / rEarth * (L_tilde * (distance + r * cos_alpha) - r * r * cos_alpha + (1.0 - cos_alpha * cos_alpha) * r * r * log((distance + L_tilde + r * cos_alpha) / ((1 + cos_alpha) * r)));
		C_4 = 1.0 / 8 / rEarth / rEarth / rEarth * ((5.0 - 3.0 * cos_alpha * cos_alpha) * (cos_alpha * pow(r, 3) * L_tilde - pow(r, 4) * cos_alpha) + 2.0 * distance * distance * L_tilde * (distance + 3.0 * r * cos_alpha) + distance * L_tilde * r * r * (5.0 + cos_alpha * cos_alpha) + 3.0 * pow(r, 4) * pow(1.0 - cos_alpha * cos_alpha, 2) * log((distance + L_tilde + r * cos_alpha) / ((1.0 + cos_alpha) * r)));
	}
	else
	{
		C_2 = (distance * (distance - 2.0 * r) * sqrt(pow(distance - r, 2))) / (2. * (distance - r)) / rEarth;
		C_4 = (distance * (distance - 2.0 * r) * sqrt(pow(distance - r, 2)) * (pow(distance, 2) - 2.0 * distance * r + 2.0 * pow(r, 2))) / (4. * (distance - r) * rEarth * rEarth * rEarth);
	}
	return mfp_prefactors[Current_Composition_Layer(r)](v) * (density_coefficients[layer][0] * C_1 + density_coefficients[layer][1] * C_2 + density_coefficients[layer][2] * C_3 + density_coefficients[layer][3] * C_4) * gram / cm / cm / cm;
}

Event Earth_Model::Sample_Next_Event(Event& current_event, obscura::DM_Particle& DM, std::mt19937& PRNG)
{
	Event new_event(current_event);
	double lambda	   = 0.0;
	double xi		   = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	double log_xi	   = log(1.0 / xi);
	double r		   = current_event.Radius();
	double v		   = current_event.Speed();
	unsigned int layer = Current_Density_Layer(r);
	while(layer < 10 && lambda < log_xi)
	{
		double time_exit			 = Time_of_Layer_Exit(new_event);
		double max_distance_in_layer = time_exit * v;
		double lambda_i				 = Lambda(new_event, max_distance_in_layer);
		if(lambda + lambda_i < log_xi)
		{
			lambda += lambda_i;
			new_event.Propagate(time_exit);
			layer = Current_Density_Layer(new_event.Radius());
		}
		else
		{
			std::function<double(double)> f = [this, lambda, &new_event, log_xi](double l) {
				return lambda + Lambda(new_event, l) - log_xi;
			};
			double l = libphysica::Find_Root(f, 0.0, max_distance_in_layer, cm);
			new_event.Propagate(l / v);
			break;
		}
	}
	return new_event;
}

void Earth_Model::Print_Summary(int mpi_rank) const
{
}

}	// namespace DaMaSCUS