#include "Underground_Distribution.hpp"

#include <functional>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"
#include "Utilities.hpp"

namespace DaMaSCUS
{

Underground_Distribution::Underground_Distribution(const Simulation_Data& simulation_data, const Simulation_Data& reference_data, unsigned int isodetection_ring, double halo_density)
: DM_Distribution("Underground distribution" + std::to_string(isodetection_ring), 0.0, simulation_data.Lowest_Speed(isodetection_ring), 1.05 * simulation_data.Highest_Speed(isodetection_ring))
{
	DD_use_eta_function						= true;
	speed_pdf								= libphysica::Perform_KDE(simulation_data.data[isodetection_ring], v_domain[0], v_domain[1]);
	DM_density								= simulation_data.Average_Weight(isodetection_ring) / reference_data.Average_Weight(isodetection_ring) * halo_density;
	std::function<double(double)> integrand = [this](double v) {
		return speed_pdf(v) / v;
	};
	std::vector<double> speeds = libphysica::Linear_Space(v_domain[0], v_domain[1], 150);
	std::vector<double> etas;
	for(auto& speed : speeds)
	{
		double eps = libphysica::Find_Epsilon(integrand, speed, v_domain[1], 1.0e-6);
		double eta = libphysica::Integrate(integrand, speed, v_domain[1], eps);
		etas.push_back(eta);
	}
	eta_function = libphysica::Interpolation(speeds, etas);
}

double Underground_Distribution::PDF_Speed(double v)
{
	if(v < v_domain[0] || v > v_domain[1])
		return 0.0;
	else
		return speed_pdf(v);
}

double Underground_Distribution::Eta_Function(double v_min)
{
	if(v_min < v_domain[0] || v_min > v_domain[1])
		return 0.0;
	else
		return eta_function(v_min);
}

std::string Underground_Distribution::Summary()
{
	return "";
}

}	// namespace DaMaSCUS