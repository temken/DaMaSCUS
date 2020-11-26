#ifndef __Underground_Distribution_hpp_
#define __Underground_Distribution_hpp_

// Headers from libphysica
#include "Numerics.hpp"

// Headers from obscura
#include "DM_Distribution.hpp"

#include "Data_Generation.hpp"

namespace DaMaSCUS
{

class Underground_Distribution : public obscura::DM_Distribution
{
  private:
	libphysica::Interpolation speed_pdf;
	libphysica::Interpolation eta_function;

  public:
	//Constructors
	Underground_Distribution(const Simulation_Data& simulation_data, const Simulation_Data& reference_data, unsigned int isodetection_ring, double halo_density);

	virtual double PDF_Speed(double v) override;
	virtual double Eta_Function(double v_min) override;

	std::string Summary();
};

}	// namespace DaMaSCUS

#endif