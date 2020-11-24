#include "Simulation_Utilities.hpp"

#include <functional>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Statistics.hpp"

namespace DaMaSCUS
{

using namespace libphysica::natural_units;

// 1. Event Class
Event::Event()
: time(0.0), position({0, 0, 0}), velocity({0, 0, 0})
{
}

Event::Event(double t, const libphysica::Vector& pos, const libphysica::Vector& vel)
: time(t), position(pos), velocity(vel)
{
}

double Event::Radius() const
{
	return position.Norm();
}

double Event::Speed() const
{
	return velocity.Norm();
}

double Event::Angular_Momentum() const
{
	return position.Cross(velocity).Norm();
}

double Event::Isodetection_Angle(const libphysica::Vector& vel_earth) const
{
	return acos(position.Normalized().Dot(vel_earth.Normalized()));
}

int Event::Isodetection_Ring(const libphysica::Vector& vel_earth, unsigned int number_of_rings) const
{
	double isodetection_angle = Isodetection_Angle(vel_earth);
	double ring_width		  = M_PI / number_of_rings;
	return isodetection_angle / ring_width;
}

void Event::Propagate(double dt)
{
	time += dt;
	position += dt * velocity;
}

Event Event::In_Units(double unit_distance, double unit_time) const
{
	return Event(libphysica::natural_units::In_Units(time, unit_time), libphysica::natural_units::In_Units(position, unit_distance), libphysica::natural_units::In_Units(velocity, unit_distance / unit_time));
}

//Overload <<
std::ostream& operator<<(std::ostream& output, const Event& event)
{
	return output << "{"
				  << event.time
				  << ","
				  << event.position
				  << ","
				  << event.velocity
				  << "}";
}

// 2. Generator of initial conditions
Event Initial_Conditions(obscura::DM_Distribution& halo_model, std::mt19937& PRNG)
{
	// 1. Asymptotic initial velocity
	// 1.1. Sample a velocity vector in the galactic rest frame
	libphysica::Vector vel_earth = dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model)->Get_Observer_Velocity();
	dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model)->Set_Observer_Velocity(libphysica::Vector({0, 0, 0}));
	std::function<double(double)> cdf = [&halo_model](double v) {
		return halo_model.CDF_Speed(v);
	};
	double u	 = libphysica::Inverse_Transform_Sampling(cdf, halo_model.Minimum_DM_Speed(), halo_model.Maximum_DM_Speed(), PRNG);
	double phi	 = libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	double theta = acos(libphysica::Sample_Uniform(PRNG, -1.0, 1.0));

	dynamic_cast<obscura::Standard_Halo_Model*>(&halo_model)->Set_Observer_Velocity(vel_earth);
	libphysica::Vector initial_velocity = libphysica::Spherical_Coordinates(u, theta, phi);

	// // 1.2 Boost the vector into the Earth's rest frame.
	initial_velocity -= vel_earth;
	// u = initial_velocity.Norm();

	// 2. Initial position
	// 2.1 Find the maximum impact parameter such that the particle still hits the Sun.
	double impact_parameter_max = rEarth;
	libphysica::Vector e_z		= (-1.0) * initial_velocity.Normalized();
	libphysica::Vector e_x({0, e_z[2], -e_z[1]});
	e_x.Normalize();
	libphysica::Vector e_y = e_z.Cross(e_x);

	// 2.2 Find a random point in the plane..
	double phi_disk						= libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI);
	double xi							= libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	double impact_parameter				= sqrt(xi) * impact_parameter_max;
	libphysica::Vector initial_position = 1.5 * rEarth * e_z + impact_parameter * (cos(phi_disk) * e_x + sin(phi_disk) * e_y);

	return Event(0.0, initial_position, initial_velocity);
}

}	// namespace DaMaSCUS