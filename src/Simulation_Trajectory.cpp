#include "Simulation_Trajectory.hpp"

#include <fstream>

// Headers from libphysica
#include "Statistics.hpp"

#include "version.hpp"

namespace DaMaSCUS
{
using namespace libphysica::natural_units;
Trajectory::Trajectory(std::vector<Event>& event_list)
: events(event_list)
{
}

unsigned int Trajectory::Number_of_Scatterings() const
{
	if(events.back().Radius() > rEarth)
		return events.size() - 2;
	else
		return events.size() - 1;
}

double Trajectory::Deepest_Depth() const
{
	double result = rEarth;
	for(unsigned int i = 0; i < events.size() - 1; i++)
	{
		double r_minimum;
		Event event_min = events[i].Point_of_Minimum_Distance();
		if(event_min.time > events[i].time && event_min.time < events[i + 1].time)
			r_minimum = event_min.Radius();
		else
			r_minimum = std::min(events[i].Radius(), events[i + 1].Radius());
		if(r_minimum < result)
			result = r_minimum;
	}
	return rEarth - result;
}

Event Trajectory::Interpolate_Trajectory(double t) const
{
	if(t < events[0].time && t > events.back().time)
	{
		std::cerr << "Error in Trajectory::Interpolate_Trajectory(): Argument t out of bound." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		unsigned int i;
		for(i = 1; i < events.size(); i++)
			if(t < events[i].time)
				break;
		libphysica::Vector vel = events[i - 1].velocity;
		libphysica::Vector pos = events[i - 1].position + (t - events[i - 1].time) * vel;
		return Event(t, pos, vel);
	}
}

std::vector<Event> Trajectory::Points_of_Depth_Crossing(double underground_depth) const
{
	std::vector<Event> result;
	double R = rEarth - underground_depth;
	for(unsigned int i = 0; i < events.size() - 1; i++)
	{
		double r_min = events[i].Point_of_Minimum_Distance().Radius();
		if(r_min < R)
		{
			double x						   = events[i].Radius();
			double v						   = events[i].Speed();
			double xv						   = events[i].position * events[i].velocity;
			double aux						   = sqrt(R * R * v * v - v * v * x * x + xv * xv);
			std::vector<double> crossing_times = {(-xv - aux) / v / v, (-xv + aux) / v / v};
			for(auto& t : crossing_times)
				if(t > 0.0 && t < events[i + 1].time - events[i].time)
					result.push_back(Event(events[i].time + t, events[i].position + t * events[i].velocity, events[i].velocity));
		}
	}
	return result;
}

void Trajectory::Save_to_File(std::string& file_path) const
{
	std::ofstream f;
	f.open(file_path);
	for(auto& event : events)
		f << event.time / sec << "\t"
		  << In_Units(event.position[0], km) << "\t" << In_Units(event.position[1], km) << "\t" << In_Units(event.position[2], km)
		  << In_Units(event.velocity[0], km / sec) << "\t" << In_Units(event.velocity[1], km / sec) << "\t" << In_Units(event.velocity[2], km / sec) << std::endl;
}

void Trajectory::Print_Summary(int mpi_rank) const
{
	if(mpi_rank == 0)
	{
		std::cout << "\nTrajectory summary:" << std::endl
				  << std::endl
				  << "\tInitial speed [km/s]:\t" << In_Units(events[0].Speed(), km / sec) << std::endl
				  << "\tFinal speed [km/s]:\t" << In_Units(events.back().Speed(), km / sec) << std::endl
				  << "\tNumber of scatterings:\t" << Number_of_Scatterings() << std::endl
				  << "\tDeepest depth [km]:\t" << In_Units(Deepest_Depth(), km) << std::endl
				  << SEPARATOR;
	}
}

void Scatter(Event& current_event, obscura::DM_Particle& DM, const obscura::Isotope& target, std::mt19937& PRNG)
{
	libphysica::Vector vel_initial	  = current_event.velocity;
	libphysica::Vector vel_normalized = vel_initial.Normalized();
	double v_initial				  = vel_initial.Norm();

	//Find cos alpha
	// double cosalpha = DM.Sample_Scattering_Angle(xi, target, vDM, IS_Angle);
	double cosalpha = 2.0 * libphysica::Sample_Uniform(PRNG, 0.0, 1.0) - 1.0;

	//Construction of n, the unit vector pointing into the direction of vfinal.
	double cosphi	= cos(libphysica::Sample_Uniform(PRNG, 0.0, 2.0 * M_PI));
	double sinphi	= sqrt(1.0 - cosphi * cosphi);
	double sinalpha = sqrt(1.0 - cosalpha * cosalpha);
	double aux		= sqrt(1.0 - vel_normalized[2] * vel_normalized[2]);
	libphysica::Vector n(
		{cosalpha * vel_normalized[0] + (sinalpha * (-vel_normalized[0] * vel_normalized[2] * cosphi + vel_normalized[1] * sinphi)) / aux,
		 cosalpha * vel_normalized[1] + (sinalpha * (-vel_normalized[1] * vel_normalized[2] * cosphi - vel_normalized[0] * sinphi)) / aux,
		 cosalpha * vel_normalized[2] + aux * cosphi * sinalpha});

	libphysica::Vector vNew = target.mass / (DM.mass + target.mass) * v_initial * n + DM.mass / (DM.mass + target.mass) * vel_initial;
	current_event.velocity	= vNew;
}

Trajectory Simulate_Trajectory(Event initial_conditions, Earth_Model& earth_model, obscura::DM_Particle& DM, std::mt19937& PRNG, double minimal_speed)
{
	std::vector<Event> events;
	events.push_back(initial_conditions);
	Event current_event = initial_conditions;
	// 1. Point of passing into the Earth
	if(initial_conditions.Radius() > rEarth)
	{
		double r	   = current_event.Radius();
		double v	   = current_event.Speed();
		double xv	   = current_event.position * current_event.velocity;
		double t_entry = (-xv - sqrt(rEarth * rEarth * v * v - v * v * r * r + xv * xv)) / v / v;
		current_event.Propagate(t_entry + micro * meter / v);
	}

	// 2. Underground propagation
	double r = current_event.Radius();
	while(r < rEarth)
	{
		current_event = earth_model.Sample_Next_Event(current_event, DM, PRNG);
		r			  = current_event.Radius();
		double v	  = current_event.Speed();
		if(r < rEarth)
		{
			obscura::Isotope target_isotope = earth_model.Sample_Target_Isotope(DM, r, v, PRNG);
			Scatter(current_event, DM, target_isotope, PRNG);
		}
		else
			current_event.Propagate(rEarth / v);

		events.push_back(current_event);
	}
	return Trajectory(events);
}

}	// namespace DaMaSCUS