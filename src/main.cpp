#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <random>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

// Headers from obscura
#include "Astronomy.hpp"
#include "Configuration.hpp"
#include "Target_Nucleus.hpp"

#include "Data_Generation.hpp"
#include "Earth_Model.hpp"
#include "Simulation_Trajectory.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;
using namespace DaMaSCUS;

int main(int argc, char* argv[])
{
	MPI_Init(NULL, NULL);
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	//Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << DAMASCUS_LOGO
			  << std::endl
			  << "MPI processes:\t" << mpi_processes << std::endl;
	////////////////////////////////////////////////////////////////////////
	std::random_device rd;
	std::mt19937 PRNG(rd());
	obscura::Import_Nuclear_Data();
	obscura::Configuration cfg(PROJECT_DIR "bin/config.cfg");
	cfg.Print_Summary();

	Earth_Model earth_model(*cfg.DM, cfg.DM_distr->Maximum_DM_Speed());
	Event test_event(0.0, libphysica::Vector({0.0, 1.0 * meter, 0.0}), libphysica::Vector({0.0, km / sec, 0.0}));

	Event ic = Initial_Conditions(*cfg.DM_distr, 1.5 * rEarth, PRNG);

	Trajectory traj = Simulate_Trajectory(ic, earth_model, *cfg.DM, PRNG);
	traj.Print_Summary();

	std::cout << traj.Points_of_Depth_Crossing(km)[0].Isodetection_Angle(libphysica::Vector({0, 230 * km / sec, 0})) / deg << std::endl;
	std::cout << traj.Points_of_Depth_Crossing(km)[1].Isodetection_Angle(libphysica::Vector({0, 230 * km / sec, 0})) / deg << std::endl;

	Simulation_Data simulation_data(100000, 1.4 * km, 200.0 * km / sec, 1);
	simulation_data.Generate_Data(*cfg.DM, earth_model, *cfg.DM_distr);
	simulation_data.Print_Summary(mpi_rank);
	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	//Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	MPI_Finalize();
	return 0;
}