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
#include "DM_Particle_Standard.hpp"
#include "Target_Nucleus.hpp"

#include "Data_Generation.hpp"
#include "Earth_Model.hpp"
#include "Simulation_Trajectory.hpp"
#include "Underground_Distribution.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;
using namespace DaMaSCUS;

int main(int argc, char* argv[])
{
	// INITIALIZATION
	// 1. MPI initialization
	MPI_Init(NULL, NULL);
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	// 2. Read in configuration file
	obscura::Configuration cfg(PROJECT_DIR "bin/config.cfg");

	// 3. Create a logfile.
	libphysica::Logger log(TOP_LEVEL_DIR "results/" + cfg.ID + "/Logfile.txt");

	// 4. Initial terminal and logfile output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	if(mpi_rank == 0)
		log << "[Started on " << ctime_start << "]" << std::endl
			<< PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			<< DAMASCUS_LOGO
			<< std::endl
			<< "MPI processes:\t" << mpi_processes << std::endl;
	cfg.Print_Summary(mpi_rank);

	//5. Random number generator
	std::random_device rd;
	std::mt19937 PRNG(rd());

	Earth_Model earth_model(*cfg.DM, cfg.DM_distr->Maximum_DM_Speed());
	////////////////////////////////////////////////////////////////////////

	unsigned int number_of_isodetection_rings = 30;
	double underground_depth				  = 1.4 * km;
	unsigned int sample_size				  = 1000;
	double v_min							  = 0.0;

	obscura::DM_Particle_SI reference_DM(cfg.DM->mass, 1.0e-80 * cm * cm);
	earth_model.Interpolate_Mean_Free_Path(reference_DM, cfg.DM_distr->Maximum_DM_Speed());
	Simulation_Data reference_data(sample_size, underground_depth, v_min, number_of_isodetection_rings);
	reference_data.Generate_Data(reference_DM, earth_model, *cfg.DM_distr);
	if(mpi_rank == 0)
		log << reference_data.Summary() << std::endl;

	earth_model.Interpolate_Mean_Free_Path(*cfg.DM, cfg.DM_distr->Maximum_DM_Speed());
	Simulation_Data simulation_data(sample_size, underground_depth, v_min, number_of_isodetection_rings);
	simulation_data.Generate_Data(*cfg.DM, earth_model, *cfg.DM_distr);
	if(mpi_rank == 0)
		log << simulation_data.Summary() << std::endl;

	unsigned int iso_ring = 0;
	Underground_Distribution distr(simulation_data, reference_data, iso_ring, cfg.DM_distr->DM_density);
	std::cout << distr.DM_density * cm * cm * cm << std::endl;
	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	MPI_Barrier(MPI_COMM_WORLD);
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	if(mpi_rank == 0)
		std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]\a" << std::endl;
	MPI_Finalize();
	return 0;
}