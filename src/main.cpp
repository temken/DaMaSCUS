#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <fstream>
#include <iostream>
#include <random>

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

// Headers from obscura
#include "Astronomy.hpp"
#include "Configuration.hpp"
#include "Target_Nucleus.hpp"

#include "Earth_Model.hpp"
#include "version.hpp"

using namespace libphysica::natural_units;
using namespace DaMaSCUS;

int main(int argc, char* argv[])
{
	//Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << std::endl;
	////////////////////////////////////////////////////////////////////////
	std::random_device rd;
	std::mt19937 PRNG(rd());
	obscura::Import_Nuclear_Data();
	obscura::Configuration cfg(PROJECT_DIR "bin/config.cfg");
	cfg.Print_Summary();

	Earth_Model earth_model(*cfg.DM, cfg.DM_distr->Maximum_DM_Speed());
	Event test_event(0.0, libphysica::Vector({0.0, 1.0 * meter, 0.0}), libphysica::Vector({0.0, km / sec, 0.0}));

	std::cout << test_event.In_Units(km, sec) << std::endl;
	std::cout << earth_model.Sample_Next_Event(test_event, *cfg.DM, PRNG).In_Units(km, sec) << std::endl;

	Event ic = Initial_Conditions(*cfg.DM_distr, PRNG);
	std::cout << ic.In_Units(km, sec) << "\t" << ic.Radius() / rEarth << std::endl;

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << std::round(1000. * durationTotal) / 1000. << "s";
	if(durationTotal > 60.0)
		std::cout << " (" << floor(durationTotal / 3600.0) << ":" << floor(fmod(durationTotal / 60.0, 60.0)) << ":" << floor(fmod(durationTotal, 60.0)) << ")]." << std::endl;
	else
		std::cout << "]" << std::endl;

	return 0;
}