#include <iostream>
#include <chrono>
#include <cmath>
#include <cstring> // for strlen

#include "version.hpp"
#include "Earth_Model.hpp"

// Headers from libphysica
#include "Natural_Units.hpp"
#include "Utilities.hpp"

// Headers from obscura
#include "Astronomy.hpp"
#include "Configuration.hpp"
#include "Target_Nucleus.hpp"

using namespace libphysica::natural_units;

int main(int argc, char *argv[])
{
	//Initial terminal output
	auto time_start = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto *ctime_start = ctime(&time_start_t);
	if (ctime_start[std::strlen(ctime_start)-1] == '\n') ctime_start[std::strlen(ctime_start)-1] = '\0';
	std::cout 	<<"[Started on " <<ctime_start<<"]" <<std::endl;
	std::cout <<PROJECT_NAME <<"-"<<PROJECT_VERSION <<"\tgit:" <<GIT_BRANCH <<"/" <<GIT_COMMIT_HASH <<std::endl <<std::endl;
	////////////////////////////////////////////////////////////////////////
	
	std::cout <<fib(10) <<std::endl;
	std::cout <<In_Units(1.0, meter/sec)<<std::endl;
	std::cout <<obscura::Fractional_Days_since_J2000(1,1,2001,12)<<std::endl;

	obscura::Import_Nuclear_Data();
	obscura::Configuration cfg(PROJECT_DIR "bin/config.cfg");
	cfg.Print_Summary();

	std::vector<double> DM_masses = libphysica::Log_Space(cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses);
	std::vector<std::vector<double>> exclusion_limits = cfg.DM_detector->Upper_Limit_Curve(*(cfg.DM), *(cfg.DM_distr), DM_masses, cfg.constraints_certainty);
	libphysica::Export_Table(PROJECT_DIR "results/" + cfg.ID + "/constraints.txt", exclusion_limits,{GeV,cm*cm});

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	auto time_end = std::chrono::system_clock::now();
	double durationTotal =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start ).count();
	std::cout 	<<"\n[Finished in "<< std::round(1000.*durationTotal)/1000.<<"s";
	if(durationTotal > 60.0)
		std::cout <<" ("<<floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<")]."<<std::endl;
	else 
		std::cout <<"]"<<std::endl;
	
	return 0;
}